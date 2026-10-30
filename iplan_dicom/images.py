"""Quantitative XBrain PNG images and checked source-patient geometry.

XBrain stores PNG samples with reversed 16-bit byte order. Its image subsystem
uses XY voxel corners and RAS-like world coordinates; pixel x runs opposite to
the linked original DICOM raster. Native pixels are retained, with an explicit
center affine, rather than changing pixel order implicitly.
"""
from __future__ import annotations

from dataclasses import dataclass
from collections import deque
from copy import deepcopy
import io
from pathlib import Path
from typing import Mapping
import xml.etree.ElementTree as ET

import numpy as np
from PIL import Image
from pydicom.dataset import Dataset

from .xbrain import FormatError, numbers, read_payload, safe_join

RAW_TO_PATIENT = np.diag([-1.0, -1.0, 1.0, 1.0])


@dataclass
class ImageVolume:
    pixels: np.ndarray  # [z,y,x], native pixel order; original stored values
    voxel_to_patient: np.ndarray  # integer voxel centers -> DICOM patient mm
    raw_to_patient: np.ndarray
    headers: tuple[Dataset, ...]
    modality: str
    rescale_slope: float
    rescale_intercept: float
    source_uid: str
    slice_affines: tuple[np.ndarray, ...] = ()
    uniform_grid: bool = True
    derivation_description: str = "Stored XBrain pixel samples reconstructed on the verified source grid"

    @property
    def spacing(self) -> np.ndarray:
        return np.linalg.norm(self.voxel_to_patient[:3, :3], axis=0)

    def slice_affine(self, index: int) -> np.ndarray:
        """Map an individual plane's (x,y,0) pixel centers to patient mm."""
        if self.slice_affines:
            return self.slice_affines[index]
        matrix = self.voxel_to_patient.copy()
        matrix[:3, 3] += index * matrix[:3, 2]
        return matrix


def decode_png(data: bytes, *, bits: int, bits_stored: int | None = None) -> np.ndarray:
    """Decode already-unpacked, non-inverted Brainlab PNG bytes."""
    if not data.startswith(b"\x89PNG\r\n\x1a\n"):
        raise FormatError("Image payload is not PNG")
    try:
        with Image.open(io.BytesIO(data)) as image:
            array = np.asarray(image).copy()
    except Exception:
        raise FormatError("Invalid PNG image payload") from None
    if array.ndim != 2 or bits not in (8, 16):
        raise FormatError("Unsupported image channels or bit depth")
    # The IHDR bit depth is independent of Pillow's in-memory signedness.
    if data[24] != bits:
        raise FormatError("PNG bit depth disagrees with image metadata")
    if bits == 16:
        array = array.astype(np.uint16).byteswap()
    else:
        array = array.astype(np.uint8)
    if bits_stored is not None:
        if not 0 < bits_stored <= bits or np.any(array >= 2**bits_stored):
            raise FormatError("PNG samples exceed source stored-bit range")
    return array


def subsystem_matrix(node: ET.Element) -> np.ndarray:
    values = numbers(node.findtext("SubSystem/Point_3D_Subsystem/_DOUBLEVEC_/_ASCII_"), counted=True)
    if len(values) != 18 or values[16:] != [2.0, -1.0]:
        raise FormatError("Unsupported image subsystem representation")
    matrix = np.array(values[:16], dtype=float).reshape(4, 4)
    if not np.allclose(matrix[3], [0, 0, 0, 1], atol=1e-10):
        raise FormatError("Image subsystem is not affine")
    return matrix


def checked_center_matrix(node: ET.Element, header: Dataset) -> np.ndarray:
    """Cross-check RAS conversion and half pixels against the original header.

    The 0.002 mm bound accommodates six-decimal original MR direction cosines
    over a 512-pixel field; it is not an image-comparison tolerance.
    """
    if node.get("headFeetOrientation") != "HeadFirst" or node.get("patientOrientation") != "Supine":
        raise FormatError("Only verified head-first supine image coordinates are supported")
    if str(getattr(header, "PatientPosition", "")) != "HFS":
        raise FormatError("Original patient position is unsupported")
    matrix = RAW_TO_PATIENT @ subsystem_matrix(node)
    matrix[:3, 3] += 0.5 * (matrix[:3, 0] + matrix[:3, 1])
    try:
        iop = np.asarray(header.ImageOrientationPatient, float)
        ipp = np.asarray(header.ImagePositionPatient, float)
        spacing = np.asarray(header.PixelSpacing, float)
        width, height = int(node.get("width", "0")), int(node.get("height", "0"))
    except (AttributeError, ValueError):
        raise FormatError("Original image geometry is incomplete") from None
    if (int(header.Columns), int(header.Rows)) != (width, height):
        raise FormatError("Image dimensions disagree with original DICOM")
    expected_origin = ipp + iop[:3] * spacing[1] * (width - 1)
    expected_x, expected_y = -iop[:3] * spacing[1], iop[3:] * spacing[0]
    for x, y in ((0, 0), (width - 1, 0), (0, height - 1), (width - 1, height - 1)):
        actual = matrix[:3, 3] + matrix[:3, 0] * x + matrix[:3, 1] * y
        expected = expected_origin + expected_x * x + expected_y * y
        if np.linalg.norm(actual - expected) > 0.002:
            raise FormatError("XBrain and original DICOM geometry disagree")
    return matrix


def validate_source_headers(headers: tuple[Dataset, ...] | list[Dataset], modality: str) -> None:
    """Reject mixed source identities before export replaces original identifiers."""
    if not headers:
        raise FormatError("Image source headers are missing")
    identities = set()
    for header in headers:
        identity = tuple(str(getattr(header, key, "")).strip() for key in
                         ("PatientID", "SeriesInstanceUID", "FrameOfReferenceUID"))
        if not all(identity):
            raise FormatError("Original image patient, series or frame identifier is missing")
        identities.add(identity)
        if str(getattr(header, "Modality", "")) != modality:
            raise FormatError("Original image modality differs from selected image metadata")
        if (str(getattr(header, "PhotometricInterpretation", "")) != "MONOCHROME2"
                or getattr(header, "SamplesPerPixel", None) != 1):
            raise FormatError("Original image photometric interpretation is unsupported")
        if ('LossyImageCompression' in header
                and str(header.LossyImageCompression) not in ('00', '01')):
            raise FormatError("Original image lossy compression history is invalid")
        try:
            allocated, stored, high_bit, representation = (int(getattr(header, key)) for key in
                ("BitsAllocated", "BitsStored", "HighBit", "PixelRepresentation"))
        except (AttributeError, TypeError, ValueError):
            raise FormatError("Original image sample encoding is incomplete") from None
        if (allocated not in (8, 16) or not 0 < stored <= allocated or
                high_bit != stored - 1 or representation not in (0, 1)):
            raise FormatError("Original image sample encoding is inconsistent or unsupported")
    if len(identities) != 1:
        raise FormatError("Image series mixes source patients, series or frames")


def load_image_volume(series: ET.Element, directory: Path, dicom_headers: Mapping[str, Dataset]) -> ImageVolume:
    """Load one selected ImageSeries; source headers are indexed by SOP UID.

    Header identities remain private in memory. This function neither logs nor
    writes them. Only complete, uniformly spaced series are accepted.
    """
    nodes = list(series.findall("Image"))
    if series.tag != "ImageSeries" or len(nodes) < 2:
        raise FormatError("A complete multi-slice ImageSeries is required")
    nodes.sort(key=lambda node: int(node.get("number", "-1")))
    numbers_seen = [int(node.get("number", "-1")) for node in nodes]
    if numbers_seen != list(range(len(nodes))) or int(series.get("imageCount", "-1")) != len(nodes):
        raise FormatError("Image numbering or count is inconsistent")
    arrays, matrices, headers = [], [], []
    metadata = set()
    for node in nodes:
        header = dicom_headers.get(node.get("dicomInstanceUID", ""))
        if header is None:
            raise FormatError("Original DICOM header is unavailable")
        if node.get("type") != "png" or node.get("modality") not in ("CT", "MR"):
            raise FormatError("Unsupported image format or modality")
        validate_source_headers([header], node.get("modality"))
        representation = int(getattr(header, "PixelRepresentation", -1))
        if representation not in (0, 1):
            raise FormatError("Original image signedness is missing")
        bits = int(node.get("bits", "0"))
        if bits != int(header.BitsAllocated):
            raise FormatError("PNG bit depth differs from original allocated sample width")
        array = decode_png(read_payload(safe_join(directory, node.get("file", "")), invert=False),
                           bits=bits, bits_stored=int(header.BitsStored))
        if representation == 1:
            if int(header.BitsStored) != bits:
                raise FormatError("Signed samples with unused high bits are not yet supported")
            array = array.view(np.int16 if bits == 16 else np.int8)
        if array.shape != (int(header.Rows), int(header.Columns)):
            raise FormatError("PNG dimensions disagree with source header")
        matrices.append(checked_center_matrix(node, header))
        arrays.append(array)
        headers.append(header)
        slope = float(getattr(header, "RescaleSlope", 1))
        intercept = float(getattr(header, "RescaleIntercept", 0))
        if not np.isfinite(slope) or not np.isfinite(intercept) or slope == 0:
            raise FormatError("Image intensity calibration is invalid")
        if node.get("modality") == "CT" and not all(hasattr(header, k) for k in ("RescaleSlope", "RescaleIntercept")):
            raise FormatError("CT HU calibration is missing")
        metadata.add((node.get("modality"), slope, intercept, bits, representation))
    if len(metadata) != 1:
        raise FormatError("Series has varying modality or intensity calibration")
    origins = np.array([matrix[:3, 3] for matrix in matrices])
    # Individual header positions are rounded. A median of rounded adjacent
    # differences accumulates directional rounding bias across a long series.
    step = (origins[-1] - origins[0]) / (len(origins) - 1)
    if np.linalg.norm(step) <= 1e-6:
        raise FormatError("Image slices have zero spacing")
    affine = matrices[0].copy()
    affine[:3, 2] = step
    for index, matrix in enumerate(matrices):
        if not np.allclose(matrix[:3, :2], affine[:3, :2], atol=1e-6, rtol=0):
            raise FormatError("Image series orientation varies")
        if np.linalg.norm(origins[index] - origins[0] - index * step) > 0.002:
            raise FormatError("Image series is not uniformly spaced")
    modality, slope, intercept, _, _ = metadata.pop()
    validate_source_headers(headers, modality)
    return ImageVolume(np.stack(arrays), affine, RAW_TO_PATIENT.copy(), tuple(headers), modality,
                       slope, intercept, series.get("uId", ""), tuple(matrices))


def load_slice_set(slice_set: ET.Element, directory: Path, case_dir: Path,
                   dicom_headers: Mapping[str, Dataset]) -> ImageVolume:
    """Reconstruct the plan's selected and individually registered image planes.

    This preserves the selected SliceSet geometry. Its planes can differ from
    the original ImageSeries and need not form one regular 3D affine grid.
    Writers must use ``slice_affine(index)`` for every output image.
    """
    from .plan import referenced_object, resolve_directory

    slices = slice_set.find("Slices")
    if slices is None:
        raise FormatError("Selected SliceSet has no planes")
    series = referenced_object(slices, "ImageSeries", directory, case_dir)
    ref = next(r for r in slices.findall("Reference") if r.get("contentType") == "ImageSeries")
    series_directory = resolve_directory(case_dir, directory, ref.get("directory", "."))
    original = load_image_volume(series, series_directory, dicom_headers)
    if (int(slice_set.get("width", "0")), int(slice_set.get("height", "0"))) != (
        original.pixels.shape[2], original.pixels.shape[1]
    ):
        raise FormatError("Selected SliceSet changes image dimensions")
    nodes = list(slices.findall("Slice"))
    if len(nodes) != int(slices.get("sliceCount", "-1")) or len(nodes) < 2:
        raise FormatError("Selected SliceSet plane count is inconsistent")
    try:
        indices = [int(node.attrib["imageIndex"]) for node in nodes]
    except (KeyError, ValueError):
        raise FormatError("Selected SliceSet image indices are invalid") from None
    if len(set(indices)) != len(indices) or min(indices) < 0 or max(indices) >= len(original.pixels):
        raise FormatError("Selected SliceSet image indices are out of bounds or repeated")
    matrices = []
    for node in nodes:
        matrix = RAW_TO_PATIENT @ subsystem_matrix(node)
        matrix[:3, 3] += .5 * (matrix[:3, 0] + matrix[:3, 1])
        norms = np.linalg.norm(matrix[:3, :2], axis=0)
        expected_norms = original.spacing[:2]
        if not np.allclose(norms, expected_norms, atol=1e-6, rtol=0):
            raise FormatError("Selected SliceSet changes pixel scale")
        if abs(np.dot(matrix[:3, 0], matrix[:3, 1])) > 1e-6:
            raise FormatError("Selected SliceSet plane axes are not orthogonal")
        matrices.append(matrix)
    origins = np.array([matrix[:3, 3] for matrix in matrices])
    affine = matrices[0].copy()
    affine[:3, 2] = (origins[-1] - origins[0]) / (len(matrices) - 1)
    uniform = all(np.allclose(matrix[:3, :2], affine[:3, :2], atol=1e-6, rtol=0)
                  and np.linalg.norm(matrix[:3, 3] - affine[:3, 3] - i * affine[:3, 2]) <= .002
                  for i, matrix in enumerate(matrices))
    return ImageVolume(original.pixels[indices], affine, original.raw_to_patient,
                       tuple(original.headers[i] for i in indices), original.modality,
                       original.rescale_slope, original.rescale_intercept, slice_set.get("uId", ""),
                       tuple(matrices), uniform)


def registration_matrix(source_uid: str, target_uid: str, plan_dir: Path) -> np.ndarray:
    """Read the selected plan's direct or uniquely shortest native registration.

    Native Fusion references are ordered source, target. A direct registration
    has precedence over a longer composed path: independently estimated native
    registrations are not necessarily cycle-consistent. Equal shortest paths
    must agree numerically, otherwise the registration is ambiguous.
    """
    from .plan import objects, vector

    if not source_uid or not target_uid:
        raise FormatError("Registration SliceSet identity is missing")
    if source_uid == target_uid:
        return np.eye(4)
    graph: dict[str, list[tuple[str, np.ndarray]]] = {}
    for fusion in objects(Path(plan_dir) / "fusion.xbrain", "Fusion"):
        refs = fusion.findall("Reference")
        if len(refs) != 2 or any(r.get("contentType") != "SliceSet" or not r.get("contentUid") for r in refs):
            raise FormatError("Fusion must reference exactly two image SliceSets")
        values = vector(fusion, "Transformation/Matrix_44")
        if values.size != 16:
            raise FormatError("Fusion matrix has invalid dimensions")
        matrix = values.reshape(4, 4)
        rotation = matrix[:3, :3]
        if (not np.isfinite(matrix).all() or not np.allclose(matrix[3], [0, 0, 0, 1], atol=1e-10)
                or not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-6, rtol=0)
                or abs(np.linalg.det(rotation) - 1) > 1e-6):
            raise FormatError("Only rigid native image fusions are supported")
        source, target = refs[0].get("contentUid"), refs[1].get("contentUid")
        for start, end, transform in ((source, target, matrix), (target, source, np.linalg.inv(matrix))):
            edges = graph.setdefault(start, [])
            if not any(end == other and np.allclose(transform, prior, atol=1e-8, rtol=0) for other, prior in edges):
                edges.append((end, transform))
    queue = deque([(source_uid, np.eye(4), (source_uid,))])
    shortest = None
    solutions = []
    visited_paths = 0
    while queue:
        current, transform, path = queue.popleft()
        hops = len(path) - 1
        if shortest is not None and hops > shortest:
            break
        if current == target_uid:
            shortest = hops
            solutions.append(transform)
            continue
        for target, edge in graph.get(current, []):
            if target not in path:
                queue.append((target, edge @ transform, path + (target,)))
                visited_paths += 1
                if visited_paths > 10000:
                    raise FormatError("Native registration graph is excessively ambiguous")
    if not solutions:
        raise FormatError("No native fusion links the ROI image grid to the plan grid")
    if any(not np.allclose(candidate, solutions[0], atol=1e-6, rtol=0) for candidate in solutions[1:]):
        raise FormatError("Native fusion has conflicting shortest registration paths")
    return solutions[0]


def register_slice_set(image: ImageVolume, source_uid: str, target_uid: str, plan_dir: Path) -> ImageVolume:
    """Place an owning ROI grid in the selected CT planning frame.

    Pixel samples stay unchanged. Registration affects every physical plane,
    including contour geometry; original header identities are retained only
    for provenance in memory. The plan's raw point transform remains unchanged.
    """
    if image.source_uid != source_uid:
        raise FormatError("Registration source differs from the loaded SliceSet")
    transform = registration_matrix(source_uid, target_uid, plan_dir)
    patient_transform = RAW_TO_PATIENT @ transform @ RAW_TO_PATIENT
    matrices = tuple(patient_transform @ image.slice_affine(index) for index in range(len(image.pixels)))
    return ImageVolume(image.pixels, patient_transform @ image.voxel_to_patient,
                       image.raw_to_patient.copy(), image.headers, image.modality,
                       image.rescale_slope, image.rescale_intercept, image.source_uid,
                       matrices, image.uniform_grid,
                       image.derivation_description if source_uid == target_uid else
                       "Stored XBrain pixel samples registered into the selected planning frame using its native rigid fusion")


def tissue_image_affine(plan_dir: Path) -> np.ndarray:
    """Return the native tissue-provider's regular physical voxel-center grid.

    Only geometry is used. Tissue-density buffer values are not substituted
    for diagnostic source CT intensities.
    """
    from .plan import objects, vector
    from .xbrain import scalar

    candidates = objects(Path(plan_dir) / "tissueprovider.xbrain", "TissueProvider")
    if len(candidates) != 1:
        raise FormatError("Regular CT requires one native tissue-provider grid")
    tissue = candidates[0]
    values = vector(tissue, "TiProP3DSubsyst/Point_3D_Subsystem")
    if len(values) != 18 or values[16:].tolist() != [2., -1.]:
        raise FormatError("Unsupported regular CT coordinate subsystem")
    matrix = values[:16].reshape(4, 4).copy()
    step = scalar(tissue, "TiProSliceDistanceZ")
    if step <= 0 or not np.isclose(np.linalg.norm(matrix[:3, 2]), 1, atol=1e-8):
        raise FormatError("Invalid regular CT slice spacing")
    matrix[:3, 2] *= step
    matrix[:3, 3] += .5 * (matrix[:3, 0] + matrix[:3, 1])
    return RAW_TO_PATIENT @ matrix


def regularize_image_volume(image: ImageVolume, target_affine: np.ndarray) -> ImageVolume:
    """Resample source CT onto a regular grid covering its complete field.

    Each target-grid normal line intersects the actual registered source
    planes. Bilinear values on those planes are interpolated linearly between
    adjacent physical intersections. This supports small independently varying
    plane rotations/translations without replacing them by an approximate
    volume affine. Crossing planes are rejected. Exterior coverage is air.
    """
    from scipy.ndimage import map_coordinates

    if image.modality != "CT" or len(image.pixels) < 2:
        raise FormatError("Regular image reconstruction currently supports multi-plane CT only")
    target = np.asarray(target_affine, float).copy()
    if target.shape != (4, 4) or not np.isfinite(target).all() or not np.allclose(target[3], [0, 0, 0, 1]):
        raise FormatError("Regular CT grid is not a finite affine")
    spacing = np.linalg.norm(target[:3, :3], axis=0)
    if np.any(spacing <= 0) or not np.allclose((target[:3, :3] / spacing).T @ (target[:3, :3] / spacing), np.eye(3), atol=1e-8):
        raise FormatError("Regular CT grid must have orthogonal axes")
    if image.rescale_slope == 0 or not np.isfinite([image.rescale_slope, image.rescale_intercept]).all():
        raise FormatError("Regular CT intensity calibration is invalid")
    height, width = image.pixels.shape[1:]
    corners = np.array([[0, 0, 0, 1], [width - 1, 0, 0, 1],
                        [0, height - 1, 0, 1], [width - 1, height - 1, 0, 1]], float).T
    source_affines = [image.slice_affine(index) for index in range(len(image.pixels))]
    inverse = np.linalg.inv(target)
    all_corners = np.concatenate([(inverse @ affine @ corners)[:3] for affine in source_affines], axis=1)
    # Only a floating-point roundoff margin is used when snapping exact integers.
    lower = np.floor(all_corners.min(axis=1) + 1e-8).astype(int)
    upper = np.ceil(all_corners.max(axis=1) - 1e-8).astype(int)
    dimensions = upper - lower + 1
    if np.any(dimensions <= 0) or np.prod(dimensions, dtype=np.int64) > 300_000_000:
        raise FormatError("Regular CT output dimensions exceed the supported limit")
    target[:3, 3] += target[:3, :3] @ lower
    out_width, out_height, out_depth = map(int, dimensions)
    yy, xx = np.indices((out_height, out_width), dtype=float)
    mappings = [np.linalg.inv(affine) @ target for affine in source_affines]
    coefficients = []
    for mapping in mappings:
        if abs(mapping[2, 2]) < 1e-6:
            raise FormatError("Source CT plane is parallel to the target sampling line")
        coefficients.append(-mapping[2, [0, 1, 3]] / mapping[2, 2])
    centers = [c[0] * (out_width - 1) / 2 + c[1] * (out_height - 1) / 2 + c[2] for c in coefficients]
    ordering = np.argsort(centers)
    limits = np.iinfo(image.pixels.dtype)
    padding = (-1024.0 - image.rescale_intercept) / image.rescale_slope
    if not limits.min <= padding <= limits.max:
        raise FormatError("Air HU cannot be represented in the source CT sample encoding")
    output = np.full((out_depth, out_height, out_width), round(padding), dtype=image.pixels.dtype)

    def projected(index):
        mapping, coeff = mappings[index], coefficients[index]
        intersections = coeff[0] * xx + coeff[1] * yy + coeff[2]
        x = mapping[0, 0] * xx + mapping[0, 1] * yy + mapping[0, 2] * intersections + mapping[0, 3]
        y = mapping[1, 0] * xx + mapping[1, 1] * yy + mapping[1, 2] * intersections + mapping[1, 3]
        # Avoid treating floating-point roundoff at an exact image edge as air.
        x[np.abs(x) < 1e-8], y[np.abs(y) < 1e-8] = 0, 0
        x[np.abs(x - (width - 1)) < 1e-8] = width - 1
        y[np.abs(y - (height - 1)) < 1e-8] = height - 1
        values = map_coordinates(image.pixels[index].astype(np.float64), [y, x], order=1, mode="constant", cval=padding)
        return intersections, values

    previous_t, previous_values = projected(ordering[0])
    for current_index in ordering[1:]:
        current_t, current_values = projected(current_index)
        gap = current_t - previous_t
        if np.any(gap <= 1e-6):
            raise FormatError("Registered source CT planes cross within the output field")
        first = max(0, int(np.ceil(previous_t.min() - 1e-8)))
        last = min(out_depth - 1, int(np.floor(current_t.max() + 1e-8)))
        for z in range(first, last + 1):
            covered = (z >= previous_t - 1e-8) & (z <= current_t + 1e-8)
            weight = np.clip((z - previous_t[covered]) / gap[covered], 0, 1)
            values = previous_values[covered] * (1 - weight) + current_values[covered] * weight
            output[z, covered] = np.clip(np.rint(values), limits.min, limits.max).astype(output.dtype)
        previous_t, previous_values = current_t, current_values
    headers = []
    center_values = np.asarray(centers)
    for z in range(out_depth):
        header = deepcopy(image.headers[int(np.argmin(abs(center_values - z)))])
        header.SliceThickness = float(spacing[2])
        headers.append(header)
    return ImageVolume(output, target, image.raw_to_patient.copy(), tuple(headers), image.modality,
                       image.rescale_slope, image.rescale_intercept, image.source_uid, (), True,
                       "Original quantitative CT resampled onto the native tissue-provider regular grid; bilinear plane sampling and physical linear inter-plane interpolation; exterior air")
