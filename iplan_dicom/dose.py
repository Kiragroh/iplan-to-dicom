"""Read stored native beam dose without recalculation or reference DICOM pixels.

AlgoDoseData arrays are little-endian float32 after BLZIP decompression and
byte inversion. The checksum is over these decoded bytes. Samples are dose
per MU; SingleMU already contains component weighting. The native algorithm
grid axes are retained by the low-level reader. The selected TissueProvider
defines their physical grid; read_dose applies the verified grid permutation.
"""
from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path
import re

import numpy as np
from scipy.ndimage import map_coordinates

from .plan import Plan, objects, read_plan, text
from .xbrain import FormatError, numbers, read_payload, resolve_inside, scalar


@dataclass
class StoredDose:
    values_gy: np.ndarray  # [algorithm k, j, i], i is fastest in native bytes
    dimensions_xyz: tuple[int, int, int]
    beam_count: int
    fractions: int
    voxel_to_patient: np.ndarray | None = None
    algorithm_versions: tuple[str, ...] = ()
    sampling: str = "native_algorithm_samples"
    support_applied: bool = False
    dose_deposition_medium: str | None = None

    def require_geometry(self) -> np.ndarray:
        if self.voxel_to_patient is None:
            raise FormatError("Stored dose decoded, but algorithm-grid patient geometry is not verified")
        return self.voxel_to_patient


def decode_dose_samples(payload: bytes, dimensions: tuple[int, int, int],
                        expected_md5: str) -> np.ndarray:
    """Validate decoded bytes, returning the unchanged native algorithm grid."""
    if len(dimensions) != 3 or any(type(n) is not int or n < 1 for n in dimensions):
        raise FormatError("Invalid stored-dose dimensions")
    count = int(np.prod(dimensions, dtype=object))
    if len(payload) != count * 4:
        raise FormatError("Stored-dose byte length disagrees with dimensions")
    if not re.fullmatch(r"[0-9a-fA-F]{32}", expected_md5):
        raise FormatError("Stored-dose checksum is invalid")
    if hashlib.md5(payload).hexdigest() != expected_md5.lower():
        raise FormatError("Stored-dose checksum mismatch")
    array = np.frombuffer(payload, dtype="<f4").reshape(dimensions[::-1])
    if not np.isfinite(array).all() or np.any(array < 0):
        raise FormatError("Stored dose contains negative or non-finite samples")
    return array


def read_stored_dose(plan: Plan, case_dir: Path) -> StoredDose:
    """Sum selected current beam records with their final MU and fractions.

    Alternative ContentVersion schemas and other plan revisions have already
    been resolved by read_plan. All selected beams must have complete matching
    AlgoDoseData grids; missing dose never becomes an implicit zero beam.
    """
    if not plan.beams or type(plan.fractions) is not int or plan.fractions < 1:
        raise FormatError("Complete beam selection and positive fractions are required")
    if plan.source.get("isPrimaryDoseAlgorithmSelected") != "false":
        raise FormatError("Only the verified stored Monte Carlo algorithm is supported")
    tissue = objects(Path(plan.source_directory) / "tissueprovider.xbrain", "TissueProvider")
    if len(tissue) != 1:
        raise FormatError("Selected dose requires one tissue-provider provenance record")
    options = tissue[0].find("XvmcOptions")
    if options is None or options.get("doseResultType") != "DoseToMedium":
        raise FormatError("Stored dose deposition medium is missing or has not been verified")
    total, common_dimensions = None, None
    versions = set()
    for beam in plan.beams:
        node = beam.source.find("AlgoDoseData")
        if node is None:
            raise FormatError("A selected treatment beam has no stored algorithm dose")
        verified_builds = {("iPlan RT Dose 4.5.5", "137"), ("iPlan RT Dose 4.5.6", "145")}
        if (node.get("productName"), node.get("buildNumber")) not in verified_builds:
            raise FormatError("Stored-dose software version has not been verified")
        versions.add(node.get("productName", "unspecified"))
        try:
            dimensions = tuple(int(node.attrib[k]) for k in ("dimX", "dimY", "dimZ"))
        except (KeyError, ValueError):
            raise FormatError("Stored-dose grid dimensions are missing or invalid") from None
        if common_dimensions is not None and dimensions != common_dimensions:
            raise FormatError("Selected beam dose grids have different dimensions")
        source = resolve_inside(case_dir, beam.source_path.parent, text(node, "BinaryFileName"))
        array = decode_dose_samples(read_payload(source), dimensions, text(node, "BinaryFileMD5"))
        if not np.isfinite(beam.mu) or beam.mu <= 0:
            raise FormatError("Selected treatment beam has invalid monitor units")
        if total is None:
            total = np.zeros(array.shape, dtype=np.float64)
            common_dimensions = dimensions
        total += array.astype(np.float64) * beam.mu
    total *= plan.fractions
    if not np.isfinite(total).all():
        raise FormatError("Accumulated dose exceeds numeric range")
    return StoredDose(total, common_dimensions, len(plan.beams), plan.fractions,
                      algorithm_versions=tuple(sorted(versions)),
                      dose_deposition_medium="DOSE_TO_MEDIUM")


def read_dose(plan_dir: Path, case_dir: Path) -> StoredDose:
    """Read a selected plan's physical dose in its planning SliceSet frame.

    Returned values use conventional [z,y,x] tissue-grid order; the affine maps
    their integer voxel centers to patient coordinates. No resampling is needed
    to preserve these stored samples. A later resampling operation must clamp
    to the edge when interpolating inside the native tissue support, as the
    source dose-volume cache does; unrestricted extrapolation is not implied.
    """
    plan = read_plan(plan_dir, case_dir)
    return attach_geometry(read_stored_dose(plan, case_dir), plan_dir)


def tissue_dose_affine(tissue, dimensions: tuple[int, int, int]) -> np.ndarray:
    """Derive conventional X,Y,Z dose centers from native algorithm dimensions.

    X is anchored at the *upper* inclusive tissue crop edge. This matters when
    the cropped source count is not divisible by the coarse-grid stride. XY
    tissue coordinates are pixel corners; Z coordinates are slice centers.
    """
    raw = numbers(tissue.findtext("TiProP3DSubsyst/Point_3D_Subsystem/_DOUBLEVEC_/_ASCII_"),
                  counted=True)
    if len(raw) != 18 or raw[16:] != [2., -1.]:
        raise FormatError("Unsupported tissue-provider coordinate subsystem")
    matrix = np.asarray(raw[:16], float).reshape(4, 4)
    if not np.allclose(matrix[3], [0, 0, 0, 1], atol=1e-10):
        raise FormatError("Tissue-provider geometry is not affine")
    pitch = np.linalg.norm(matrix[:3, :3], axis=0)
    if np.any(pitch <= 0) or not np.isclose(pitch[2], 1, atol=1e-9):
        raise FormatError("Unsupported tissue-provider slice-axis scale")
    unit = matrix[:3, :3] / pitch
    if not np.allclose(unit.T @ unit, np.eye(3), atol=1e-8):
        raise FormatError("Tissue-provider axes are not orthogonal")
    spacing_z = scalar(tissue, "TiProSliceDistanceZ")
    requested = scalar(tissue, "XvmcOptions/SpatialResolution")
    if spacing_z <= 0 or requested != 2.0:
        raise FormatError("Only the verified 2-mm Monte Carlo grid mode is supported")
    if tissue.find("XvmcOptions") is None:
        raise FormatError("Monte Carlo dose geometry metadata is missing")
    matrix[:3, 2] *= spacing_z
    pitch[2] = spacing_z
    stride = np.maximum(1, np.floor(requested / pitch + .5)).astype(int)
    try:
        low = np.array([int(tissue.attrib["tiPro" + axis + "Min"]) for axis in "XYZ"])
        high = np.array([int(tissue.attrib["tiPro" + axis + "Max"]) for axis in "XYZ"])
        size = np.array([int(tissue.attrib[key]) for key in ("tiProWidth", "tiProHeight", "tiProDepth")])
    except (KeyError, ValueError):
        raise FormatError("Tissue-provider crop metadata is incomplete") from None
    if np.any(low < 0) or np.any(high < low) or np.any(high >= size):
        raise FormatError("Tissue-provider crop is out of bounds")
    expected_xyz = (high - low + 1) // stride
    # Native algorithm dimY is tissue Z; dimZ is tissue Y.
    if tuple(expected_xyz[[0, 2, 1]]) != tuple(dimensions):
        raise FormatError("Stored dose dimensions disagree with tissue crop and resolution")
    coarse = np.eye(4)
    coarse[:3, :3] = np.diag([-stride[0], stride[1], stride[2]])
    coarse[:3, 3] = [high[0] + 1 - stride[0] / 2,
                      low[1] + stride[1] / 2, low[2]]
    return np.diag([-1., -1., 1., 1.]) @ matrix @ coarse


def attach_geometry(dose: StoredDose, plan_dir: Path) -> StoredDose:
    tissue = objects(Path(plan_dir) / "tissueprovider.xbrain", "TissueProvider")
    if len(tissue) != 1:
        raise FormatError("Selected plan must have one tissue-provider geometry")
    affine = tissue_dose_affine(tissue[0], dose.dimensions_xyz)
    # [algorithm k,j,i] -> [tissue z,y,x] = [algorithm j,k,i].
    dims = dose.dimensions_xyz
    return StoredDose(dose.values_gy.transpose(1, 0, 2), (dims[0], dims[2], dims[1]),
                      dose.beam_count, dose.fractions, affine, dose.algorithm_versions,
                      dose_deposition_medium=dose.dose_deposition_medium)


@dataclass
class TissueSupport:
    stored_values: np.ndarray  # native uint16 HU+1024; zero denotes empty tissue
    voxel_to_patient: np.ndarray
    crop_min_xyz: tuple[int, int, int] | None = None
    crop_max_xyz: tuple[int, int, int] | None = None


def load_tissue_support(plan_dir: Path, case_dir: Path) -> TissueSupport:
    """Read the selected calculation's checksum-verified fine tissue lattice."""
    found = objects(Path(plan_dir) / "tissueprovider.xbrain", "TissueProvider")
    if len(found) != 1:
        raise FormatError("Selected plan must have one tissue provider")
    tissue = found[0]
    try:
        dims = tuple(int(tissue.attrib[k]) for k in ("tiProWidth", "tiProHeight", "tiProDepth"))
        min_hu, max_hu = int(tissue.attrib["tiProMinHu"]), int(tissue.attrib["tiProMaxHu"])
        low = tuple(int(tissue.attrib["tiPro" + axis + "Min"]) for axis in "XYZ")
        high = tuple(int(tissue.attrib["tiPro" + axis + "Max"]) for axis in "XYZ")
    except (KeyError, ValueError):
        raise FormatError("Tissue-provider dimensions or intensity range are missing") from None
    if min(dims) < 1 or (min_hu, max_hu) != (-1024, 3071):
        raise FormatError("Tissue-provider stored intensity range has not been verified")
    if any(lo < 0 or hi < lo or hi >= dim for lo, hi, dim in zip(low, high, dims)):
        raise FormatError("Tissue-provider calculation crop is invalid")
    source = resolve_inside(case_dir, Path(plan_dir), text(tissue, "BinaryFileName"))
    payload = read_payload(source)
    if len(payload) != int(np.prod(dims, dtype=object))*2:
        raise FormatError("Tissue-provider byte count disagrees with dimensions")
    expected = text(tissue, "BinaryFileMD5").lower()
    if not re.fullmatch(r"[0-9a-f]{32}", expected) or hashlib.md5(payload).hexdigest() != expected:
        raise FormatError("Tissue-provider checksum mismatch")
    data = np.frombuffer(payload, dtype="<u2").reshape(dims[::-1])
    if np.any(data > max_hu-min_hu):
        raise FormatError("Tissue-provider intensities exceed declared range")
    raw = numbers(tissue.findtext("TiProP3DSubsyst/Point_3D_Subsystem/_DOUBLEVEC_/_ASCII_"), counted=True)
    if len(raw) != 18 or raw[16:] != [2., -1.]:
        raise FormatError("Unsupported tissue-provider coordinate subsystem")
    matrix = np.asarray(raw[:16], float).reshape(4, 4)
    distance = scalar(tissue, "TiProSliceDistanceZ")
    if distance <= 0 or not np.allclose(matrix[3], [0, 0, 0, 1], atol=1e-10):
        raise FormatError("Invalid tissue-provider geometry")
    matrix[:3, 2] *= distance
    matrix[:3, 3] += .5*(matrix[:3, 0]+matrix[:3, 1])
    return TissueSupport(data, np.diag([-1., -1., 1., 1.]) @ matrix, low, high)


def sample_tissue_support(support: TissueSupport, points_xyz: np.ndarray) -> np.ndarray:
    """Return source support at patient points using nearest fine tissue voxel."""
    points = np.asarray(points_xyz, dtype=float)
    if points.ndim != 2 or points.shape[0] != 3 or not np.isfinite(points).all():
        raise FormatError("Tissue sampling requires finite 3-by-N patient coordinates")
    inverse = np.linalg.inv(support.voxel_to_patient)
    coordinates = inverse[:3, :3] @ points + inverse[:3, 3, None]
    present = map_coordinates(support.stored_values, coordinates[::-1], order=0,
                              mode="constant", cval=0, prefilter=False) > 0
    if support.crop_min_xyz is not None and support.crop_max_xyz is not None:
        low, high = np.asarray(support.crop_min_xyz)[:, None], np.asarray(support.crop_max_xyz)[:, None]
        present &= np.all((coordinates >= low-.5) & (coordinates < high+.5), axis=0)
    return present


def resample_to_tissue_grid(dose: StoredDose, plan_dir: Path, case_dir: Path) -> StoredDose:
    """Sample first, then mask on the source-derived fine calculation lattice.

    No manufacturer DICOM grid or pixels are used. Clipping the coarse algorithm
    samples *before* interpolation would wrongly remove dose near the body edge.
    The selected tissue buffer provides both the target lattice and its support.
    """
    if dose.support_applied:
        raise FormatError("Dose support has already been applied")
    support = load_tissue_support(plan_dir, case_dir)
    return resample_on_support(dose, support)


def resample_on_support(dose: StoredDose, support: TissueSupport) -> StoredDose:
    """Pure source-grid resampling primitive, also exercised by synthetic tests."""
    transform = np.linalg.inv(dose.require_geometry()) @ support.voxel_to_patient
    zz, yy, xx = support.stored_values.shape
    row, col = np.indices((yy, xx), dtype=float)
    crop_xy = np.ones((yy, xx), dtype=bool)
    z_low, z_high = 0, zz-1
    if support.crop_min_xyz is not None and support.crop_max_xyz is not None:
        low, high = support.crop_min_xyz, support.crop_max_xyz
        crop_xy = (col >= low[0]) & (col <= high[0]) & (row >= low[1]) & (row <= high[1])
        z_low, z_high = low[2], high[2]
    ij = np.vstack((col.ravel(), row.ravel(), np.zeros(xx*yy), np.ones(xx*yy)))
    output = np.zeros((zz, yy, xx), dtype=np.float64)
    coarse_shape_xyz = np.array(dose.values_gy.shape[::-1])[:, None]
    for z in range(zz):
        if z < z_low or z > z_high:
            continue
        ij[2] = z
        coordinates = (transform @ ij)[:3]
        # A cropped fine tissue box can exceed the actual coarse calculation
        # cells after stride truncation. Native caches are zero there even when
        # tissue is present. Clamp only in the half-cell rim, not beyond it.
        cell_domain = np.all((coordinates >= -.5) &
                             (coordinates < coarse_shape_xyz-.5), axis=0).reshape(yy, xx)
        sampled = map_coordinates(dose.values_gy, coordinates[::-1], order=1,
                                   mode="nearest", prefilter=False).reshape(yy, xx)
        output[z] = np.where((support.stored_values[z] > 0) & crop_xy & cell_domain, sampled, 0.)
    return StoredDose(output, (xx, yy, zz), dose.beam_count, dose.fractions,
                      support.voxel_to_patient.copy(), dose.algorithm_versions,
                      "trilinear_on_source_tissue_grid_then_source_support", True,
                      dose.dose_deposition_medium)
