"""Strict decoding of the observed XBrain planar run-length mask format."""
from __future__ import annotations

from dataclasses import dataclass
import struct
from typing import Iterator

import numpy as np
from skimage.measure import find_contours

from .images import ImageVolume
from .xbrain import FormatError


@dataclass
class MaskVolume:
    runs: tuple[np.ndarray, ...]
    width: int
    height: int
    xy_scale: float
    transform: np.ndarray

    @property
    def shape(self) -> tuple[int, int, int]:
        return len(self.runs), self.height, self.width

    def plane(self, index: int) -> np.ndarray:
        runs = self.runs[index]
        return np.repeat(np.arange(len(runs), dtype=np.uint32) % 2, runs.astype(np.int64)).astype(bool).reshape(self.height, self.width)

    def iter_planes(self) -> Iterator[np.ndarray]:
        for index in range(len(self.runs)):
            yield self.plane(index)

    @property
    def voxel_count(self) -> int:
        return sum(int(runs[1::2].sum(dtype=np.uint64)) for runs in self.runs)


def decode_voxel_mask(data: bytes, *, base_width: int, base_height: int, xy_scale: float) -> MaskVolume:
    """Parse decompressed and XOR-inverted mask bytes, retaining sparse planes.

    Header: integer width byte 4, real width byte 8, uint32 matrix count16,
    16 float64 affine values, uint32 slice count. Each slice is uint32 run count
    and uint32 run lengths alternating background and foreground.
    """
    if len(data) < 138 or data[:2] != b"\x04\x08":
        raise FormatError("Unsupported or truncated voxel mask header")
    if struct.unpack_from("<I", data, 2)[0] != 16:
        raise FormatError("Voxel mask matrix size is unsupported")
    transform = np.frombuffer(data, dtype="<f8", offset=6, count=16).reshape(4, 4).copy()
    if not np.isfinite(transform).all() or not np.allclose(transform, np.eye(4), atol=1e-12):
        raise FormatError("Non-identity mask transforms are not yet verified")
    if not np.isfinite(xy_scale) or xy_scale <= 0 or base_width <= 0 or base_height <= 0:
        raise FormatError("Invalid mask dimensions or scaling")
    width, height = int(base_width * xy_scale), int(base_height * xy_scale)
    if width != base_width * xy_scale or height != base_height * xy_scale or width * height > 64_000_000:
        raise FormatError("Unsupported mask dimensions")
    count = struct.unpack_from("<I", data, 134)[0]
    if not 0 < count <= 10000:
        raise FormatError("Invalid mask slice count")
    position, planes = 138, []
    for _ in range(count):
        if position + 4 > len(data):
            raise FormatError("Truncated mask plane header")
        run_count = struct.unpack_from("<I", data, position)[0]
        position += 4
        if not 0 < run_count <= width * height + 1 or position + run_count * 4 > len(data):
            raise FormatError("Invalid or truncated mask run data")
        runs = np.frombuffer(data, dtype="<u4", offset=position, count=run_count).copy()
        position += run_count * 4
        if int(runs.sum(dtype=np.uint64)) != width * height:
            raise FormatError("Mask run lengths disagree with dimensions")
        planes.append(runs)
    if position != len(data):
        raise FormatError("Unexpected trailing mask data")
    return MaskVolume(tuple(planes), width, height, xy_scale, transform)


def mask_to_contours(mask: MaskVolume, image: ImageVolume) -> list[tuple[int, np.ndarray]]:
    """Return (slice index, Nx3 patient-mm polygon), using CLOSEDPLANAR_XOR.

    Padding preserves objects touching a field edge. The half-pixel mapping
    preserves the same physical corners at finer XY mask resolution. All rings,
    including holes, must be written with DICOM CLOSEDPLANAR_XOR semantics.
    """
    if mask.shape[0] != image.pixels.shape[0] or (mask.width, mask.height) != (
        image.pixels.shape[2] * mask.xy_scale, image.pixels.shape[1] * mask.xy_scale
    ):
        raise FormatError("Mask and selected image grids do not match")
    if not np.allclose(mask.transform, np.eye(4), atol=1e-12):
        raise FormatError("Non-identity mask contour mapping is unsupported")
    contours = []
    for z, plane in enumerate(mask.iter_planes()):
        if not plane.any():
            continue
        for ring in find_contours(np.pad(plane.astype(np.uint8), 1), 0.5, fully_connected="low"):
            ring = ring - 1
            if len(ring) < 4 or not np.allclose(ring[0], ring[-1]):
                raise FormatError("Mask contour is not closed")
            ring = ring[:-1]  # CLOSEDPLANAR_XOR implies closure; avoid duplicate vertex
            voxels = np.column_stack(((ring[:, 1] + .5) / mask.xy_scale - .5,
                                      (ring[:, 0] + .5) / mask.xy_scale - .5,
                                      np.zeros(len(ring)), np.ones(len(ring))))
            contours.append((z, (image.slice_affine(z) @ voxels.T).T[:, :3]))
    return contours
