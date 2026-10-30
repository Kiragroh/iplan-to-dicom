"""Patient-free frozen-runtime check, writing only an explicitly selected report."""
from __future__ import annotations

import io
import json
from pathlib import Path
import tempfile
import tkinter as tk

import numpy as np
import pydicom
from PIL import Image
from scipy.ndimage import map_coordinates
from skimage.measure import find_contours

from . import __version__
from .dicom import DicomContext, dose_dataset, verify_readback, write_dataset
from .gui import ConverterWindow
from .images import decode_png


def run(report_path: Path) -> None:
    report_path = Path(report_path)
    if report_path.exists():
        raise ValueError('Self-test report already exists')
    checks = []
    data = np.arange(16, dtype=np.uint16).reshape(4, 4)
    png = io.BytesIO()
    Image.fromarray(data).save(png, format='PNG')
    decoded = decode_png(png.getvalue(), bits=16, bits_stored=16)
    assert decoded.shape == data.shape
    checks.append('png_decoder')
    mask = np.zeros((10, 10))
    mask[2:7, 2:7] = 1
    assert len(find_contours(mask, .5)) == 1
    checks.append('contour_backend')
    assert np.isfinite(map_coordinates(mask, np.array([[3.5], [3.5]]), order=1)).all()
    checks.append('resampling_backend')
    with tempfile.TemporaryDirectory(prefix='iplan-selftest-') as temp:
        ctx = DicomContext.from_source(pydicom.Dataset())
        ds = dose_dataset(ctx, np.ones((3, 4, 5)), np.eye(4), pydicom.uid.generate_uid())
        path = Path(temp) / 'synthetic-dose.dcm'
        write_dataset(ds, path)
        verify_readback(ds, pydicom.dcmread(path))
    checks.append('dicom_pixels_and_metadata')
    root = tk.Tk()
    root.withdraw()
    window = ConverterWindow(root)
    root.update_idletasks()
    root.update()
    root.destroy()
    checks.append('native_gui')
    report_path.parent.mkdir(parents=True, exist_ok=True)
    with report_path.open('x', encoding='utf-8') as stream:
        json.dump({'version': __version__, 'status': 'passed', 'checks': checks}, stream, indent=2)
