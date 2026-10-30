# Prototype v0.1.0 — release evidence

Verified locally on Windows x64, 2026-10-30:

- **165 synthetic tests passed**, including ZIP traversal/collision/CRC/limit
  checks, batch failure isolation, resume integrity and GUI state handling.
- A private native test case was wrapped in an ordinary ZIP, extracted through
  the public input adapter and converted with the final source code. Result:
  **750 DICOM instances**, including 191 regular viewing CT slices, selected
  source CT, complete original CT, registered/original MR, 18 ROIs, nine beams,
  RTPLAN and RTDOSE. Every instance passed data-element/pixel readback.
- Independent DICOM IOD validation against edition **2025e** checked all
  **750 instances with zero reported errors**; see
  [the machine-readable report](iod-validation.json).
- **QuickLookRT viewer 0.2.38** loaded the final RTPLAN through its actual WPF
  viewer, with a CT volume, 18 ROI overlays, nine beams and one associated dose.
  Reference/frame checks and in-memory rendering passed. The focus was within
  0.363 mm of the isocenter, on the nearest CT plane.
  [Viewer report](quicklook-validation.json).
- The final single-file EXE exited successfully from its patient-free self-test:
  PNG decoding, contour backend, physical resampling, DICOM pixel/metadata
  roundtrip and construction of the native GUI.

The executable bundles Python 3.13 and the
[recorded dependencies](portable-dependencies.txt). It also contains the MIT
license and [third-party notices](third-party-licenses.txt). The release asset's
SHA-256 is provided in `SHA256SUMS.txt`.

Testing was performed on the development workstation, not a clean Windows
installation. The private reference and viewer tests are not reproducible from
the public repository alone because patient data are intentionally excluded.
Native reference dose differences, resampling boundaries, source version limits
and the blocked identity-conflict case remain documented in
[validation status](validation.md). No live PACS import or clinical acceptance
is claimed.
