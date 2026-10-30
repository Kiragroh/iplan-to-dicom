# Validation status

This is a prototype for local archival conversion and technical review. Successful
writing, a passing IOD check or a viewer load does not establish clinical
interchangeability with a manufacturer export.

## Automated checks

The synthetic suite covers decoding, corruption rejection, coordinate transforms,
image resampling, contour holes, registration direction, native reference graphs,
patient reconciliation, beams, dose checksums/scaling/support, DICOM references
and serialization, organized batch output, verified resume, GUI selection state,
and safe ZIP extraction. CI runs these without private input data.

Every conversion checks internal references and reads back every written data
element and decoded pixel array. Empty DICOM Type 2 values are compared
semantically. A partial result retains an INCOMPLETE marker and does not become
a completed export.

Independent IOD validation uses `dicom-validator` against an explicitly selected
DICOM edition. The reusable check prints tag numbers and anonymous counts only:

```powershell
python tools/check_iods.py "D:\DICOM-Export\Patient_Example\Plan_Example" --standard-cache "D:\DICOM-spec-cache" --edition 2025e
```

Install `dicom-validator` in a separate validation environment. Downloaded public
standard text is the only network input to that tool. `--sample-series` checks one
image per series plus all RT objects; omit it for the full export.

## Private reference comparisons

Seven paired native cases/manufacturer DICOM exports were examined. The public
reports contain aggregate metrics, no source images, patient names or identifiers.

- Geometry comparison reconstructed 115 ROIs across six cases. One case stops
  because an imported MR's identity conflicts with the selected CT; that check
  is not bypassed.
- A selected CT in the held-out case omits ten intermediate acquired slices,
  creating an approximately 11-mm gap. A source-plane area-only centroid metric
  was misleading for its external contour. The supplemental physical-slab
  calculation is documented separately; the original comparison was retained.
- Native beam geometry/MU comparisons distinguish compatible plan candidates,
  but some current variants share geometry and differ in stored dose. There is
  no automatic choice of the final plan. Native reference RTPLAN lacked MLC
  positions, so those positions have source-level checks but no independent
  manufacturer-DICOM equivalence check.
- Direct reconstruction of stored dose plus source support matched all 1,850,200
  native voxels in the held-out 4.5.6 case within quantization. Another 4.5.6 case
  has one nearest-neighbor boundary discrepancy. Five other cases retain 4.5.5
  stored dose while the reference export reports 4.5.6; their dose differences
  remain unresolved and are not corrected from the reference.
- The actual exported dose uses a different source-derived lattice. Interpolating
  its already masked values onto the manufacturer's lattice introduces additional
  boundary differences, with observed maxima up to about **5.04 Gy** in this
  sample. This is explicitly not a claim of manufacturer-equivalent dose files.

See [geometry results](geometry-validation.json), [geometry interpretation](geometry-findings.md),
[dose results](dose-validation.json) and [dose interpretation](dose-findings.md).
These numerical comparisons are observations, not clinical acceptance criteria.

## Viewer acceptance

QuickLookRT 0.2.38's local WPF viewer test loaded the primary CT volume and linked
RTPLAN from the prototype's flat plan folder. Its automatic isocenter positioning
selected the nearest CT plane. The final export also passed associated dose,
18-ROI overlay, nine-beam, frame/reference and in-memory rendering checks.
The final portable build's checks are recorded in the
[release evidence](release-evidence.md).

The source CT and registered MR contour images live in subfolders. Recursive
discovery is needed to browse those additional series; all referenced objects are
included when importing the complete plan folder into another DICOM system.

## Unverified boundaries

Other iPlan versions, treatment techniques, machines and patient orientations;
vendor-specific ZIP layouts; a clean Windows machine without the development
environment; and live PACS ingestion all need separate testing. Installer files
alone are not compatibility evidence. Export reports and batch failures must be
reviewed before any archive migration is treated as complete.
