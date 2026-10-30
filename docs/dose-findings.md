# Stored iPlan dose: anonymous technical findings

Investigation date: 2026-10-30. Sources were read only. Development cases 1–6
were inspected. Reserve case 7 was evaluated only after the affine formula was
fixed from case 1, and confirmed it without any formula changes. The indices are transient sorted
case/plan positions, not patient identifiers. No source filenames, identities,
UIDs, patient coordinates, or dose pixels are included here.

## Confirmed storage and summation

- Dose is present even when a case has no `.dosevol` files. All six development
  cases contain `AlgoDoseData` records in native conformal beam/arc objects.
- `AlgoDoseData` provides `dimX`, `dimY`, `dimZ`, `BinaryFileName`, and
  `BinaryFileMD5`. The binary reference is relative to the **referring file's
  directory**. Identical binary filenames occur in different plan directories;
  a case-wide basename lookup selected wrong data in two case-2 beams.
- Decode BLZIP (5-byte marker, big-endian uncompressed length, zlib from byte 9),
  then invert every decoded byte with XOR 255. Interpret as little-endian
  IEEE-754 float32. Byte length must be exactly `4*dimX*dimY*dimZ`.
- `BinaryFileMD5` matches the decompressed, byte-inverted payload, not the
  compressed file or the un-inverted decompressed payload. This checksum was
  verified for every beam in the successfully decoded plans below.
- All verified arrays are finite and nonnegative. Native `dimX` is the fastest
  varying array index; retaining native order gives shape `(dimZ,dimY,dimX)`.
  These algorithm axes must not be mistaken for patient X/Y/Z.
- Select the highest supported `ContentVersion` serialization schema, not all
  schema alternatives. For example, case 2 contains duplicate schema-4 and
  schema-3 serializations of the same four arcs. Historical revisions are
  resolved separately through the unversioned `Revision` pointer.
- Selected plan dose is consistent with `sum(beam_array * SingleMU) * fractions`.
  `SingleMU` already contains component weighting; multiplying component weights
  again would be incorrect. Fractions come from the selected coordinator's
  fractionated target objects. Missing beam dose is an error, not a zero beam.
- Case 2 has four arcs and five static beams. With five fractions, stored-grid
  maximum is 31.242069 Gy; its native RTDOSE maximum is 31.229500 Gy. The different
  grids prevent a direct equality claim. This supports the physical scaling;
  it does not alone establish patient geometry.

## Checked development plans

All rows below passed decoded length, MD5, finite/nonnegative, uniform dimension,
and complete beam summation checks. Plan indices are zero-based current-directory
ordering and are not a recommendation of the clinically intended revision.

| Case | Plan | Beams | Fractions | Algorithm dimensions X,Y,Z | Stored-grid max Gy |
|---:|---:|---:|---:|---|---:|
| 1 | 0 | 8 | 30 | 111,93,111 | 60.022747 |
| 1 | 2 | 9 | 30 | 111,93,111 | 59.900220 |
| 1 | 3 | 9 | 30 | 111,93,111 | 60.237674 |
| 2 | 0 | 9 | 5 | 126,88,125 | 31.242069 |
| 3 | 0 | 11 | 5 | 92,81,90 | 31.352305 |
| 3 | 1 | 10 | 5 | 92,81,90 | 31.235741 |
| 3 | 2 | 10 | 5 | 92,81,90 | 31.264859 |
| 4 | 0 | 8 | 5 | 100,95,103 | 31.335579 |
| 5 | 0 | 8 | 5 | 116,100,116 | 31.290159 |
| 5 | 1 | 8 | 5 | 116,100,116 | 31.278589 |

Other current directories include setup/incomplete plans without stored dose or
unambiguous fractionation. The final per-case comparison below includes all seven
cases, using plan-scoped machine definitions; unsupported variants fail explicitly.

## Separate dose-volume caches

Case 1 has two `.dosevol` / `.dosevolbin` pairs. Each binary is
130*129*95*4 = 6,372,600 decoded bytes, with the same inversion and float32
encoding. Their individual maxima are 56.929535 and 5.112919 Gy; their voxelwise
sum maximum is 59.171554 Gy. Metadata describes group-level Monte Carlo volumes
and invalid fraction counters (-1), not a self-contained complete selected plan.
Their maximum differs from the current per-beam plans and native RTDOSE
(60.043831 Gy). Therefore these caches cannot silently replace the selected
plan's complete beam-dose graph. Their metadata uniquely associates both caches with current plan 3 by PlanName,
RTPlanName and component names, with matching TissueProviderMd5. They are
resampled group-dose caches rather than the native algorithm lattice; therefore
their lower sampled maximum does not imply a different physical dose.

## Verified source geometry

The selected `TissueProvider` contains an orthogonal fitted image grid in the
selected planning SliceSet's frame. Original CT DICOM coordinates cannot
substitute for that planning frame. Relevant fields are
`TiProP3DSubsyst/Point_3D_Subsystem`, `TiProSliceDistanceZ`, crop limits
`tiProXMin/Max`, `tiProYMin/Max`, `tiProZMin/Max`, and
`XvmcOptions/SpatialResolution`.

For the observed 2-mm Monte Carlo mode, reduction factors are the nearest positive
integer to requested physical resolution divided by tissue-grid pitch. The
TissueProvider matrix's third column is a unit direction and is multiplied by
`TiProSliceDistanceZ`. Expected algorithm dimensions are floor(crop voxel count /
stride), permuted from tissue X,Y,Z to algorithm X,Z,Y. This relation is checked
before export. XY tissue subsystem positions are pixel corners; Z positions are
slice centers.

The source-defined mapping is:

```
TP_x = xmax + 1 - sx/2 - sx*i
TP_y = ymin + sy/2 + sy*k
TP_z = zmin + sz*j
```

Here `(i,j,k)` are stored algorithm indices and `sx/sy/sz` the reduction factors.
Reversed X is anchored at the **upper inclusive crop edge**, not the lower edge.
These differ when crop length is not divisible by stride; this explains the
initial case-1 one-pixel offset. Apply the TissueProvider matrix and the established
raw-to-patient sign conversion `diag(-1,-1,+1,1)` afterward. The returned dose
volume is transposed to conventional tissue `[z,y,x]` order without resampling.

Independent source cache confirmation (case 1, plan 3):

- Static group: summation times final MU and fractions, sampled on the explicit
  `.dosevol` affine, reproduces 5,340 cache samples above 1 Gy with RMSE
  7.5e-8 Gy and maximum difference 4.1e-7 Gy.
- Arc group: 432,683 samples above 0.1 Gy match with RMSE 1.59e-7 Gy and maximum
  difference 6.01e-6 Gy. Trilinear interpolation must clamp at the edge of the
  algorithm center lattice within the native tissue support. Using zero outside
  the center lattice caused the initially unexplained boundary residual.
- A diagnostic least-squares check recovered all seven stored arc MU coefficients
  within 1.3e-7 MU. No fitted coefficients, offsets, or reference pixels are used
  in the final reconstruction.

## Independent native RTDOSE comparison

`tools/validate_dose.py` creates `docs/dose-validation.json` from read-only source
and reference data. The bound is declared as half the native DoseGridScaling plus
8 float32 machine epsilons times maximum dose. This follows native integer
quantization and float32 arithmetic; it is not tuned to observed differences.
Center-lattice interior, positive reference samples, high-dose samples and
outside-center-lattice samples are reported separately.

| Case | Closest plan index | Source version | Positive-inside RMSE Gy | Max error Gy | Quantization-bound positive match |
|---:|---:|---|---:|---:|---|
| 1 | 3 | 4.5.6 | 0.00026643 | 0.00046672 | yes, 1,180,390 samples |
| 2 | 0 | 4.5.5 | 0.04774444 | 0.60634799 | no |
| 3 | 2 | 4.5.5 | 0.04250318 | 0.60572850 | no |
| 4 | 0 | 4.5.5 | 0.03883008 | 0.65823376 | no |
| 5 | 1 | 4.5.5 | 0.04242208 | 0.76898740 | no |
| 6 | 0 | 4.5.5 | 0.04168218 | 0.64443409 | no |
| 7 | 2 | 4.5.6 | 0.00026681 | 0.00046711 | yes, 949,342 samples |

All seven manufacturer reference RTDOSE objects report software version 4.5.6.
Only the matching source plans in cases 1 and 7 contain algorithm payloads created
by 4.5.6 build 145. Cases 2–6 retain 4.5.5 build 137 payloads. This version mismatch
is verified; recalculation or a changed Monte Carlo implementation is a plausible
explanation of their residual, **not a proved cause**. Their stored 4.5.5 dose is
preserved exactly; no recalculation or replacement with reference dose occurs.

Case 1 plans 2 and 3 have identical treatment geometry/MU but different algorithm
payloads. Dose comparison uniquely distinguishes plan 3 as the positive-dose
reference match. The held-out case 7 similarly distinguishes plan 2 from plan 0.
The formula was unchanged for that held-out check.

### Source support and full-volume comparison

The calculation TissueProvider binary is BLZIP/XOR255 little-endian uint16,
shape `[depth,height,width]`, storing HU+1024. Its BinaryFileMD5 matches the decoded
bytes. Its zero samples define empty tissue. The fine TP center affine is the TP
matrix with average Z spacing and half an XY voxel added to its corner origin.
Nearest-voxel sampling of this **source** buffer, followed by `stored_value > 0`,
reproduces native dose support. The calculation crop is independently enforced;
all seven selected buffers have zero positive tissue samples outside that crop.
The algorithm lattice imposes a second, narrower boundary: its half-cell extent.
Stride truncation can leave positive tissue beyond that extent even inside the
TP crop. In case 1, 4,630 explicit source-cache points are tissue-supported but
outside the coarse half-cell extent; every cache value is zero, while unrestricted
clamping would produce up to 3.025194 Gy. The final exporter therefore intersects
tissue/crop support with native coarse coordinates in `[-0.5, n-0.5)` per axis.
The upper face is excluded as a conservative half-open cell convention; no
empirical epsilon is applied at exact boundary ties.

Using the original external ROI's nearest mask is unsuitable: case 1 would wrongly
erase 5,206 positive manufacturer samples, including high-dose boundary values.
The calculation TP buffer has already been resampled from that ROI/CT and is the
relevant support source.

After dose interpolation and source TP support are evaluated at native reference
locations:

- Held-out case 7: **all 1,850,200 reference voxels** pass the quantization bound;
  RMSE 0.00019112 Gy, maximum absolute error 0.00046711 Gy, no support mismatch.
- Case 1: 3,161,087 of 3,161,088 voxels pass. The sole support-boundary difference
  has 0.587499 Gy source dose against reference zero. Its source TP x-coordinate
  lies only 0.000001641 voxel below a half-integer boundary. This is consistent
  with numerical nearest-voxel rounding ambiguity; no arbitrary epsilon or
  reference-derived correction is applied. Full-volume RMSE is 0.00036837 Gy.
- Cases 2–6 retain the measured version-related discrepancy. Source support
  exactly matches their native nonzero support after both the TissueProvider and
  algorithm-cell domains are enforced. No positive native samples are erased.
  Dose-value provenance differences remain reported; exact manufacturer
  equivalence is not claimed.

The JSON separates raw algorithm-dose reconstruction, source-supported native-grid
comparison, and the actual default output-grid comparison.

## Default export grid and its limits

`resample_to_tissue_grid` samples the reconstructed algorithm dose onto the
checksum-verified **fine source TissueProvider lattice**, then applies support on
that target lattice. It does not zero coarse dose samples before interpolation.
The full TP extent is retained, with explicit zero beyond either the calculation
crop or the algorithm half-cell extent.
This source-derived output does not depend on the manufacturer's DICOM grid or
pixels. Metadata records the interpolation and support operation.

A binary support edge is discontinuous. Interpolating the already masked fine
output again onto the manufacturer's different grid cannot exactly reproduce
masking only after interpolation at the manufacturer's voxel centers. This is an
output-grid difference, separately measured below; it must not be hidden by the
excellent direct source-dose-plus-support reconstruction check.

| Case | Fine masked output sampled on native grid: inside RMSE Gy | Maximum absolute error Gy |
|---:|---:|---:|
| 1 | 0.057097 | 5.037336 |
| 2 | 0.049024 | 2.836265 |
| 3 | 0.034745 | 0.855698 |
| 4 | 0.036814 | 1.691451 |
| 5 | 0.033579 | 1.406109 |
| 6 | 0.038205 | 1.093632 |
| 7 | 0.034132 | 1.153022 |

These are observations, not acceptance tolerances. In cases 1 and 7, the original
stored dose and source support have independent near/exact manufacturer evidence,
while the different masked export lattice still has boundary interpolation error.
No clinical interchangeability is inferred from these checks.

## Implementation status

`iplan_dicom/dose.py` provides checked payload decoding, complete selected beam
summation, verified TissueProvider geometry and source-supported fine-grid
resampling. `read_stored_dose` retains native algorithm order without an affine;
`attach_geometry` or `read_dose` returns conventional `[z,y,x]` samples with the
patient affine. `load_tissue_support`, `sample_tissue_support` and
`resample_to_tissue_grid` expose the source support operations.

Fourteen synthetic tests pass, including upper-edge crop remainders, byte layout,
checksum/length failure, invalid values, MU/fraction scaling, duplicate-basename
locality, missing beam dose, unverified algorithms, interpolation-before-clipping,
preventing extrapolation outside the calculation crop and algorithm cell extent,
and rejecting missing/unknown dose-medium provenance. Supported metadata is
currently Monte Carlo selection, 2-mm requested resolution, and observed iPlan
4.5.5/build 137 or 4.5.6/build 145 payloads. The selected TissueProvider must
explicitly declare `XvmcOptions doseResultType="DoseToMedium"`; missing values and
other media are rejected. `StoredDose.dose_deposition_medium` is unset until that
source check succeeds, then remains `DOSE_TO_MEDIUM` through geometry and
resampling.

DICOM Dose Deposition Calculation Medium (3004,0084), VR CS, uses the enumerated
value `DOSE_TO_MEDIUM`, not `MEDIUM`. This is specified by the
[NEMA RT Dose Module](https://dicom.nema.org/medical/dicom/2025e/output/chtml/part03/sect_C.8.8.3.html).
The installed pydicom 2.4.5 dictionary lacks this newer keyword; writers must add
the verified numeric tag explicitly rather than assuming keyword support.

Other modes and incomplete selected
beam graphs fail explicitly. Technical reconstruction is not clinical release.
