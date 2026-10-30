# XBrain image and voxel-mask format findings

Local read-only inspection, 2026-10-30. Format development started on the first
sorted source case. The seventh case was reserved until the final cohort
validation. No identifiers or source filenames are reproduced here.

## Image samples and physical geometry

- ImageSeries/Image `file` references a BLZIP payload containing a standard PNG
  signature after decompression **without XOR inversion**. The matching `.dcm`
  files examined are original metadata headers without PixelData.
- The 16-bit PNG samples require a byte swap after a normal Pillow decode.
  Three CT sample planes then have no samples above the original 12-bit range;
  uncorrected PNG samples have hundreds of thousands of out-of-range samples.
  No intensity offset is added to stored samples. HU use original DICOM
  RescaleSlope and RescaleIntercept.
- SubSystem/Point_3D_Subsystem/_DOUBLEVEC_/_ASCII_ has count18, a row-major 4x4
  affine, and trailing values `2,-1`. The first two axes describe pixel corners,
  not centers. Add half of each in-plane column to obtain pixel centers.
- For observed HFS data, XBrain world to original DICOM patient coordinates is
  diag(-1,-1,+1,1). The resulting native pixel (0,0) position equals original
  DICOM pixel (last-column,0), so native x order is reversed. All four image
  series' first planes match this relation: CT exactly; three MR series within
  0.000062 mm (rounded original direction cosines). The implementation checks
  all four image corners for every plane within 0.002 mm before accepting.
- Keeping native pixel order and writing its explicit DICOM orientation avoids
  an implicit pixel reversal. The stored affine maps voxel centers to patient
  millimetres; its z step is derived from the actual slice origins, not the
  unit third column of SubSystem.
- The selected plan SliceSet can contain a subset of the original image series
  and **individual registered transforms for every plane**. Its imageIndex
  references select source pixels; its own per-plane SubSystem sets physical
  geometry. Using only the original ImageSeries affine would misplace both CT
  and masks. `load_slice_set` preserves these planes through `slice_affines`;
  CT writers and mask contour mapping must use `slice_affine(index)`.
- Original source CT and native plan-export CT differ in origin, slice count,
  and through-plane sampling. Source CT has 191 slices at 1 mm. Reference CT has
  193 slices at approximately 0.9958 mm. Exporting the original quantitative
  selected registered planes is distinct from reproducing the manufacturer's
  plan resampling. First-case selected CT191 loads successfully and is marked
  `uniform_grid=False`; it must not be represented by a single regular affine.
- MR ImageSeries uses its original scanner coordinates. Optional original MR
  remains in separate source-frame groupings. ROI-owned MR masks require the
  owning SliceSet and the selected plan's explicit native fusion into its CT
  planning frame, as detailed below.

## Voxel mask binary representation

After BLZIP decompression **with XOR inversion**:

| Offset | Type | Meaning |
|---|---|---|
| 0 | uint8 | integer byte width: 4 |
| 1 | uint8 | floating-point byte width: 8 |
| 2 | uint32 little endian | matrix element count: 16 |
| 6 | 16 float64 little endian | row-major affine; identity in inspected masks |
| 134 | uint32 little endian | number of planes |
| 138 onward | repeated plane blocks | uint32 run count followed by run lengths |

Run lengths are uint32 little endian and alternate background, foreground,
background, starting with background. Empty planes are one run covering the
whole plane. Runs cover an entire flattened plane in row-major order.
XYScalingFactor increases resolution: scale2 means twice the width and height,
with half the original in-plane spacing. It does not downsample. The physical
base-grid coordinate for a mask index is `(mask_index + .5)/scale - .5`.

Fifteen masks in the initial object document were parsed. Every plane sums to
the expected raster area, every mask has 191 planes, and parsing consumes every
byte exactly. Scale1 masks have 512x512 planes and scale2 masks 1024x1024 planes.
The implementation keeps RLE planes in memory and expands one plane at a time.

Contour extraction uses padded marching squares at 0.5 and retains all rings.
The RTSTRUCT writer must use CLOSEDPLANAR_XOR for holes and disconnected rings.
  Physical mapping uses each selected image plane's center affine. Nonidentity mask
transforms and mismatched image grids fail explicitly.

The binary structure and dimensions are established. Anatomical mask agreement
with manufacturer RTSTRUCT, exact volume conventions, and quantitative pixel
agreement with the manufacturer's resampled export are still separate checks.
The source's stored volume and a direct voxel count need not be identical;
neither is silently substituted for the other.

## First-case RTSTRUCT comparison

Four source masks were paired by exact structure name with the manufacturer's
RTSTRUCT (names stay in protected memory). Source centroid uses foreground
voxel centers transformed by the selected SliceSet's individual plane affine.
Reference centroid uses polygon area weighting on exported axial planes.

| Anonymous matched mask | XY scale | Foreground voxels | Centroid distance (mm) | Largest absolute bounding-extrema difference (mm) |
|---|---:|---:|---:|---:|
| 1 | 1 | 2115157 | 0.0256 | 0.544 |
| 2 | 2 | 76952 | 0.0488 | 0.281 |
| 3 | 2 | 58892 | 0.0335 | 0.649 |
| 4 | 2 | 1662 | 0.1217 | 0.440 |

Bounding extrema compare source **voxel centers** with reference contour
**boundaries**, so half-pixel extent and native resampling are intentionally
not conflated with a registration error. These results support the decoded
run ordering, XY supersampling, handedness and selected-plane registration.

Reference polygon-area times median exported-plane spacing gives volumes
1223207.7, 11111.7, 8484.3 and 232.4 mm3. Source voxel count times original
in-plane area and nominal 1 mm thickness gives 1220961.3, 11105.0, 8498.8 and
239.8 mm3. These are different volume conventions; the small fourth mask is
particularly sensitive to contour resampling. No numerical acceptance threshold
was fitted from these differences. Dice/surface metrics and exhaustive mask
comparison remain outstanding.

## MR-owned structures and native fusion

The first case's full structure coordinator contains 18 masks, including three
stored on MR grids (one on 26 planes, two on 13 planes). Applying the CT191
geometry to these arrays is invalid. Each VoxelObject reference directory has
an owning `sliceset.xbrain`; its imageIndex values select source image samples.
The selected plan's own `fusion.xbrain` defines the registration graph.

Fusion references are ordered source, target; the matrix maps the first
SliceSet world into the second. The loader uses a direct registration where
available, otherwise a uniquely determined shortest composed path. Independent
native direct registrations need not equal longer composed paths. Conflicting
equally short paths or missing links fail explicitly. Only rigid transforms
are accepted. Patient-space application is `RAS_to_LPS @ Fusion @ LPS_to_RAS`.

Direction was checked independently against the three MR-owned native RS
structures. Unregistered centroid errors were 25.53, 22.13 and 22.15 mm;
forward registration gave 0.270, 0.175 and 0.067 mm; inverse registration gave
51.19, 44.15 and 44.15 mm. Registered MR planes and their RTSTRUCT contours must
share the selected CT planning frame. Original MR output, when included, keeps
its original frame and is a separate series.

## Synthetic checks

Thirteen tests cover 16-bit PNG byte order and range, 8-bit PNG preservation,
half-pixel geometry and reflection, alternating RLE with empty planes, invalid
run sums/truncation/trailing data, supersampled contour positions, and masks
touching an image boundary, individual registered plane positions, fusion
composition, direct-registration precedence, and missing/ambiguous links. All four
original first-case series loaded successfully: CT191 unsigned16 and MR96,
MR13, MR26 signed16. These are software checks, not clinical validation.

## Additional regular CT for viewers

The irregular original registered CT planes remain authoritative for source
contour references. Some viewers require constant orientation and regular
spacing. `regularize_image_volume` creates an additional derived CT series on
the native tissue-provider's orthogonal lattice, expanded to cover all original
CT plane corners. Only that grid's geometry is reused; tissue-density samples
are not substituted for the original CT pixels.

For each target-grid normal line, the algorithm finds its physical intersections
with adjacent registered source planes, samples bilinearly within each plane,
and interpolates linearly between those intersections. It rejects crossing
planes. Values outside acquisition coverage use -1024 HU padding; stored values
are rounded to their existing integer encoding. The original image volume is
unchanged. The result carries a specific derivation description and remains
separate from the original registered series and its structure references.

Three synthetic checks establish exact identity-grid preservation, interpolation
of a known linear field across displaced/nonuniform planes, and rejection of
crossing planes. A live first-case run produced a regular 191x514x514 uint16
volume, range0..4095, with orthogonal axes and source-derived spacing. It took
7.56 seconds for resampling. This is not a claim to reproduce the manufacturer's
different regularized CT raster or interpolation implementation.

## Seven-case validation

`tools/validate_geometry.py` produced `geometry-validation.json`, preserving
explicit source errors and ambiguities. Six cases completed all structure and
original-MR checks: 115 of 115 source ROIs matched the native reference names.
The fourth case stopped at the imported-MR identity reconciliation check; it is
not counted as a complete conversion or structure-validation success. Plan MU
and control-point matching was ambiguous in three cases, recorded with all
matching indices and the chosen source pointer ordering.

The held-out seventh case exposed a limitation in the predeclared area-only
centroid/median-spacing volume proxies. Its selected CT skips ten original
images in one interval, producing an 11 mm plane gap. An external-ROI area-only
centroid differs by about 4.19 mm from the regular native reference. An explicit
supplementary analysis weighting source RLE voxel centers by physical slab
thickness reduces the distance to 0.0466 mm and gives a volume-proxy difference
of +0.0865%. The original metrics remain unchanged in the JSON; this is labeled
as an explanatory post-hoc diagnostic, not a revised acceptance threshold.
No export mapping was fitted to this result. The regular viewer reconstruction
interpolates this physical gap and must report it; the omitted original image
samples remain distinct from those interpolated values.

## Native tissue-support comparison

`tools/validate_tissue_support.py` and `tissue-support-validation.json` inspect
the source tissue-provider buffer without using it as diagnostic CT. Its
nonzero nearest-neighbor support explains the manufacturer's dose clipping.
Held-out case7 has exact support agreement at all 949342 positive native dose
voxels; full-grid maximum dose difference is 0.0004671 Gy, consistent with the
native stored-dose quantization. Case1 has one support-boundary mismatch: its
continuous tissue x index lies 0.00000164 voxel from a half-integer boundary.
That exception remains reported; no fitted numerical epsilon was applied.
