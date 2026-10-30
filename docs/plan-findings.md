# Native plan mapping evidence

Local read-only investigation and implementation, 2026-10-30. Source identifiers,
patient names, source paths, native UIDs and treatment names are intentionally
excluded. Numerical comparisons below concern anonymized development case 1.
Case 7 remains held out.

## References, revision and schema selection

- Plan folders contain `rtplan.xbrain`, `linacisoccomponent.xbrain`,
  `confarcelement.xbrain`, `confbeamelement.xbrain` and
  `rtobjectcoordinator.xbrain`. Use the current unversioned file's explicit
  `Revision` pointer; do not choose a maximum filename suffix.
- `File/Content/ContentVersion` siblings are alternative serialization schemas.
  The examined element collections contain schemas 4 and 3 with the same beams.
  Concatenating both would duplicate treatment and dose. Select one supported
  schema. Examined versions: RTPlan 4, RTObjectCoordinator 4,
  LINACIsocComponent 3, conformal elements 3/4, Isocenter 1.
- A `Reference` has a `contentType`, `contentUid`, and a directory relative to
  the containing file's folder. Legitimate references include `..` components.
  Resolve these and then enforce containment within the selected case.
- RTPlan points to the primary SliceSet and TissueProvider. A component points
  to RTPlan and Isocenter. Each beam/arc points to its LINACIsocComponent.
  RTObjectCoordinator/RTObject points to each included VoxelObject; this is the
  plan's explicit structure membership, rather than every mask in a case.
- Each treatment element must have exactly one LINACIsocComponent reference.
  Its normalized directory must be the selected plan directory, and the
  resolved object must match that local component. A matching UID in another
  directory does not override the reference location.
- Every populated current `*element.xbrain` collection is inventoried. Unknown
  collection names or object types reject the plan; they are never silently
  left out. Empty collections remain permitted.

## Verified parameter mapping

| Native field | DICOM mapping and evidence |
|---|---|
| Component `componentSortId`, element `sortId` | Sort components first, then elements within each component; generate consecutive BeamNumber values. |
| `ConfArcElement/ConformalField` | One control point per stored field state for dynamic conformal arcs. |
| `desiredGantryRotationDirection=counterclockwise` | Reverse the stored field sequence, use `1 - CumMeterSetWeight`, and set direction `CC`. |
| `desiredGantryRotationDirection=clockwise` | Preserve stored field order and weights, direction `CW`. |
| `ConfBeamElement/ConformalField` | Single static geometry becomes two DICOM control points with weights 0 and 1; native weight -1 is a sentinel. |
| `GantryAngle`, `TableAngle`, `CollimatorAngle` | Degrees modulo 360; become GantryAngle, PatientSupportAngle and BeamLimitingDeviceAngle. |
| `SingleMU` | Exact BeamMeterset. Component RelativeWeighting must not be applied again. |
| Component `DoseRate` | MU/s; multiply by 60 for DICOM MU/min. |
| `JawsState/X1,X2,Y1,Y2` | Millimetres, ASYMX and ASYMY bank order unchanged. |
| `MultileafState/X1,X2` | Count-prefixed vectors; bank 1 followed by bank 2, with machine-defined MLCX boundaries. |
| `Isocenter/NamedPoint/Point/Point_3D` | Count-prefixed homogeneous 4-vector. `diag(-1,-1,+1,1)` maps the examined raw world coordinates to DICOM LPS. |
| RTObject with `hasFractionation=true` | `numberOfFractions` and DosePerFraction give fractionation. Conflicting target fractionations stop export. |

The examined native manufacturer export has **no MLC device or leaf positions**.
It contains ASYMX/ASYMY and an empty APERTURE block descriptor. The machine
configuration's `exportBlockCode=true` is consistent with that export style.
The reconstruction preserves the available native MLC positions; it cannot
claim their agreement with an absent reference representation.

`SingleDose` is not numerically equal to native DICOM `BeamDose` for all beams.
It is retained as source information but is not written as BeamDose. The
calculation point and normalization semantics remain unresolved.

## Machine mapping

Machine `.xmach` files use BLZIP plus byte inversion and XML prefixes that may
be undeclared (`rtmachine`, and `xvmc` for dose profiles). They are independent
documents rather than File/Content wrappers.

Each plan owns its machine snapshot in its own subtree. Development case 6
contains different snapshots with reused identifiers. Indexing all machines
across a case incorrectly conflates those revisions; resolution is scoped to
the selected plan's subtree.

`LINACIsocComponent.machineConfiguration` identifies MachineProfileFileable.
Its `linacId` identifies LINACParameters; `BlockingDeviceId` entries identify
JawsParameters and MultileafParameters. The latter identifies
MultileafGeometry through `multileafGeometry`.

- LINACParameters `exportId` equals native TreatmentMachineName.
- MachineProfileFileable NominalEnergy equals native NominalBeamEnergy.
- LINACParameters SID equals the native SourceAxisDistance.
- Examined geometry has 26 leaf pairs, `isocentric=true`, `Orientation=0`,
  `swapBanksOnExport=false`, and contiguous LeafBeginY/LeafEndY intervals.
  DefinitionSID equals source-axis distance. These checks are required before
  writing the MLC representation.
- The implemented convention is Generic.IEC1217. Other conventions, bank
  transforms, non-isocentric geometries and unsupported treatment types stop
  the plan export rather than substituting geometry.
- Every declared BlockingDeviceId must resolve uniquely. The supported device
  set is exactly one JawsParameters and one MultileafParameters record.
  Both jaw axes must be available and asymmetric, with zero orientation.
  Extra accessories, missing devices, symmetric jaws and unavailable axes
  currently stop export rather than generating assumed jaw/accessory values.
- The component's `machineMD5` is not an MD5 of the whole compressed file,
  decompressed file or decoded XML in the checked instance. Its exact
  serialization/hash scope remains unresolved; it is not falsely presented
  as a verified file integrity check.

## First-case numerical comparison and selection limitation

The matched current treatment contains nine beams. Control-point counts are
9, 4, 9, 9, 9, 9, 4, 2, 2. All metersets, gantry angles, jaws and isocenters match
the native reference numerically; maximum cumulative weight difference is
3.1e-15. Fraction count agrees. Machine energy, source-axis distance and dose
rate agree. This comparison does not validate absent reference MLC geometry.

Four current plan candidates were found. Two have different treatment
parameters; two share exactly the same nine-beam treatment geometry, including
every raw leaf position. These latter candidates also reference the same
SliceSet and have equal fractionation, but all nine stored dose payload hashes
differ. Their RTPlan records differ only in two dosimetry approval child fields.
Therefore RTPLAN parameter matching alone does **not** identify a unique dose
revision. Native RTDOSE comparison is required to disambiguate that choice.

## Implementation and limits

`iplan_dicom/plan.py` exposes `discover_plans`, `read_plan`, and
`populate_rtplan_dataset`. The neutral Plan retains each beam's native element
and source pointer for dose reconstruction. Source RTPlan and source directory
remain in memory for graph resolution. Newly written DICOM is explicitly
UNAPPROVED; no native clinical approval is copied.

Synthetic tests check nonuniform control-point weights under direction reversal,
MLC geometry ordering, static-field duplication, vector count failure, reference
containment, alternate-schema selection, MU preservation and homogeneous
isocenter transformation. The first-case comparison is technical evidence, not
clinical authorization or a complete DICOM conformance audit.

Structural checks of the first six cases also found `arcType=conformal` with
constant MLC aperture and multiple angular control points. This representation
is preserved as gantry motion with the stored fixed aperture; a conformal arc
with varying MLC positions is rejected. A non-treatment candidate in case 2
lacks a unique fractionation and is rejected. Case 7 was not opened for this
mapping work.

After the specification review, the structural check covers 21 current plan
candidates in the first six cases: 17 pass, four reject. The four are the
fractionation case above and three case-6 machine snapshots that explicitly
declare `availableX=false`. Their absent X jaw must not be fabricated as ASYMX.
The two 10-beam dynamic candidates of case 6 remain supported. Static fields
must contain the verified native `CumMeterSetWeight=-1` sentinel; a missing or
arbitrary value is rejected before generating DICOM weights 0 and 1.

The focused plan test suite contains 29 passing synthetic cases, including
reference cardinality/location, unresolved and extra blocking devices, jaw
availability/symmetry, unknown treatment collections and static sentinels.
