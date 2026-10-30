# Reconstruction notes

The converter consumes native case directories, not the separately supplied
manufacturer DICOM exports. Those exports are used only by comparison tools.

Typical logical layout (names below are illustrative):

```text
case/
  patient.xbrain
  image-series/
    imageseries.xbrain[,revision]
    original DICOM headers and stored PNG image payloads
  plan-variant/
    rtplan.xbrain[,revision]
    sliceset.xbrain[,revision]
    fusion.xbrain[,revision]
    rtobjectcoordinator.xbrain[,revision]
    tissueprovider.xbrain[,revision]
    treatment components, voxel masks and stored algorithm dose
```

Directories and payload basenames are resolved through the native reference
graph; this sketch is not a fixed directory-name contract. A current Revision
pointer is followed explicitly. ContentVersion branches represent alternative
schema serializations and are not concatenated or treated as plan history.

BLZIP contains a length-prefixed zlib stream. XML and numerical payloads use
inverted bytes; embedded PNG does not. PNG sample byte order and signed MR
encoding are verified against original DICOM headers. Every generated object
receives a new DICOM SOP instance identifier.

The selected plan's per-slice geometry is authoritative. It can differ from the
original scanner geometry and need not form a regular volume. Contours are
reconstructed from voxel masks on their owning CT/MR grids. MR-owned masks use
the selected plan's native rigid registration graph. All contour rings are
preserved with `CLOSEDPLANAR_XOR` semantics.

A separate regular CT series supports volume viewers. It samples the original
image intensities on an orthogonal source-derived lattice; it does not substitute
the calculation's tissue-density array for acquired CT pixels. Selected-plan
gaps are interpolated and reported, while both selected source geometry and the
full acquired CT series remain in separate folders.

RTPLAN includes source beams, MU, fractionation, control points, isocenters, jaws
and MLC positions. Unverified beam types, machine configurations and patient
orientations fail explicitly. Native SingleDose is not assumed to be DICOM
BeamDose. Plans are exported as `UNAPPROVED`.

RTDOSE is reconstructed from checksum-verified stored per-MU Monte Carlo dose,
multiplied by the selected beam MU and number of fractions. DoseToMedium must be
explicit in the source. Fine-grid output interpolates before applying source
tissue support and the verified calculation domain. It is not recalculated.

Detailed evidence and remaining differences:

- [Images, masks and registrations](geometry-findings.md)
- [Plan parameters](plan-findings.md)
- [Dose reconstruction and export-grid effects](dose-findings.md)

## ZIP containers

An ordinary ZIP can contain one native case or one collection of case folders,
including enclosing wrapper directories. Extraction preserves the relative graph
and uses an owned temporary directory. No proprietary installer or program is
executed. Encrypted ZIPs, links, unsafe paths, conflicting names and ambiguous
collections are rejected. Limits are 100,000 entries, 20 GiB expanded total and
2 GiB per member. For larger archives use an already extracted folder collection.

Other vendor archive layouts have not been verified. A ZIP whose contents use a
different serialization or hierarchy may fail even if the archive itself opens.
Normal completion, source changes and normal GUI close clean up extracted files;
forced termination can leave temporary files in the system temporary directory.
