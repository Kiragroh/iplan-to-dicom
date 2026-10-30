![iPlan to DICOM — cranial stereotactic radiosurgery](assets/banner.png)

# iPlan → DICOM

Local conversion of legacy Brainlab iPlan/XBrain cases into linked **CT, MR,
RTSTRUCT, RTPLAN and RTDOSE** files. Export one selected plan for review or process
a collection of patient folders with a resumable batch run.

**Development release.** Supported source variants are checked explicitly. This
is an independent archival/research converter, not a clinically validated export
or a treatment delivery tool. It reconstructs stored data; it does not recalculate
dose. See [validation and limitations](docs/validation.md).

## Start

**Windows without installation:** download the
[portable GUI](https://github.com/Kiragroh/iplan-to-dicom/releases/download/v0.1.0/iPlan-to-DICOM.exe)
and try a test case. It supports single-plan and batch processing from folders or
ordinary ZIP archives. See the [portable guide](docs/portable.md) and
[release checks](docs/release-evidence.md). Other vendor ZIP layouts remain
unverified; a representative archive may need additional support.

For source installation:

Python 3.11 or newer is required. On Windows, use a Python installation with Tk.

```powershell
git clone https://github.com/Kiragroh/iplan-to-dicom.git
cd iplan-to-dicom
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install -e .
python -m iplan_dicom gui
```

The native desktop interface offers **Ausgewählten Plan exportieren** and **Alle
Fälle und Pläne exportieren**. Choose the original case folder (containing
`patient.xbrain`) or its parent collection, then an output root outside the source
case. The source files are read only.

## One plan

```powershell
python -m iplan_dicom inspect "D:\Legacy\ExampleCase"
python -m iplan_dicom convert "D:\Legacy\ExampleCase" "D:\DICOM-Export" --plan-index 0
```

Select the zero-based index reported by `inspect`. Each current native plan
variant is listed separately; the converter does not guess which is final.
`--without-mr` omits optional original MR series, but retains any image series
required by structure references.

## A complete collection

```powershell
python -m iplan_dicom batch "D:\Legacy" "D:\DICOM-Export" --resume
python -m iplan_dicom batch "D:\Legacy-Cases.zip" "D:\DICOM-Export" --resume
```

Batch discovers case folders directly below the selected root and processes their
current plans independently. An unsupported or inconsistent plan is recorded as
failed while other plans continue. A malformed case inventory fails that case.
The JSON manifest records counts and reasons; a partial failure returns exit code
2. No result is silently marked complete.

`--resume` skips an existing export only when the converter version, selected
plan, options, complete source fingerprint and every output file hash match.
Changed or incomplete results are never overwritten. Use a new output root for a
new conversion version or changed source data.

## Export layout

```text
DICOM-Export/
  Patient_<ID>/
    <CT-acquisition-date>__Plan_<label>__<stable-key>/
      CT_00001.dcm ...
      RTSTRUCT.dcm
      RTPLAN.dcm
      RTDOSE.dcm
      CT_SourceGeometry_<CT-date>/       # if a regular viewing CT was needed
      CT_Original_<CT-date>_01/           # complete acquired scanner-frame CT
      MR_Registered_<MR-date>_01/        # required contour source images
      MR_Original_<MR-date>_01/          # optional original scanner frame
      conversion-report.json
      batch-provenance.json
  batch-manifest-<run>.json
```

Every plan folder is self-contained, including the images referenced by its
structures. Copy the **whole plan folder** when sharing a plan. Names use the
planning CT's acquisition date, with series/study dates as documented fallbacks;
an unknown date remains `DateUnknown`. Other imaging series use their own dates.
Exports contain patient information and are not anonymized.

The primary CT and RT objects share a folder for QuickLookRT's default directory
scan. Open `RTPLAN.dcm` or a CT image. Use recursive folder discovery to include
the additional image series. A regular viewing CT is derived from the original
CT when the native planning slices do not form a regular volume; exact source
geometry remains available separately.

For an archive migration, import each complete plan folder through the receiving
PACS/DICOM system's existing importer and reconcile the batch manifest with its
import report. This project exports files; it does not send C-STORE requests or
alter a live PACS database. Acceptance of RT objects depends on the receiving
system.

## Supported input and safeguards

- Verified XBrain/BLZIP decoding, revision references and checksums.
- Head-first supine CT/MR and explicitly resolved native image registrations.
- Static and dynamic photon plans with verified jaws and 26-pair MLC data in the
  investigated source variants; unsupported treatment components fail explicitly.
- Stored Monte Carlo dose payloads from iPlan 4.5.5/build 137 and
  4.5.6/build 145, with the verified 2-mm calculation mode. Dose is summed with
  native MU and fractionation, then sampled on the source tissue grid.
  Source DoseToMedium provenance is required and preserved.
- Complete reference graph and DICOM pixel/data-element readback before the
  staging folder becomes a completed export.
- Imported image identifiers are reconciled only through verified native
  references and matching nonempty normalized name, birth date and sex.
  Conflicts stop the export for review.

The input parser is deliberately conservative. Other versions, patient positions,
machines, algorithms and source variants are not implied to be supported by the
availability of an installer. Original software and license files are not
distributed with this project.

## Development

```powershell
python -m pip install -e ".[test]"
python -m pytest -q
```

Tests use synthetic data. Private comparison datasets are not part of the
repository. [Format notes](docs/format.md) explain the reconstruction and
[validation](docs/validation.md) separates automated tests, native reference
comparisons and viewer acceptance.

The [banner](assets/banner.md) uses GPT Image-generated synthetic head stereotaxy
illustrations. It is not a clinical screenshot. This project is not affiliated
with Brainlab; product names identify the source format.
