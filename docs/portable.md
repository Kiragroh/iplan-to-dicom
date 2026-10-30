# Portable Windows GUI

Download **iPlan-to-DICOM.exe** from the
[prototype release](https://github.com/Kiragroh/iplan-to-dicom/releases/tag/v0.1.0).
It is a standalone 64-bit Windows executable with Python and its required
libraries bundled. No Python installation, administrator rights or installer is
needed. Initial startup unpacks its runtime into the Windows temporary directory.

1. Start the executable and choose **Ordner wählen** or **ZIP wählen**.
2. For a single test, select a case and explicitly choose its plan variant.
3. Choose an output root outside the source case and select **Ausgewählten Plan
   exportieren**. For a collection choose **Alle Fälle und Pläne exportieren**.
4. Review the result and warnings. Open the generated `RTPLAN.dcm` in QuickLookRT,
   or import the complete plan folder into your test DICOM viewer.

The GUI supports folders and ordinary ZIPs containing the supported native case
layout. A different proprietary ZIP layout has not been tested; a representative
test archive is needed to establish compatibility. The converter does not modify
the selected source. Exports retain patient information.

Please try the prototype on test copies and report the generic failure message
and source software version. Do not attach patient data to public GitHub issues.
This release has no clinical approval or treatment-delivery certification.

Source code is available for adaptation and separate batch pipelines. The CLI
supports the same ZIP/folder processing and resume checks as the GUI.

## Build and verify

Use a clean Python 3.13 x64 environment on Windows:

```powershell
python -m venv build-env
build-env\Scripts\python -m pip install -r requirements-build.txt
.\tools\build_portable.ps1 -Python .\build-env\Scripts\python.exe
```

Build output is `dist/iPlan-to-DICOM.exe`. The release includes a SHA-256 file and
a frozen dependency inventory. The compiled runtime has an optional patient-free
check that exercises PNG decoding, contours, resampling, DICOM roundtrip and Tk:

```powershell
.\iPlan-to-DICOM.exe --self-test "$env:TEMP\iplan-self-test.json"
```

Use a report filename that does not already exist. This confirms bundled runtime
components; it does not replace a real case comparison or a clean-machine test.
