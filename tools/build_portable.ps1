param([string]$Python = 'python')
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
Push-Location -LiteralPath $projectRoot
try {
    & $Python tools/collect_licenses.py
    if ($LASTEXITCODE -ne 0) { throw 'License collection failed' }
    & $Python -m PyInstaller --noconfirm --clean --onefile --windowed --name iPlan-to-DICOM --collect-submodules skimage.measure --collect-data skimage --add-data 'docs/third-party-licenses.txt;.' --add-data 'LICENSE;.' BrainScanConverter.pyw
    if ($LASTEXITCODE -ne 0) { throw 'Portable build failed' }
    Get-FileHash -Algorithm SHA256 -LiteralPath (Join-Path $projectRoot 'dist\iPlan-to-DICOM.exe')
} finally {
    Pop-Location
}
