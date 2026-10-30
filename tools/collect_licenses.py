"""Collect installed dependency notices for the portable redistribution."""
from importlib.metadata import distribution
from pathlib import Path
import sys

packages = ('numpy', 'scipy', 'pydicom', 'Pillow', 'scikit-image', 'imageio',
            'networkx', 'lazy-loader', 'packaging', 'tifffile', 'pyinstaller')
sections = ['Third-party notices for the portable iPlan-to-DICOM prototype.\n']
python_license = Path(sys.base_prefix) / 'LICENSE.txt'
if python_license.is_file():
    sections.append('Python runtime\n' + python_license.read_text(encoding='utf-8', errors='replace'))
for package in packages:
    dist = distribution(package)
    notices = []
    for file in dist.files or ():
        if any(part.lower().startswith(('license', 'copying', 'notice')) for part in file.parts):
            path = Path(dist.locate_file(file))
            if path.is_file() and path.suffix.lower() not in ('.py', '.pyc', '.png', '.so', '.dll'):
                notices.append(str(file) + '\n' + path.read_text(encoding='utf-8', errors='replace'))
    if not notices:
        raise SystemExit('No redistribution license found for ' + package)
    sections.append('\n' + '=' * 72 + '\n' + package + ' ' + dist.version + '\n' + '\n\n'.join(notices))
target = Path(__file__).resolve().parents[1] / 'docs' / 'third-party-licenses.txt'
target.write_text('\n\n'.join(sections), encoding='utf-8')
print('Dependency notices collected:', len(packages))
