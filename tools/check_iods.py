"""Independent DICOM IOD check; reports tags/counts without patient values.

Install the optional dicom-validator package in a separate validation environment.
The standard edition is pinned explicitly; first use downloads public DICOM specs.
"""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import warnings

import pydicom
from dicom_validator.spec_reader.edition_reader import EditionReader
from dicom_validator.validator.iod_validator import IODValidator
from dicom_validator.validator.error_handler import NullValidationResultHandler


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('export', type=Path)
    parser.add_argument('--standard-cache', type=Path, required=True)
    parser.add_argument('--edition', default='2025e')
    parser.add_argument('--sample-series', action='store_true', help='Check one instance per image series plus every RT object')
    parser.add_argument('--report', type=Path)
    args = parser.parse_args()
    logging.disable(logging.CRITICAL)
    warnings.filterwarnings('ignore', module='pydicom')
    info = EditionReader(str(args.standard_cache)).dicom_info_for_edition(args.edition)
    if info is None:
        raise SystemExit('DICOM standard could not be loaded')
    result = {'edition': args.edition, 'sampled': args.sample_series, 'checked': 0, 'failed': 0, 'by_modality': {}}
    seen = set()
    for path in sorted(args.export.rglob('*.dcm')):
        ds = pydicom.dcmread(path)
        modality = str(ds.get('Modality', 'UNKNOWN'))
        series = str(ds.get('SeriesInstanceUID', ''))
        if args.sample_series and modality in ('CT', 'MR') and series in seen:
            continue
        seen.add(series)
        checked = IODValidator(ds, info, error_handler=NullValidationResultHandler()).validate()
        result['checked'] += 1
        failed = checked.status.name != 'Passed'
        result['failed'] += int(failed)
        group = result['by_modality'].setdefault(modality, {'checked': 0, 'failed': 0, 'errors': []})
        group['checked'] += 1
        group['failed'] += int(failed)
        for module, errors in checked.module_errors.items():
            for tag, error in errors.items():
                entry = {'module': module, 'tag': str(tag), 'code': error.code.name, 'type': str(error.type)}
                if entry not in group['errors']:
                    group['errors'].append(entry)
    if result['checked'] == 0:
        raise SystemExit('No DICOM instances found')
    serialized = json.dumps(result, indent=2)
    if args.report:
        args.report.write_text(serialized, encoding='utf-8')
    print(serialized)
    raise SystemExit(int(result['failed'] > 0))


if __name__ == '__main__':
    try:
        main()
    except Exception:
        print(json.dumps({'status': 'failed', 'error': 'IOD validation failed; source details omitted'}), flush=True)
        raise SystemExit(1)
