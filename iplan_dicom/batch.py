"""Organized local exports, verified resume and failure-isolated batch runs."""
from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
import hashlib
import json
from pathlib import Path
import re
import uuid

from . import __version__
from .images import validate_source_headers
from .plan import discover_plans, read_plan, referenced_object, text
from .service import convert_case, load_headers, selected_slice_set
from .source import prepare_source
from .xbrain import FormatError

PROVENANCE = 'batch-provenance.json'


def _hash_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def _files(root: Path) -> list[Path]:
    root = root.resolve()
    paths = sorted(p for p in root.rglob('*') if p.is_file())
    if any(not p.resolve().is_relative_to(root) for p in paths):
        raise FormatError('File reference escapes the selected directory')
    return paths


def _source_state(case: Path) -> dict[str, tuple[int, int]]:
    result = {}
    for path in _files(case):
        stat = path.stat()
        result[path.relative_to(case).as_posix()] = (stat.st_size, stat.st_mtime_ns)
    return result


def source_fingerprint(case: Path) -> tuple[str, dict]:
    """Hash complete case content; names affect the digest but are never logged."""
    case = Path(case).resolve()
    before = _source_state(case)
    digest = hashlib.sha256()
    for relative in sorted(before):
        digest.update(relative.encode('utf-8'))
        digest.update(b'\0')
        digest.update(_hash_file(case / relative).encode('ascii'))
        digest.update(b'\n')
    if _source_state(case) != before:
        raise FormatError('Source changed during fingerprint verification')
    return digest.hexdigest(), before


def discover_cases(source: Path) -> list[Path]:
    source = Path(source).resolve()
    if not source.is_dir():
        raise FormatError('Source is not a case directory or case collection')
    cases = ([source] if (source / 'patient.xbrain').is_file() else
             sorted(p for p in source.iterdir() if p.is_dir() and (p / 'patient.xbrain').is_file()))
    if not cases:
        raise FormatError('No native cases were found')
    return cases


def _component(value: str, *, limit: int) -> str:
    value = value.strip()
    if not value:
        raise FormatError('Output identity or plan label is missing')
    cleaned = re.sub(r'[<>:"/\\|?*\x00-\x1f]', '_', value).rstrip(' .')[:limit].rstrip(' .')
    if not cleaned:
        cleaned = 'Value'
    if cleaned != value:
        cleaned += '_' + hashlib.sha256(value.encode('utf-8')).hexdigest()[:8]
    return cleaned


def _ct_date(header) -> str:
    for field in ('AcquisitionDate', 'AcquisitionDateTime', 'SeriesDate', 'StudyDate'):
        value = str(header.get(field, '')).strip()
        if value:
            try:
                if field == 'AcquisitionDateTime':
                    value = value[:8]
                if len(value) != 8:
                    raise ValueError
                return datetime.strptime(value, '%Y%m%d').strftime('%Y-%m-%d')
            except ValueError:
                raise FormatError('Selected CT date metadata is invalid') from None
    return 'DateUnknown'


@dataclass(frozen=True)
class ExportTarget:
    path: Path
    plan_key: str


def organized_target(case: Path, output_root: Path, plan_index: int) -> ExportTarget:
    """Use only selected native-plan and linked original CT metadata for naming."""
    case, output_root = Path(case).resolve(), Path(output_root).resolve()
    if output_root == case or case in output_root.parents:
        raise FormatError('Output must be outside the source case')
    folders = discover_plans(case)
    if type(plan_index) is not int or not 0 <= plan_index < len(folders):
        raise FormatError('Select an existing plan index')
    plan = read_plan(folders[plan_index], case)
    grid, directory = selected_slice_set(plan, case)
    slices = grid.find('Slices')
    if slices is None:
        raise FormatError('Selected planning grid has no images')
    series = referenced_object(slices, 'ImageSeries', directory, case)
    nodes = sorted(series.findall('Image'), key=lambda n: int(n.get('number', '-1')))
    headers = load_headers(case)
    try:
        indices = [int(s.attrib['imageIndex']) for s in slices.findall('Slice')]
        if not indices or min(indices) < 0 or max(indices) >= len(nodes):
            raise ValueError
        selected = [headers[nodes[i].attrib['dicomInstanceUID']] for i in indices]
    except (KeyError, ValueError):
        raise FormatError('Selected planning CT references are incomplete') from None
    validate_source_headers(selected, 'CT')
    patient = str(selected[0].PatientID)
    label = text(plan.source, 'PlanName')
    key_material = json.dumps([patient, plan.source_directory.relative_to(case).as_posix(),
                               plan.source.get('uId', '')], ensure_ascii=False)
    key = hashlib.sha256(key_material.encode('utf-8')).hexdigest()
    path = (output_root / ('Patient_' + _component(patient, limit=64)) /
            (_ct_date(selected[0]) + '__Plan_' + _component(label, limit=60) + '__' + key[:12])).resolve()
    if not path.is_relative_to(output_root):
        raise FormatError('Organized output escapes the selected output root')
    return ExportTarget(path, key)


def _output_hashes(output: Path) -> dict[str, str]:
    return {p.relative_to(output).as_posix(): _hash_file(p) for p in _files(output)
            if p.relative_to(output).as_posix() != PROVENANCE}


def _verified_resume(target: ExportTarget, fingerprint: str, include_mr: bool) -> dict:
    try:
        saved = json.loads((target.path / PROVENANCE).read_text(encoding='utf-8'))
        expected = {'version': __version__, 'source_fingerprint': fingerprint,
                    'plan_key': target.plan_key, 'include_mr': include_mr}
        if not isinstance(saved, dict) or any(saved.get(k) != v for k, v in expected.items()):
            raise ValueError
        hashes = _output_hashes(target.path)
        if (not hashes or hashes != saved.get('output_sha256') or
                'conversion-report.json' not in hashes or
                any('INCOMPLETE' in name.upper() for name in hashes)):
            raise ValueError
        report = json.loads((target.path / 'conversion-report.json').read_text(encoding='utf-8'))
        if report.get('status') != 'converted_not_clinically_validated':
            raise ValueError
    except (OSError, ValueError, TypeError, KeyError):
        raise FormatError('Existing export failed source or output integrity verification; it was not overwritten') from None
    return {'status': 'resumed_verified', 'files': report.get('files', len(hashes)),
            'plan_index': report.get('plan_index')}


def _convert_plan_directory(case: Path, output_root: Path, *, plan_index: int, include_mr: bool = True,
                 resume: bool = False, progress=None, fingerprint_cache: tuple | None = None) -> dict:
    case = Path(case).resolve()
    target = organized_target(case, output_root, plan_index)
    fingerprint, state = fingerprint_cache or source_fingerprint(case)
    if _source_state(case) != state:
        raise FormatError('Source changed since batch fingerprint verification')
    if target.path.exists():
        if not resume:
            raise FormatError('Organized export already exists; use verified resume or a new output root')
        report = _verified_resume(target, fingerprint, include_mr)
        if _source_state(case) != state:
            raise FormatError('Source changed during resume verification')
        return report
    report = convert_case(case, target.path, plan_index=plan_index, include_mr=include_mr, progress=progress)
    if report.get('status') != 'converted_not_clinically_validated':
        raise FormatError('Conversion did not confirm a completed export')
    if _source_state(case) != state:
        (target.path / 'BATCH_INCOMPLETE.txt').write_text('Source changed during conversion; provenance is unverified.\n', encoding='utf-8')
        raise FormatError('Source changed during conversion; export provenance is unverified')
    provenance = {'version': __version__, 'source_fingerprint': fingerprint, 'plan_key': target.plan_key,
                  'include_mr': include_mr, 'output_sha256': _output_hashes(target.path)}
    with (target.path / PROVENANCE).open('x', encoding='utf-8') as stream:
        json.dump(provenance, stream, indent=2)
    return report


def _convert_batch_directory(source: Path, output_root: Path, *, include_mr: bool = True,
                  resume: bool = False, progress=None) -> dict:
    cases, output_root = discover_cases(source), Path(output_root).resolve()
    report_progress = progress or (lambda message: None)
    if any(output_root == case or case in output_root.parents for case in cases):
        raise FormatError('Batch output must be outside every source case')
    entries = []
    for case_index, case in enumerate(cases):
        report_progress(f'Fall-Index {case_index}: Aktuelle Pläne und Quellen prüfen')
        try:
            folders = discover_plans(case)
            if not folders:
                raise FormatError('Case has no current treatment plans')
            fingerprint = source_fingerprint(case)
        except Exception as exc:
            entries.append({'case_index': case_index, 'plan_index': None, 'status': 'failed',
                            'reason': str(exc) if isinstance(exc, FormatError) else 'Source inspection failed'})
            continue
        for plan_index in range(len(folders)):
            plan_progress = lambda message, ci=case_index, pi=plan_index: report_progress(
                f'Fall-Index {ci}, Plan-Index {pi}: {message}')
            plan_progress('Export oder Fortsetzen prüfen')
            try:
                report = _convert_plan_directory(case, output_root, plan_index=plan_index, include_mr=include_mr,
                                      resume=resume, progress=plan_progress, fingerprint_cache=fingerprint)
                entries.append({'case_index': case_index, 'plan_index': plan_index,
                                'status': 'resumed_verified' if report.get('status') == 'resumed_verified' else 'converted',
                                'files': report.get('files', 0)})
                plan_progress('Vorhandenen Export verifiziert' if report.get('status') == 'resumed_verified' else 'Export abgeschlossen')
            except Exception as exc:
                entries.append({'case_index': case_index, 'plan_index': plan_index, 'status': 'failed',
                                'reason': str(exc) if isinstance(exc, FormatError) else 'Conversion failed; no complete result confirmed'})
                plan_progress('Export fehlgeschlagen; nächster Plan wird geprüft')
    failed = sum(entry['status'] == 'failed' for entry in entries)
    report = {'converter_version': __version__, 'status': 'partial_failure' if failed else 'completed',
              'cases': len(cases), 'converted': sum(e['status'] == 'converted' for e in entries),
              'resumed_verified': sum(e['status'] == 'resumed_verified' for e in entries),
              'failed': failed, 'entries': entries}
    output_root.mkdir(parents=True, exist_ok=True)
    manifest_name = 'batch-manifest-' + datetime.now().strftime('%Y%m%dT%H%M%S') + '-' + uuid.uuid4().hex[:8] + '.json'
    report['manifest'] = manifest_name
    with (output_root / manifest_name).open('x', encoding='utf-8') as stream:
        json.dump(report, stream, indent=2)
    return report


def select_case(root: Path, case_index: int | None = None) -> Path:
    cases = discover_cases(root)
    if case_index is None:
        if len(cases) != 1:
            raise FormatError('Source contains multiple cases; select a case index or use batch')
        return cases[0]
    if type(case_index) is not int or not 0 <= case_index < len(cases):
        raise FormatError('Select an existing case index')
    return cases[case_index]


def convert_plan(case: Path, output_root: Path, *, plan_index: int, include_mr: bool = True,
                 resume: bool = False, progress=None, fingerprint_cache: tuple | None = None,
                 case_index: int | None = None) -> dict:
    """Convert one folder/ZIP case, keeping any private extraction alive until done."""
    with prepare_source(case) as source:
        source.check_output(output_root)
        selected = select_case(source.root, case_index)
        return _convert_plan_directory(selected, output_root, plan_index=plan_index, include_mr=include_mr,
                                       resume=resume, progress=progress, fingerprint_cache=fingerprint_cache)


def convert_batch(source: Path, output_root: Path, *, include_mr: bool = True,
                  resume: bool = False, progress=None) -> dict:
    """Convert a folder or generic native-layout ZIP; always clean owned extraction."""
    with prepare_source(source) as opened:
        opened.check_output(output_root)
        return _convert_batch_directory(opened.root, output_root, include_mr=include_mr,
                                        resume=resume, progress=progress)
