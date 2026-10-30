import json
from pathlib import Path
from types import SimpleNamespace
import xml.etree.ElementTree as ET
import zipfile

import pytest
from pydicom.dataset import Dataset

import iplan_dicom.batch as batch
from iplan_dicom.xbrain import FormatError


def case(tmp_path, name='case'):
    root = tmp_path / name
    root.mkdir()
    (root / 'patient.xbrain').write_text('synthetic source', encoding='utf-8')
    return root


def fake_convert(case, output, *, plan_index, **kwargs):
    output.mkdir(parents=True)
    (output / 'image.dcm').write_bytes(b'SYNTHETIC-PIXELS')
    report = {'status': 'converted_not_clinically_validated', 'files': 1, 'plan_index': plan_index}
    (output / 'conversion-report.json').write_text(json.dumps(report), encoding='utf-8')
    return report


def fixture_target(monkeypatch):
    monkeypatch.setattr(batch, 'organized_target', lambda source, output, index:
                        batch.ExportTarget(output / 'Patient_SYNTHETIC' / f'Plan_{index}', f'key-{index}'))
    monkeypatch.setattr(batch, 'convert_case', fake_convert)


def test_case_discovery_is_bounded_to_immediate_children(tmp_path):
    first = case(tmp_path, 'first')
    other = tmp_path / 'unrelated'
    other.mkdir()
    case(other, 'nested')
    assert batch.discover_cases(first) == [first]
    assert batch.discover_cases(tmp_path) == [first]


def test_source_fingerprint_tracks_names_and_bytes(tmp_path):
    source = case(tmp_path)
    original, _ = batch.source_fingerprint(source)
    (source / 'patient.xbrain').write_text('different source', encoding='utf-8')
    changed, _ = batch.source_fingerprint(source)
    assert changed != original
    (source / 'patient.xbrain').rename(source / 'renamed.xbrain')
    assert batch.source_fingerprint(source)[0] != changed


def test_organized_target_uses_selected_ct_and_explicit_date_priority(tmp_path, monkeypatch):
    source = case(tmp_path)
    folder = source / 'native-plan'
    folder.mkdir()
    plan = ET.fromstring('<RTPlan uId="synthetic-plan"><PlanName><_ASCII_>Synthetic plan</_ASCII_></PlanName></RTPlan>')
    grid = ET.fromstring('<SliceSet><Slices><Slice imageIndex="1"/></Slices></SliceSet>')
    series = ET.fromstring('<ImageSeries><Image number="0" dicomInstanceUID="unused"/><Image number="1" dicomInstanceUID="chosen"/></ImageSeries>')
    header = Dataset()
    header.PatientID = 'SYNTHETIC'
    header.Modality = 'CT'
    header.SeriesInstanceUID, header.FrameOfReferenceUID = '1.2.3', '1.2.4'
    header.BitsAllocated = header.BitsStored = 16
    header.HighBit, header.PixelRepresentation = 15, 0
    header.PhotometricInterpretation, header.SamplesPerPixel = 'MONOCHROME2', 1
    header.AcquisitionDate, header.SeriesDate, header.StudyDate = '20200102', '20200101', '20191231'
    monkeypatch.setattr(batch, 'discover_plans', lambda _: [folder])
    monkeypatch.setattr(batch, 'read_plan', lambda *args: SimpleNamespace(source=plan, source_directory=folder))
    monkeypatch.setattr(batch, 'selected_slice_set', lambda *args: (grid, folder))
    monkeypatch.setattr(batch, 'referenced_object', lambda *args: series)
    monkeypatch.setattr(batch, 'load_headers', lambda *args: {'chosen': header})
    output = tmp_path / 'out'
    first = batch.organized_target(source, output, 0)
    assert first.path.parent.name == 'Patient_SYNTHETIC'
    assert first.path.name.startswith('2020-01-02__Plan_Synthetic plan__')
    assert first == batch.organized_target(source, output, 0)
    header.AcquisitionDateTime = '20200103131415.125+0100'
    assert batch.organized_target(source, output, 0).path.name.startswith('2020-01-02__')
    del header.AcquisitionDate
    assert batch.organized_target(source, output, 0).path.name.startswith('2020-01-03__')
    del header.AcquisitionDateTime
    assert batch.organized_target(source, output, 0).path.name.startswith('2020-01-01__')
    del header.SeriesDate
    assert batch.organized_target(source, output, 0).path.name.startswith('2019-12-31__')
    del header.StudyDate
    assert batch.organized_target(source, output, 0).path.name.startswith('DateUnknown__')
    with pytest.raises(FormatError, match='outside'):
        batch.organized_target(source, source / 'out', 0)


def test_verified_resume_checks_source_output_and_options(tmp_path, monkeypatch):
    source = case(tmp_path)
    fixture_target(monkeypatch)
    output = tmp_path / 'out'
    report = batch.convert_plan(source, output, plan_index=0)
    assert report['status'] == 'converted_not_clinically_validated'
    assert batch.convert_plan(source, output, plan_index=0, resume=True)['status'] == 'resumed_verified'
    with pytest.raises(FormatError, match='already exists'):
        batch.convert_plan(source, output, plan_index=0)
    with pytest.raises(FormatError, match='integrity'):
        batch.convert_plan(source, output, plan_index=0, resume=True, include_mr=False)
    (source / 'patient.xbrain').write_text('changed source', encoding='utf-8')
    with pytest.raises(FormatError, match='integrity'):
        batch.convert_plan(source, output, plan_index=0, resume=True)


@pytest.mark.parametrize('change', ['modify', 'delete', 'extra', 'missing_provenance'])
def test_resume_rejects_any_output_inventory_or_hash_change(tmp_path, monkeypatch, change):
    source = case(tmp_path)
    fixture_target(monkeypatch)
    output = tmp_path / 'out'
    batch.convert_plan(source, output, plan_index=0)
    leaf = batch.organized_target(source, output, 0).path
    if change == 'modify':
        (leaf / 'image.dcm').write_bytes(b'MODIFIED')
    elif change == 'delete':
        (leaf / 'image.dcm').unlink()
    elif change == 'extra':
        (leaf / 'extra.dcm').write_bytes(b'EXTRA')
    else:
        (leaf / batch.PROVENANCE).unlink()
    before = {p.name: p.read_bytes() for p in leaf.iterdir()}
    with pytest.raises(FormatError, match='integrity'):
        batch.convert_plan(source, output, plan_index=0, resume=True)
    assert {p.name: p.read_bytes() for p in leaf.iterdir()} == before


def test_source_mutation_during_conversion_prevents_resume_marker(tmp_path, monkeypatch):
    source = case(tmp_path)
    fixture_target(monkeypatch)
    def changing_convert(case, output, **kwargs):
        result = fake_convert(case, output, **kwargs)
        (case / 'patient.xbrain').write_text('changed during conversion', encoding='utf-8')
        return result
    monkeypatch.setattr(batch, 'convert_case', changing_convert)
    output = tmp_path / 'out'
    with pytest.raises(FormatError, match='Source changed'):
        batch.convert_plan(source, output, plan_index=0)
    leaf = batch.organized_target(source, output, 0).path
    assert (leaf / 'BATCH_INCOMPLETE.txt').exists()
    assert not (leaf / batch.PROVENANCE).exists()


def test_batch_continues_after_individual_failure_and_saves_manifest(tmp_path, monkeypatch):
    source = case(tmp_path)
    fixture_target(monkeypatch)
    monkeypatch.setattr(batch, 'discover_plans', lambda _: [Path('a'), Path('b'), Path('c')])
    def one_failure(case, output, *, plan_index, **kwargs):
        if plan_index == 1:
            raise FormatError('Unsupported synthetic plan')
        return fake_convert(case, output, plan_index=plan_index, **kwargs)
    monkeypatch.setattr(batch, 'convert_case', one_failure)
    output = tmp_path / 'out'
    progress = []
    report = batch.convert_batch(source, output, progress=progress.append)
    assert (report['converted'], report['failed']) == (2, 1)
    assert report['status'] == 'partial_failure'
    assert [e['plan_index'] for e in report['entries']] == [0, 1, 2]
    assert json.loads((output / report['manifest']).read_text()) == report
    assert 'SYNTHETIC' not in json.dumps(report)
    assert any('Fall-Index 0, Plan-Index 1' in line and 'fehlgeschlagen' in line for line in progress)
    assert any('Fall-Index 0, Plan-Index 2' in line and 'abgeschlossen' in line for line in progress)


def test_cli_returns_nonzero_for_batch_partial_failure(monkeypatch, capsys):
    import iplan_dicom.__main__ as cli
    monkeypatch.setattr(cli, 'convert_batch', lambda *args, **kwargs: {'status': 'partial_failure', 'failed': 1})
    assert cli.main(['batch', 'synthetic-source', 'synthetic-output']) == 2
    assert json.loads(capsys.readouterr().out)['failed'] == 1


def test_cli_unexpected_errors_do_not_print_patient_values_or_tracebacks(monkeypatch, capsys):
    import iplan_dicom.__main__ as cli
    def failed(*args, **kwargs):
        raise ValueError('SYNTHETIC_PRIVATE_SOURCE_VALUE')
    monkeypatch.setattr(cli, 'convert_plan', failed)
    assert cli.main(['convert', 'synthetic-source', 'synthetic-output', '--plan-index', '0']) == 4
    captured = capsys.readouterr()
    assert 'SYNTHETIC_PRIVATE_SOURCE_VALUE' not in captured.err
    assert 'Traceback' not in captured.err
    assert json.loads(captured.err)['status'] == 'stopped'


def test_zip_batch_resumes_after_fresh_extraction_and_cleans_temporary_sources(tmp_path, monkeypatch):
    path = tmp_path / 'source.zip'
    with zipfile.ZipFile(path, 'w') as archive:
        archive.writestr('wrapper/case/patient.xbrain', 'synthetic')
    original = path.read_bytes()
    fixture_target(monkeypatch)
    seen = []
    def plans(source):
        seen.append(source)
        assert (source / 'patient.xbrain').is_file()
        return [source / 'plan']
    monkeypatch.setattr(batch, 'discover_plans', plans)
    report = batch.convert_batch(path, tmp_path / 'out')
    assert report['converted'] == 1
    assert all(not folder.exists() for folder in seen)
    second = batch.convert_batch(path, tmp_path / 'out', resume=True)
    assert second['resumed_verified'] == 1
    assert all(not folder.exists() for folder in seen)
    assert path.read_bytes() == original


def test_cli_inspects_selected_case_inside_zip_and_cleans_it(tmp_path, monkeypatch, capsys):
    import iplan_dicom.__main__ as cli
    path = tmp_path / 'cases.zip'
    with zipfile.ZipFile(path, 'w') as archive:
        archive.writestr('wrapper/first/patient.xbrain', 'one')
        archive.writestr('wrapper/second/patient.xbrain', 'two')
    seen = []
    def inspect(source):
        seen.append(source)
        assert (source / 'patient.xbrain').read_text() == 'two'
        return {'plans': []}
    monkeypatch.setattr(cli, 'inspect_case', inspect)
    assert cli.main(['inspect', str(path), '--case-index', '1']) == 0
    assert all(not source.exists() for source in seen)
    assert json.loads(capsys.readouterr().out) == {'plans': []}


def test_zip_batch_failure_still_cleans_owned_extraction(tmp_path, monkeypatch):
    path = tmp_path / 'source.zip'
    with zipfile.ZipFile(path, 'w') as archive:
        archive.writestr('case/patient.xbrain', 'synthetic')
    fixture_target(monkeypatch)
    seen = []
    def plans(source):
        seen.append(source)
        return [source / 'plan']
    def failed(*args, **kwargs):
        raise FormatError('Synthetic conversion failure')
    monkeypatch.setattr(batch, 'discover_plans', plans)
    monkeypatch.setattr(batch, 'convert_case', failed)
    report = batch.convert_batch(path, tmp_path / 'out')
    assert report['status'] == 'partial_failure' and report['failed'] == 1
    assert seen and all(not source.exists() for source in seen)
