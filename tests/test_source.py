from pathlib import Path
import stat
import warnings
import zipfile

import pytest

from iplan_dicom.source import prepare_source, ZipLimits
from iplan_dicom.xbrain import FormatError


def archive(tmp_path, entries):
    path = tmp_path / 'synthetic.zip'
    with zipfile.ZipFile(path, 'w', compression=zipfile.ZIP_DEFLATED) as output:
        for name, value in entries:
            with warnings.catch_warnings():
                warnings.filterwarnings('ignore', message='Duplicate name:')
                output.writestr(name, value)
    return path


def track_temporary(monkeypatch):
    import iplan_dicom.source as source
    real, created = source.tempfile.TemporaryDirectory, []
    def create(*args, **kwargs):
        result = real(*args, **kwargs)
        created.append(Path(result.name))
        return result
    monkeypatch.setattr(source.tempfile, 'TemporaryDirectory', create)
    return created


def test_folder_is_not_copied_or_deleted(tmp_path):
    source = tmp_path / 'case'
    source.mkdir()
    (source / 'patient.xbrain').write_text('synthetic')
    with prepare_source(source) as opened:
        assert opened.root == source
        assert not opened.is_zip
    assert (source / 'patient.xbrain').is_file()


def test_nested_single_case_preserves_relative_graph_and_cleans_up(tmp_path):
    path = archive(tmp_path, [('wrapper/case/patient.xbrain', 'synthetic'),
                              ('wrapper/case/plan/rtplan.xbrain', r'..\images\sliceset.xbrain'),
                              ('wrapper/case/images/sliceset.xbrain', 'source')])
    original = path.read_bytes()
    with prepare_source(path) as opened:
        assert opened.is_zip
        assert opened.root.name == 'case'
        root = opened.root
        assert (root / 'plan/rtplan.xbrain').read_text() == r'..\images\sliceset.xbrain'
        assert (root / 'images/sliceset.xbrain').read_text() == 'source'
        assert {member.relative_to(root).as_posix(): member.read_bytes()
                for member in root.rglob('*') if member.is_file()} == {
                    'patient.xbrain': b'synthetic',
                    'plan/rtplan.xbrain': b'..\\images\\sliceset.xbrain',
                    'images/sliceset.xbrain': b'source'}
        with pytest.raises(FormatError, match='temporary'):
            opened.check_output(root.parent / 'output')
    assert not root.exists()
    assert path.read_bytes() == original


def test_nested_collection_is_unwrapped_without_rearranging_cases(tmp_path):
    path = archive(tmp_path, [('outer/wrapper/one/patient.xbrain', 'one'),
                              ('outer/wrapper/two/patient.xbrain', 'two')])
    with prepare_source(path) as opened:
        assert opened.root.name == 'wrapper'
        assert sorted(p.name for p in opened.root.iterdir()) == ['one', 'two']


@pytest.mark.parametrize('name', ['../escape', 'case/../../escape', '/absolute',
    'C:/escape', 'C:escape', r'\\server\share\escape', r'case\..\escape',
    'case/file:stream', 'case/CON.txt', 'case/COM¹.txt', 'case/CONOUT$',
    'case/trailing.', 'case/trailing ', 'case//empty'])
def test_zip_rejects_unsafe_paths_before_extraction(tmp_path, monkeypatch, name):
    created = track_temporary(monkeypatch)
    path = archive(tmp_path, [('case/patient.xbrain', 'synthetic'), (name, 'unsafe')])
    with pytest.raises(FormatError):
        prepare_source(path)
    assert created and all(not root.exists() for root in created)
    assert not (tmp_path / 'escape').exists()


@pytest.mark.parametrize('names', [('case/patient.xbrain', 'case/patient.xbrain'),
    ('case/patient.xbrain', 'CASE/other'), ('case/patient.xbrain', 'case/PATIENT.XBRAIN'),
    ('case/patient.xbrain', 'case'), ('case/patient.xbrain', 'case/é'),
])
def test_duplicate_case_and_file_directory_collisions(tmp_path, names):
    entries = [(name, 'synthetic') for name in names]
    if names[-1] == 'case/é':
        entries.append(('case/e\u0301', 'collision'))
    path = archive(tmp_path, entries)
    with pytest.raises(FormatError, match='collid|duplicate|conflicting'):
        prepare_source(path)


def test_symlink_and_special_members_are_rejected(tmp_path):
    for kind in (stat.S_IFLNK, stat.S_IFIFO):
        info = zipfile.ZipInfo('case/link')
        info.create_system = 3
        info.external_attr = (kind | 0o777) << 16
        path = archive(tmp_path, [('case/patient.xbrain', 'synthetic'), (info, '../outside')])
        with pytest.raises(FormatError, match='special'):
            prepare_source(path)


@pytest.mark.parametrize('limits', [ZipLimits(max_files=1), ZipLimits(max_expanded_bytes=4),
                                   ZipLimits(max_member_bytes=4)])
def test_expansion_and_file_count_limits_cleanup(tmp_path, monkeypatch, limits):
    created = track_temporary(monkeypatch)
    path = archive(tmp_path, [('case/patient.xbrain', 'synthetic'), ('case/other', 'more')])
    with pytest.raises(FormatError, match='limit'):
        prepare_source(path, limits=limits)
    assert all(not root.exists() for root in created)


def test_unsupported_vendor_layout_is_explicit_and_cleanup_occurs(tmp_path, monkeypatch):
    created = track_temporary(monkeypatch)
    path = archive(tmp_path, [('vendor/data.unknown', 'synthetic')])
    with pytest.raises(FormatError, match='Other vendor archive layouts are not yet verified'):
        prepare_source(path)
    assert all(not root.exists() for root in created)


def test_extraction_cleanup_on_consumer_failure(tmp_path):
    path = archive(tmp_path, [('patient.xbrain', 'synthetic')])
    with pytest.raises(RuntimeError):
        with prepare_source(path) as opened:
            root = opened.root
            raise RuntimeError('synthetic consumer failure')
    assert not root.exists()


def test_explicit_directory_followed_by_child_is_valid(tmp_path):
    path = archive(tmp_path, [('case/', ''), ('case/patient.xbrain', 'synthetic')])
    with prepare_source(path) as opened:
        assert (opened.root / 'patient.xbrain').is_file()


def test_corrupt_payload_is_rejected_and_partial_extraction_removed(tmp_path, monkeypatch):
    created = track_temporary(monkeypatch)
    path = tmp_path / 'corrupt.zip'
    with zipfile.ZipFile(path, 'w', compression=zipfile.ZIP_STORED) as output:
        output.writestr('case/patient.xbrain', b'SYNTHETIC_RAW_PAYLOAD')
    path.write_bytes(path.read_bytes().replace(b'SYNTHETIC_RAW_PAYLOAD', b'CORRUPTED_RAW_PAYLOAD'))
    with pytest.raises(FormatError, match='extraction failed'):
        prepare_source(path)
    assert all(not root.exists() for root in created)
