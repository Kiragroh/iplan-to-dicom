from pathlib import Path

import pytest

from iplan_dicom.service import checked_output
from iplan_dicom.xbrain import FormatError


def test_output_must_be_new_and_outside_source(tmp_path):
    source = tmp_path / 'source'
    source.mkdir()
    assert checked_output(source, tmp_path / 'result') == tmp_path / 'result'
    with pytest.raises(FormatError):
        checked_output(source, source / 'converted')
    existing = tmp_path / 'existing'
    existing.mkdir()
    with pytest.raises(FormatError):
        checked_output(source, existing)
