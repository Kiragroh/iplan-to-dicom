import numpy as np
import pydicom
import pytest
from pydicom.dataset import Dataset
from pydicom.sequence import Sequence
from pydicom.uid import generate_uid

from iplan_dicom.dicom import (DicomContext, dose_dataset, format_decimal_values,
                               verify_readback, write_dataset, validate_references)
from iplan_dicom.xbrain import FormatError


def test_roundtrip_checks_decoded_pixels_and_modified_geometry(tmp_path):
    ctx = DicomContext.from_source(Dataset())
    ds = dose_dataset(ctx, np.arange(24).reshape(2, 3, 4), np.eye(4), generate_uid())
    path = tmp_path / 'synthetic.dcm'
    write_dataset(ds, path)
    saved = pydicom.dcmread(path)
    verify_readback(ds, saved)
    saved.ImagePositionPatient[0] = 5
    with pytest.raises(FormatError, match='serialization'):
        verify_readback(ds, saved)


def test_nested_decimal_strings_fit_dicom_limit(tmp_path):
    ctx = DicomContext.from_source(Dataset())
    ds = ctx.new_dataset('RTPLAN')
    ds.TableTopVerticalPosition = ''
    beam = Dataset()
    beam.FinalCumulativeMetersetWeight = 0.3333333333333333333
    ds.BeamSequence = Sequence([beam])
    format_decimal_values(ds)
    assert len(str(beam.FinalCumulativeMetersetWeight)) <= 16
    assert float(beam.FinalCumulativeMetersetWeight) == pytest.approx(1 / 3, abs=1e-14)
    path = tmp_path / 'synthetic-plan.dcm'
    write_dataset(ds, path)
    saved = pydicom.dcmread(path)
    assert saved.TableTopVerticalPosition is None
    verify_readback(ds, saved)


def test_rt_frame_mismatch_and_file_meta_class_are_rejected():
    ctx = DicomContext.from_source(Dataset())
    plan = ctx.new_dataset('RTPLAN')
    dose = dose_dataset(ctx, np.ones((2, 2, 2)), np.eye(4), str(plan.SOPInstanceUID))
    dose.FrameOfReferenceUID = generate_uid()
    assert 'RT objects use different planning frames' in validate_references([plan, dose])
    dose.file_meta.MediaStorageSOPClassUID = plan.SOPClassUID
    assert 'File-meta class identifier mismatch' in validate_references([plan, dose])
