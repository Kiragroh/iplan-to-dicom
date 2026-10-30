import numpy as np
import pydicom
import pytest
from pydicom.dataset import Dataset
from pydicom.uid import generate_uid

from iplan_dicom.dicom import DicomContext, dose_dataset, validate_references, write_dataset
from iplan_dicom.xbrain import FormatError


def context():
    ds = Dataset()
    ds.PatientID = 'SYNTHETIC'
    ds.PatientName = 'Synthetic^Only'
    return DicomContext.from_source(ds)


def test_dose_physical_positions_and_roundtrip(tmp_path):
    ctx = context()
    values = np.arange(24, dtype=float).reshape(2, 3, 4) / 3
    affine = np.diag([2., 3., -4., 1.])
    affine[:3, 3] = [10, -20, 30]
    plan_uid = generate_uid()
    ds = dose_dataset(ctx, values, affine, plan_uid)
    out = tmp_path / 'dose.dcm'
    write_dataset(ds, out)
    rd = pydicom.dcmread(out)
    assert np.allclose(rd.pixel_array * float(rd.DoseGridScaling), values, atol=1e-8)
    assert list(rd.ImagePositionPatient) == [10, -20, 30]
    assert list(rd.PixelSpacing) == [3, 2]
    assert list(rd.GridFrameOffsetVector) == [0, -4]
    assert rd.ReferencedRTPlanSequence[0].ReferencedSOPInstanceUID == plan_uid
    assert rd.file_meta.MediaStorageSOPInstanceUID == rd.SOPInstanceUID


@pytest.mark.parametrize('value', [-1, float('nan'), float('inf')])
def test_dose_rejects_invalid_values(value):
    with pytest.raises(FormatError):
        dose_dataset(context(), np.full((2, 2, 2), value), np.eye(4), generate_uid())


def test_dose_rejects_shear():
    affine = np.eye(4)
    affine[0, 1] = .2
    with pytest.raises(FormatError):
        dose_dataset(context(), np.ones((2, 2, 2)), affine, generate_uid())


def test_dangling_plan_reference_is_detected():
    ctx = context()
    plan = ctx.new_dataset('RTPLAN')
    dose = dose_dataset(ctx, np.ones((2, 2, 2)), np.eye(4), str(plan.SOPInstanceUID))
    assert validate_references([plan, dose]) == []
    assert any('reference' in x.lower() for x in validate_references([dose]))


def test_duplicate_instance_is_detected():
    ctx = context()
    ds = ctx.new_dataset('CT')
    assert any('duplicate' in x.lower() for x in validate_references([ds, ds]))


def test_verified_dose_medium_is_written_explicitly():
    ds = dose_dataset(context(), np.ones((2, 2, 2)), np.eye(4), generate_uid(),
                      deposition_medium='DOSE_TO_MEDIUM')
    assert ds[0x30040084].VR == 'CS'
    assert ds[0x30040084].value == 'DOSE_TO_MEDIUM'
    with pytest.raises(FormatError):
        dose_dataset(context(), np.ones((2, 2, 2)), np.eye(4), generate_uid(), deposition_medium='MEDIUM')
