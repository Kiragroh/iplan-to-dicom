import hashlib
from pathlib import Path
from types import SimpleNamespace
import xml.etree.ElementTree as ET

import numpy as np
import pytest

from iplan_dicom.dose import (StoredDose, TissueSupport, decode_dose_samples,
    read_stored_dose, tissue_dose_affine, resample_on_support, sample_tissue_support)
from iplan_dicom.xbrain import FormatError, INVERT


def test_float_endian_axis_order_and_checksum():
    data = np.arange(24, dtype="<f4").tobytes()
    actual = decode_dose_samples(data, (4, 3, 2), hashlib.md5(data).hexdigest())
    assert actual.shape == (2, 3, 4)
    assert actual[1, 2, 3] == 23
    with pytest.raises(FormatError, match="checksum mismatch"):
        decode_dose_samples(data, (4, 3, 2), "0" * 32)
    with pytest.raises(FormatError, match="byte length"):
        decode_dose_samples(data[:-4], (4, 3, 2), hashlib.md5(data[:-4]).hexdigest())


@pytest.mark.parametrize("sample", [-1.0, float("nan"), float("inf")])
def test_invalid_dose_values_fail(sample):
    data = np.array([sample], dtype="<f4").tobytes()
    with pytest.raises(FormatError, match="negative or non-finite"):
        decode_dose_samples(data, (1, 1, 1), hashlib.md5(data).hexdigest())


def beam(directory, values, mu):
    directory.mkdir()
    data = np.asarray(values, dtype="<f4").tobytes()
    (directory / "dose.xbrain").write_bytes(data.translate(INVERT))
    source = ET.fromstring(
        '<ConfBeamElement><AlgoDoseData dimX="2" dimY="1" dimZ="1" '
        'productName="iPlan RT Dose 4.5.6" buildNumber="145">'
        '<BinaryFileName><_ASCII_>dose.xbrain</_ASCII_></BinaryFileName>'
        f'<BinaryFileMD5><_ASCII_>{hashlib.md5(data).hexdigest()}</_ASCII_></BinaryFileMD5>'
        '</AlgoDoseData></ConfBeamElement>')
    return SimpleNamespace(source=source, source_path=directory / "confbeamelement.xbrain", mu=mu)


def source_plan(directory, beams, fractions=5):
    (directory / 'tissueprovider.xbrain').write_text('''<File><Content count="1">
        <ContentVersion value="3"><TissueProvider>
        <XvmcOptions doseResultType="DoseToMedium"/>
        </TissueProvider></ContentVersion></Content></File>''')
    return SimpleNamespace(beams=beams, fractions=fractions, source_directory=directory,
                           source=ET.Element('RTPlan', isPrimaryDoseAlgorithmSelected='false'))


def test_final_mu_fractionation_and_reference_locality(tmp_path):
    # Identical filenames in two directories must resolve at the referring beam.
    plan = source_plan(tmp_path, [beam(tmp_path / "a", [.1, .2], 10),
                                  beam(tmp_path / "b", [.3, .4], 20)])
    result = read_stored_dose(plan, tmp_path)
    np.testing.assert_allclose(result.values_gy, [[[35, 50]]], rtol=1e-7)
    assert result.beam_count == 2
    assert result.dose_deposition_medium == 'DOSE_TO_MEDIUM'
    with pytest.raises(FormatError, match="geometry is not verified"):
        result.require_geometry()


def test_missing_beam_dose_is_not_zero(tmp_path):
    plan = source_plan(tmp_path, [SimpleNamespace(source=ET.Element("ConfBeamElement"))])
    with pytest.raises(FormatError, match="no stored algorithm dose"):
        read_stored_dose(plan, tmp_path)


def test_crop_remainder_anchors_reversed_x_to_upper_edge():
    tissue = ET.fromstring('''<TissueProvider tiProXMin="3" tiProXMax="12"
        tiProYMin="2" tiProYMax="10" tiProZMin="1" tiProZMax="6"
        tiProWidth="14" tiProHeight="12" tiProDepth="8">
        <TiProP3DSubsyst><Point_3D_Subsystem><_DOUBLEVEC_><_ASCII_>
        18,0.75,0,0,10,0,-0.75,0,20,0,0,-1,30,0,0,0,1,2,-1
        </_ASCII_></_DOUBLEVEC_></Point_3D_Subsystem></TiProP3DSubsyst>
        <TiProSliceDistanceZ><_DOUBLE_><_ASCII_>1</_ASCII_></_DOUBLE_></TiProSliceDistanceZ>
        <XvmcOptions><SpatialResolution><_DOUBLE_><_ASCII_>2</_ASCII_></_DOUBLE_></SpatialResolution></XvmcOptions>
        </TissueProvider>''')
    affine = tissue_dose_affine(tissue, (3, 3, 3))
    # X corner 13 minus half of 3, rather than lower corner plus full-cell extent.
    np.testing.assert_allclose(affine @ [0, 0, 0, 1], [-18.625, -17.375, 29, 1])
    np.testing.assert_allclose(affine[:3, :3], np.diag([2.25, 2.25, -2]))
    with pytest.raises(FormatError, match="dimensions disagree"):
        tissue_dose_affine(tissue, (4, 3, 3))


def test_unverified_algorithm_metadata_cannot_be_exported(tmp_path):
    item = beam(tmp_path / 'beam', [.1, .2], 10)
    plan = source_plan(tmp_path, [item])
    item.source.find('AlgoDoseData').set('buildNumber', '999')
    with pytest.raises(FormatError, match="software version"):
        read_stored_dose(plan, tmp_path)
    item.source.find('AlgoDoseData').set('buildNumber', '145')
    plan.source.set('isPrimaryDoseAlgorithmSelected', 'true')
    with pytest.raises(FormatError, match="Monte Carlo algorithm"):
        read_stored_dose(plan, tmp_path)


@pytest.mark.parametrize('medium', ['DoseToWater', 'unknown', ''])
def test_dose_medium_is_checked_not_assumed(tmp_path, medium):
    plan = source_plan(tmp_path, [beam(tmp_path / 'beam', [.1, .2], 10)])
    path = tmp_path / 'tissueprovider.xbrain'
    path.write_text(path.read_text().replace('DoseToMedium', medium))
    with pytest.raises(FormatError, match='deposition medium'):
        read_stored_dose(plan, tmp_path)


def test_tissue_mask_is_applied_after_dose_interpolation():
    # The first coarse point is outside tissue but contributes to interpolation
    # at the first fine tissue point; coarse pre-clipping would halve its dose.
    coarse_affine = np.diag([2., 1., 1., 1.])
    dose = StoredDose(np.array([[[10., 30.]]]), (2, 1, 1), 1, 1, coarse_affine)
    support = TissueSupport(np.array([[[0, 1024, 1024]]], dtype=np.uint16), np.eye(4))
    fine = resample_on_support(dose, support)
    np.testing.assert_allclose(fine.values_gy, [[[0., 20., 30.]]])
    assert fine.support_applied
    assert fine.dimensions_xyz == (3, 1, 1)
    inside = sample_tissue_support(support, np.array([[.4, .6, 1.7, 5.], [0]*4, [0]*4]))
    assert inside.tolist() == [False, True, True, False]


def test_no_dose_extrapolation_outside_source_calculation_crop():
    dose = StoredDose(np.array([[[10., 20.]]]), (2, 1, 1), 1, 1, np.eye(4))
    support = TissueSupport(np.full((1, 1, 5), 1024, dtype=np.uint16), np.eye(4),
                            (1, 0, 0), (2, 0, 0))
    fine = resample_on_support(dose, support)
    np.testing.assert_allclose(fine.values_gy, [[[0., 20., 0., 0., 0.]]])
    accepted = sample_tissue_support(support, np.array([[0., 1., 2., 3.], [0]*4, [0]*4]))
    assert accepted.tolist() == [False, True, True, False]


def test_edge_clamp_stays_inside_coarse_half_cell_domain():
    coarse_affine = np.diag([2., 1., 1., 1.])
    dose = StoredDose(np.array([[[10., 20.]]]), (2, 1, 1), 1, 1, coarse_affine)
    target_affine = np.eye(4)
    target_affine[0, 3] = -.5
    support = TissueSupport(np.full((1, 1, 5), 1024, dtype=np.uint16), target_affine,
                            (0, 0, 0), (4, 0, 0))
    result = resample_on_support(dose, support)
    # Edge clamping is valid at -0.5 and 2.5 mm; 3.5 mm exceeds the final cell.
    np.testing.assert_allclose(result.values_gy, [[[10., 12.5, 17.5, 20., 0.]]])
