"""Synthetic checks for metadata that would otherwise be erased on export."""
from copy import deepcopy
from types import SimpleNamespace
import xml.etree.ElementTree as ET

import numpy as np
import pytest
from PIL import Image
from pydicom.dataset import Dataset

from iplan_dicom.dicom import DicomContext
from iplan_dicom.images import ImageVolume, load_image_volume
from iplan_dicom.service import load_rois
from iplan_dicom.writers import image_datasets
from iplan_dicom.xbrain import FormatError


def header():
    ds = Dataset()
    ds.PatientID = 'SYNTHETIC'
    ds.SeriesInstanceUID = '1.2.3'
    ds.FrameOfReferenceUID = '1.2.4'
    ds.Modality = 'CT'
    ds.PhotometricInterpretation = 'MONOCHROME2'
    ds.SamplesPerPixel = 1
    ds.PatientPosition = 'HFS'
    ds.Rows = ds.Columns = 2
    ds.ImagePositionPatient = [0, 0, 0]
    ds.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]
    ds.PixelSpacing = [1, 1]
    ds.BitsAllocated, ds.BitsStored, ds.HighBit = 16, 12, 11
    ds.PixelRepresentation = 0
    ds.RescaleSlope, ds.RescaleIntercept = 1, -1024
    ds.SliceThickness = 2.7
    return ds


def source_series(tmp_path):
    series = ET.Element('ImageSeries', imageCount='2', uId='synthetic')
    headers = {}
    for index in range(2):
        ds = header()
        ds.ImagePositionPatient = [0, 0, index * 3]
        headers[str(index)] = ds
        matrix = np.array([[1, 0, 0, -1.5], [0, -1, 0, .5],
                           [0, 0, 1, index * 3], [0, 0, 0, 1.]])
        node = ET.SubElement(series, 'Image', number=str(index), width='2', height='2',
                             headFeetOrientation='HeadFirst', patientOrientation='Supine',
                             dicomInstanceUID=str(index), type='png', modality='CT', bits='16',
                             file=f'{index}.png')
        values = ET.SubElement(ET.SubElement(ET.SubElement(ET.SubElement(node, 'SubSystem'),
                              'Point_3D_Subsystem'), '_DOUBLEVEC_'), '_ASCII_')
        values.text = ','.join(map(str, [18, *matrix.ravel(), 2, -1]))
        Image.fromarray(np.full((2, 2), 500, dtype='u2').byteswap()).save(tmp_path / f'{index}.png')
    return series, headers


@pytest.mark.parametrize(('field', 'value'), [('PatientID', 'DIFFERENT'),
    ('SeriesInstanceUID', '1.2.30'), ('FrameOfReferenceUID', '1.2.40'), ('Modality', 'MR'),
    ('PatientID', ''), ('SeriesInstanceUID', ''), ('FrameOfReferenceUID', '')])
def test_loader_rejects_mixed_or_missing_source_identity(tmp_path, field, value):
    series, headers = source_series(tmp_path)
    setattr(headers['1'], field, value)
    with pytest.raises(FormatError):
        load_image_volume(series, tmp_path, headers)


@pytest.mark.parametrize(('field', 'value'), [('BitsAllocated', 8), ('HighBit', 15),
    ('BitsStored', 17), ('RescaleSlope', 0), ('RescaleSlope', float('inf')),
    ('RescaleIntercept', float('nan'))])
def test_loader_rejects_invalid_encoding_and_calibration(tmp_path, field, value):
    series, headers = source_series(tmp_path)
    for ds in headers.values():
        setattr(ds, field, value)
    with pytest.raises(FormatError):
        load_image_volume(series, tmp_path, headers)


def test_writer_preserves_acquired_thickness_not_unit_subsystem_axis(tmp_path):
    series, headers = source_series(tmp_path)
    image = load_image_volume(series, tmp_path, headers)
    result = image_datasets(DicomContext.from_source(headers['0']), image)
    assert [float(ds.SliceThickness) for ds in result] == [2.7, 2.7]
    assert np.array_equal(result[0].pixel_array, image.pixels[0])
    del headers['0'].SliceThickness
    result = image_datasets(DicomContext.from_source(headers['0']), image)
    assert result[0].SliceThickness == ''


def test_writer_rejects_truncating_pixel_encoding():
    ds = header()
    ds.BitsAllocated = 8  # Stored=12 must not silently cast 500 to 244.
    image = ImageVolume(np.full((1, 2, 2), 500, dtype='u2'), np.eye(4), np.eye(4),
                        (ds,), 'CT', 1, -1024, 'synthetic')
    with pytest.raises(FormatError):
        image_datasets(DicomContext.from_source(ds), image)


def test_writer_rejects_patient_relabeling():
    original = header()
    other = deepcopy(original)
    other.PatientID = 'DIFFERENT'
    image = ImageVolume(np.zeros((1, 2, 2), dtype='u2'), np.eye(4), np.eye(4),
                        (original,), 'CT', 1, -1024, 'synthetic')
    with pytest.raises(FormatError):
        image_datasets(DicomContext.from_source(other), image)


@pytest.mark.parametrize('value', ['MONOCHROME1', 'RGB', ''])
def test_unknown_photometric_is_not_silently_replaced(tmp_path, value):
    series, headers = source_series(tmp_path)
    headers['0'].PhotometricInterpretation = value
    with pytest.raises(FormatError, match='photometric'):
        load_image_volume(series, tmp_path, headers)


def test_lossy_history_is_preserved_not_reclassified(tmp_path):
    series, headers = source_series(tmp_path)
    headers['0'].LossyImageCompression = '01'
    headers['0'].LossyImageCompressionRatio = 3.5
    headers['0'].LossyImageCompressionMethod = 'ISO_10918_1'
    image = load_image_volume(series, tmp_path, headers)
    output = image_datasets(DicomContext.from_source(headers['0']), image)
    assert output[0].LossyImageCompression == '01'
    assert float(output[0].LossyImageCompressionRatio) == 3.5
    assert output[0].LossyImageCompressionMethod == 'ISO_10918_1'
    assert 'LossyImageCompression' not in output[1]


def test_registered_roi_checks_patient_before_registration(tmp_path, monkeypatch):
    import iplan_dicom.images as images
    import iplan_dicom.service as service
    primary, own = header(), header()
    own.PatientID = 'DIFFERENT'
    primary_image = SimpleNamespace(headers=(primary,), source_uid='ct')
    own_image = SimpleNamespace(headers=(own,), source_uid='mr')
    coordinator = ET.fromstring('<RTObjectCoordinator><RTObject><Reference '
        'contentType="VoxelObject" contentUid="roi" directory="mr"/></RTObject></RTObjectCoordinator>')
    grids = ET.Element('SliceSet')
    monkeypatch.setattr(service, 'objects', lambda path, tag: [coordinator] if tag == 'RTObjectCoordinator' else [grids])
    monkeypatch.setattr(service, 'referenced_object', lambda *args: ET.Element('VoxelObject'))
    monkeypatch.setattr(images, 'load_slice_set', lambda *args: own_image)
    monkeypatch.setattr(images, 'register_slice_set', lambda *args: pytest.fail('Registration must not relabel another patient'))
    plan = SimpleNamespace(source_directory=tmp_path)
    def reject_identity(*args):
        raise FormatError('Source image identity review required')
    graph = SimpleNamespace(reconcile=reject_identity)
    with pytest.raises(FormatError, match='identity review required'):
        load_rois(plan, tmp_path, primary_image, headers={}, image_cache={tmp_path: primary_image}, native_graph=graph)
