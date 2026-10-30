from copy import deepcopy
from types import SimpleNamespace
import xml.etree.ElementTree as ET

import numpy as np
import pytest
from pydicom.dataset import Dataset

from iplan_dicom.identity import NativeImageGraph
from iplan_dicom.images import ImageVolume
from iplan_dicom.xbrain import FormatError


def write_objects(path, nodes):
    root = ET.Element('File')
    content = ET.SubElement(root, 'Content', count=str(len(nodes)))
    version = ET.SubElement(content, 'ContentVersion', value='1')
    version.extend(deepcopy(nodes))
    path.write_bytes(ET.tostring(root))


@pytest.fixture
def native_images(tmp_path):
    images = []
    for series_index, modality in enumerate(('CT', 'MR')):
        directory = tmp_path / modality
        directory.mkdir()
        series = ET.Element('ImageSeries', uId=modality + '-series', imageCount='2')
        headers = []
        grid = ET.Element('SliceSet', uId=modality + '-grid')
        slices = ET.SubElement(grid, 'Slices', sliceCount='2')
        ET.SubElement(slices, 'Reference', contentType='ImageSeries', contentUid=series.get('uId'), directory='.')
        for index in range(2):
            header = Dataset()
            header.PatientID = modality + '-ID'
            header.PatientName = 'Synthetic^Person'
            header.PatientBirthDate = '20000101'
            header.PatientSex = 'O'
            header.SeriesInstanceUID = f'1.2.{series_index + 1}'
            header.FrameOfReferenceUID = f'1.3.{series_index + 1}'
            header.SOPInstanceUID = f'1.4.{series_index + 1}.{index + 1}'
            header.Modality = modality
            header.BitsAllocated = header.BitsStored = 16
            header.HighBit = 15
            header.PixelRepresentation = 0
            header.PhotometricInterpretation = 'MONOCHROME2'
            header.SamplesPerPixel = 1
            headers.append(header)
            ET.SubElement(series, 'Image', number=str(index), modality=modality,
                          dicomInstanceUID=str(header.SOPInstanceUID))
            ET.SubElement(slices, 'Slice', imageIndex=str(index))
        write_objects(directory / 'imageseries.xbrain', [series])
        write_objects(directory / 'sliceset.xbrain', [grid])
        images.append(ImageVolume(np.zeros((2, 2, 2), dtype='u2'), np.eye(4), np.eye(4),
                                  tuple(headers), modality, 1, 0, modality + '-grid'))
    plan_node = ET.Element('RTPlan')
    ET.SubElement(plan_node, 'Reference', contentType='SliceSet', contentUid='CT-grid', directory='CT')
    fusion = ET.Element('Fusion')
    for modality in ('MR', 'CT'):
        ET.SubElement(fusion, 'Reference', contentType='SliceSet', contentUid=modality + '-grid', directory=modality)
    transform = ET.SubElement(ET.SubElement(ET.SubElement(ET.SubElement(fusion, 'Transformation'),
                              'Matrix_44'), '_DOUBLEVEC_'), '_ASCII_')
    transform.text = ','.join(map(str, [16, *np.eye(4).ravel()]))
    write_objects(tmp_path / 'fusion.xbrain', [fusion])
    graph = NativeImageGraph(SimpleNamespace(source=plan_node, source_directory=tmp_path), tmp_path)
    return graph, images[0], images[1]


def test_verified_native_image_reconciliation_copies_identity_only(native_images):
    graph, ct, mr = native_images
    result = graph.reconcile(mr, ct)
    assert all(h.PatientID == 'CT-ID' for h in result.headers)
    assert all(h.PatientID == 'MR-ID' for h in mr.headers)
    assert result.pixels is mr.pixels
    assert result.headers[0].SOPInstanceUID == mr.headers[0].SOPInstanceUID
    assert result.headers[0].FrameOfReferenceUID == mr.headers[0].FrameOfReferenceUID
    assert graph.reconciled_series == {'MR-grid'}


@pytest.mark.parametrize(('field', 'value'), [('PatientName', 'Other^Person'),
    ('PatientName', ''), ('PatientBirthDate', ''), ('PatientBirthDate', '20000102'),
    ('PatientSex', ''), ('PatientSex', 'F')])
def test_demographic_disagreement_or_missing_data_stops(native_images, field, value):
    graph, ct, mr = native_images
    setattr(mr.headers[1], field, value)
    with pytest.raises(FormatError, match='identity review required'):
        graph.reconcile(mr, ct)
    assert mr.headers[0].PatientID == 'MR-ID'
    assert not graph.reconciled_series


def test_checks_every_canonical_plane_and_normalizes_only_spelling_case_whitespace(native_images):
    graph, ct, mr = native_images
    mr.headers[0].PatientName = 'synthetic ^ person'
    assert graph.reconcile(mr, ct).headers[0].PatientID == 'CT-ID'
    ct.headers[1].PatientBirthDate = '19990101'
    with pytest.raises(FormatError, match='identity review required'):
        graph.reconcile(mr, ct)


def test_arbitrary_series_and_replaced_original_sops_cannot_be_reconciled(native_images):
    graph, ct, mr = native_images
    original_uid = mr.source_uid
    mr.source_uid = 'unlinked'
    with pytest.raises(FormatError, match='not owned'):
        graph.reconcile(mr, ct)
    mr.source_uid = original_uid
    mr.headers[0].SOPInstanceUID = '9.8.7'
    with pytest.raises(FormatError, match='SOP references'):
        graph.reconcile(mr, ct)


def test_original_mr_series_also_requires_exact_graph_ownership(native_images):
    graph, ct, mr = native_images
    mr.source_uid = 'MR-series'
    result = graph.reconcile(mr, ct)
    assert result.headers[0].PatientID == ct.headers[0].PatientID
    assert len(graph.series) == 2
