import numpy as np
from pydicom.dataset import Dataset

from iplan_dicom.dicom import DicomContext, validate_references
from pydicom.uid import generate_uid
from iplan_dicom.writers import image_datasets, structure_dataset, ROI
from iplan_dicom.images import ImageVolume


def test_images_contours_share_physical_frame_and_references():
    header = Dataset()
    header.PatientID = 'SYNTHETIC'
    header.Modality = 'CT'
    header.PhotometricInterpretation = 'MONOCHROME2'
    header.SamplesPerPixel = 1
    header.SeriesInstanceUID = generate_uid()
    header.FrameOfReferenceUID = generate_uid()
    header.BitsAllocated = 16
    header.BitsStored = 12
    header.HighBit = 11
    header.PixelRepresentation = 0
    header.PatientPosition = 'HFS'
    affine = np.diag([-2., 3., 4., 1.])
    affine[:3, 3] = [10, 20, -30]
    image = ImageVolume(np.arange(24, dtype=np.uint16).reshape(2, 3, 4), affine, np.eye(4),
                        (header, header), 'CT', 1., -1024., 'synthetic')
    ctx = DicomContext.from_source(header)
    ct = image_datasets(ctx, image, series_number=10)
    assert len(ct) == 2
    assert list(ct[1].ImagePositionPatient) == [10, 20, -26]
    assert list(ct[0].ImageOrientationPatient) == [-1, 0, 0, 0, 1, 0]
    assert np.array_equal(ct[0].pixel_array, image.pixels[0])
    polygon = np.array([[10, 20, -30], [8, 20, -30], [8, 23, -30], [10, 23, -30]])
    rs = structure_dataset(ctx, [ROI('Synthetic ROI', (255, 0, 0), [(0, polygon)])], ct)
    assert validate_references([*ct, rs]) == []
    assert rs.ROIContourSequence[0].ContourSequence[0].ContourGeometricType == 'CLOSEDPLANAR_XOR'
    assert rs.ROIContourSequence[0].ContourSequence[0].ContourImageSequence[0].ReferencedSOPInstanceUID == ct[0].SOPInstanceUID


def test_signed_mr_and_individual_plane_positions():
    header = Dataset()
    header.PatientID = 'SYNTHETIC'
    header.Modality = 'MR'
    header.PhotometricInterpretation = 'MONOCHROME2'
    header.SamplesPerPixel = 1
    header.SeriesInstanceUID = generate_uid()
    header.FrameOfReferenceUID = generate_uid()
    header.BitsAllocated = header.BitsStored = 16
    header.HighBit = 15
    header.PixelRepresentation = 1
    affines = [np.eye(4), np.eye(4)]
    affines[1][:3, 3] = [.2, .1, 2.7]
    image = ImageVolume(np.array([[[-5, 100]], [[200, -20]]], dtype='i2'), np.eye(4), np.eye(4),
                        (header, header), 'MR', 1., 0., 'synthetic', tuple(affines), False)
    output = image_datasets(DicomContext.from_source(header), image)
    assert np.array_equal(output[0].pixel_array, image.pixels[0])
    assert np.allclose(output[1].ImagePositionPatient, [.2, .1, 2.7])


def test_roi_can_reference_registered_mr_in_shared_frame():
    ctx = DicomContext.from_source(Dataset())
    ct, mr = ctx.new_dataset('CT'), ctx.new_dataset('MR')
    for ds in (ct, mr):
        ds.ImagePositionPatient = [0, 0, 0]
        ds.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]
    points = np.array([[0, 0, 0], [1, 0, 0], [1, 1, 0]])
    roi = ROI('MR-defined ROI', (0, 255, 0), [(0, points)], images=[mr])
    rs = structure_dataset(ctx, [roi], [ct])
    assert rs.ROIContourSequence[0].ContourSequence[0].ContourImageSequence[0].ReferencedSOPInstanceUID == mr.SOPInstanceUID
    assert len(rs.ReferencedFrameOfReferenceSequence[0].RTReferencedStudySequence[0].RTReferencedSeriesSequence) == 2
    assert validate_references([ct, mr, rs]) == []
