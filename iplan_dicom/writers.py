"""Image and structure writers; inputs are already in one verified frame."""
from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from pydicom.dataset import Dataset
from pydicom.sequence import Sequence
from pydicom.uid import generate_uid, ImplicitVRLittleEndian

from .dicom import DicomContext, dsnum, grid_geometry
from .images import ImageVolume, validate_source_headers
from .xbrain import FormatError


def image_datasets(ctx: DicomContext, image: ImageVolume, *, series_number: int = 10,
                   frame_uid: str | None = None) -> list[Dataset]:
    if len(image.headers) != len(image.pixels):
        raise FormatError('Image header and pixel counts differ')
    validate_source_headers(image.headers, image.modality)
    if str(image.headers[0].PatientID) != str(ctx.patient.get('PatientID', '')):
        raise FormatError('Image source patient differs from the selected export patient')
    if (not np.isfinite(image.rescale_slope) or not np.isfinite(image.rescale_intercept)
            or image.rescale_slope == 0):
        raise FormatError('Image intensity calibration is invalid')
    series_uid = generate_uid()
    result = []
    copy_fields = (
        'AcquisitionDate', 'AcquisitionTime', 'ContentDate', 'ContentTime', 'AcquisitionNumber',
        'PatientPosition', 'KVP', 'WindowCenter', 'WindowWidth',
        'ScanningSequence', 'SequenceVariant', 'ScanOptions', 'MRAcquisitionType',
        'RepetitionTime', 'EchoTime', 'InversionTime', 'EchoTrainLength', 'EchoNumbers',
        'MagneticFieldStrength', 'FlipAngle', 'ImagingFrequency', 'ImagedNucleus',
        'NumberOfAverages', 'PercentSampling', 'PercentPhaseFieldOfView', 'PixelBandwidth',
        'ReceiveCoilName', 'TransmitCoilName', 'AcquisitionMatrix', 'InPlanePhaseEncodingDirection',
        'SAR', 'dBdt', 'TemporalPositionIdentifier', 'NumberOfTemporalPositions',
        'LossyImageCompression', 'LossyImageCompressionRatio', 'LossyImageCompressionMethod',
    )
    for index, (pixels, header) in enumerate(zip(image.pixels, image.headers)):
        affine = image.slice_affine(index)
        origin, orientation, spacing, zstep = grid_geometry(affine)
        ds = ctx.new_dataset(image.modality, series_uid=series_uid)
        if frame_uid is not None:
            ds.FrameOfReferenceUID = frame_uid
        for key in copy_fields:
            if hasattr(header, key):
                setattr(ds, key, deepcopy(getattr(header, key)))
        ds.SeriesNumber = series_number
        ds.InstanceNumber = index + 1
        ds.SeriesDescription = 'Converted iPlan ' + image.modality
        ds.ImageType = ['DERIVED', 'SECONDARY', 'AXIAL']
        ds.DerivationDescription = image.derivation_description
        ds.AcquisitionNumber = getattr(header, 'AcquisitionNumber', '')
        ds.SamplesPerPixel = 1
        ds.PhotometricInterpretation = 'MONOCHROME2'
        ds.Rows, ds.Columns = pixels.shape
        ds.BitsAllocated = int(header.BitsAllocated)
        ds.BitsStored = int(header.BitsStored)
        ds.HighBit = int(header.HighBit)
        ds.PixelRepresentation = int(header.PixelRepresentation)
        if ds.PixelRepresentation not in (0, 1) or ds.BitsAllocated not in (8, 16):
            raise FormatError('Output image sample encoding is unsupported')
        lower = -(2 ** (ds.BitsStored - 1)) if ds.PixelRepresentation else 0
        upper = 2 ** (ds.BitsStored - ds.PixelRepresentation)
        if np.any(pixels < lower) or np.any(pixels >= upper):
            raise FormatError('Output pixel values exceed stored bit range')
        dtype = ('<i' if ds.PixelRepresentation else '<u') + str(ds.BitsAllocated // 8)
        ds.PixelData = np.asarray(pixels, dtype=dtype).tobytes()
        ds.ImagePositionPatient = [dsnum(v) for v in origin]
        ds.ImageOrientationPatient = [dsnum(v) for v in orientation]
        ds.PixelSpacing = [dsnum(v) for v in spacing]
        # A subsystem's third axis is a direction, not the acquired thickness.
        # Preserve the source value; DICOM Type 2 permits unknown thickness.
        thickness = getattr(header, 'SliceThickness', None)
        if thickness is None or str(thickness).strip() == '':
            ds.SliceThickness = ''
        else:
            try:
                thickness = float(thickness)
            except (TypeError, ValueError):
                raise FormatError('Original image slice thickness is invalid') from None
            if not np.isfinite(thickness) or thickness <= 0:
                raise FormatError('Original image slice thickness is invalid')
            ds.SliceThickness = dsnum(thickness)
        if image.modality == 'CT':
            ds.KVP = getattr(header, 'KVP', '')
            ds.RescaleIntercept = dsnum(image.rescale_intercept)
            ds.RescaleSlope = dsnum(image.rescale_slope)
            ds.RescaleType = 'HU'
        elif image.rescale_slope != 1 or image.rescale_intercept != 0:
            ds.RescaleIntercept = dsnum(image.rescale_intercept)
            ds.RescaleSlope = dsnum(image.rescale_slope)
        result.append(ds)
    return result


@dataclass
class ROI:
    name: str
    color: tuple[int, int, int]
    contours: list[tuple[int, np.ndarray]]
    interpreted_type: str = ''
    images: list[Dataset] | None = None
    source_grid: Path | None = field(default=None, repr=False)


def sop_reference(ds: Dataset) -> Dataset:
    ref = Dataset()
    ref.ReferencedSOPClassUID = ds.SOPClassUID
    ref.ReferencedSOPInstanceUID = ds.SOPInstanceUID
    return ref


def structure_dataset(ctx: DicomContext, rois: list[ROI], images: list[Dataset]) -> Dataset:
    if not images or not rois:
        raise FormatError('Structures require images and at least one ROI')
    if len({str(d.SeriesInstanceUID) for d in images}) != 1:
        raise FormatError('Contour image references must use one selected series')
    ds = ctx.new_dataset('RTSTRUCT')
    # DS contour arrays can legitimately exceed the 64 KiB explicit-VR length.
    ds.file_meta.TransferSyntaxUID = ImplicitVRLittleEndian
    ds.is_implicit_VR = True
    ds.SeriesNumber = 100
    ds.StructureSetLabel = 'IPLAN_EXPORT'
    ds.StructureSetDate = ds.InstanceCreationDate
    ds.StructureSetTime = ds.InstanceCreationTime
    ds.StructureSetDescription = 'Reconstructed from stored iPlan voxel masks'
    ds.ApprovalStatus = 'UNAPPROVED'
    series_items = {}
    for group in [images] + [roi.images for roi in rois if roi.images is not None]:
        uid = str(group[0].SeriesInstanceUID)
        if any(str(image.FrameOfReferenceUID) != ctx.frame_uid for image in group):
            raise FormatError('ROI images must share the selected planning frame')
        if uid not in series_items:
            series = Dataset()
            series.SeriesInstanceUID = uid
            series.ContourImageSequence = Sequence([sop_reference(image) for image in group])
            series_items[uid] = series
    study = Dataset()
    study.ReferencedSOPClassUID = '1.2.840.10008.3.1.2.3.2'
    study.ReferencedSOPInstanceUID = ctx.study_uid
    study.RTReferencedSeriesSequence = Sequence(list(series_items.values()))
    frame = Dataset()
    frame.FrameOfReferenceUID = ctx.frame_uid
    frame.RTReferencedStudySequence = Sequence([study])
    ds.ReferencedFrameOfReferenceSequence = Sequence([frame])
    ds.StructureSetROISequence = Sequence()
    ds.ROIContourSequence = Sequence()
    ds.RTROIObservationsSequence = Sequence()
    for number, roi in enumerate(rois, 1):
        definition = Dataset()
        definition.ROINumber = number
        definition.ReferencedFrameOfReferenceUID = ctx.frame_uid
        definition.ROIName = roi.name
        definition.ROIGenerationAlgorithm = 'AUTOMATIC'
        definition.ROIGenerationDescription = 'Polygon extraction from stored voxel mask'
        ds.StructureSetROISequence.append(definition)
        contours = Dataset()
        contours.ReferencedROINumber = number
        contours.ROIDisplayColor = [int(v) for v in roi.color]
        contours.ContourSequence = Sequence()
        roi_images = roi.images if roi.images is not None else images
        for plane, vertices in roi.contours:
            points = np.asarray(vertices, dtype=float)
            if not 0 <= plane < len(roi_images) or points.ndim != 2 or points.shape[1] != 3 or len(points) < 3 or not np.isfinite(points).all():
                raise FormatError('Invalid contour geometry or image reference')
            image = roi_images[plane]
            normal = np.cross(np.asarray(image.ImageOrientationPatient[:3], float), np.asarray(image.ImageOrientationPatient[3:], float))
            distances = (points - np.asarray(image.ImagePositionPatient, float)) @ normal
            if np.max(np.abs(distances)) > 0.002:
                raise FormatError('Contour does not lie on its referenced image plane')
            item = Dataset()
            item.ContourGeometricType = 'CLOSEDPLANAR_XOR'
            item.NumberOfContourPoints = len(points)
            item.ContourData = [dsnum(v) for v in points.ravel()]
            item.ContourImageSequence = Sequence([sop_reference(image)])
            contours.ContourSequence.append(item)
        ds.ROIContourSequence.append(contours)
        observation = Dataset()
        observation.ObservationNumber = number
        observation.ReferencedROINumber = number
        observation.RTROIInterpretedType = roi.interpreted_type
        observation.ROIInterpreter = ''
        ds.RTROIObservationsSequence.append(observation)
    return ds
