"""DICOM dataset context, physical dose encoding and reference checks."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np
from pydicom.dataset import Dataset, FileDataset, FileMetaDataset
from pydicom.sequence import Sequence
from pydicom.uid import (
    CTImageStorage, MRImageStorage, RTDoseStorage, RTPlanStorage,
    RTStructureSetStorage, ExplicitVRLittleEndian, generate_uid,
)
from pydicom.valuerep import format_number_as_ds
from pydicom.multival import MultiValue

from . import __version__
from .xbrain import FormatError

SOP_CLASSES = {'CT': CTImageStorage, 'MR': MRImageStorage, 'RTDOSE': RTDoseStorage,
               'RTPLAN': RTPlanStorage, 'RTSTRUCT': RTStructureSetStorage}


def dsnum(number: float) -> str:
    return format_number_as_ds(float(number))


def format_decimal_values(ds: Dataset) -> None:
    """Serialize computed decimal strings within DICOM's 16-character limit."""
    for element in ds:
        if element.VR == 'SQ':
            for item in element.value:
                format_decimal_values(item)
        elif element.VR == 'DS' and element.value is not None:
            values = element.value if isinstance(element.value, (list, tuple, MultiValue)) else [element.value]
            formatted = [dsnum(value) if str(value).strip() else '' for value in values]
            element.value = formatted if isinstance(element.value, (list, tuple, MultiValue)) else formatted[0]


@dataclass
class DicomContext:
    patient: dict
    study: dict
    study_uid: str
    frame_uid: str

    @classmethod
    def from_source(cls, source: Dataset) -> 'DicomContext':
        patient = {key: getattr(source, key, '') for key in
                   ('PatientName', 'PatientID', 'PatientBirthDate', 'PatientSex')}
        study = {key: getattr(source, key, '') for key in
                 ('StudyDate', 'StudyTime', 'ReferringPhysicianName', 'StudyID', 'AccessionNumber')}
        return cls(patient, study, generate_uid(), generate_uid())

    def new_dataset(self, modality: str, *, series_uid: str | None = None) -> FileDataset:
        if modality not in SOP_CLASSES:
            raise FormatError('Unsupported output modality')
        uid = generate_uid()
        meta = FileMetaDataset()
        meta.MediaStorageSOPClassUID = SOP_CLASSES[modality]
        meta.MediaStorageSOPInstanceUID = uid
        meta.TransferSyntaxUID = ExplicitVRLittleEndian
        meta.ImplementationClassUID = '2.25.208049006731873114005784153381237577921'
        meta.ImplementationVersionName = 'IPLDCM_010'
        ds = FileDataset(None, {}, file_meta=meta, preamble=b'\0' * 128)
        ds.is_little_endian = True
        ds.is_implicit_VR = False
        ds.SpecificCharacterSet = 'ISO_IR 192'
        for key, value in (self.patient | self.study).items():
            setattr(ds, key, value)
        ds.SOPClassUID = meta.MediaStorageSOPClassUID
        ds.SOPInstanceUID = uid
        ds.StudyInstanceUID = self.study_uid
        ds.SeriesInstanceUID = series_uid or generate_uid()
        ds.FrameOfReferenceUID = self.frame_uid
        ds.PositionReferenceIndicator = ''
        ds.Modality = modality
        if modality.startswith('RT'):
            ds.OperatorsName = ''
        ds.SeriesNumber = 1
        ds.InstanceNumber = 1
        ds.Manufacturer = 'BrainScan-Action'
        ds.ManufacturerModelName = 'iPlan DICOM Converter'
        ds.SoftwareVersions = __version__
        now = datetime.now().astimezone()
        ds.InstanceCreationDate = now.strftime('%Y%m%d')
        ds.InstanceCreationTime = now.strftime('%H%M%S.%f')
        ds.TimezoneOffsetFromUTC = now.strftime('%z')
        return ds


def grid_geometry(affine: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    """Return origin, IOP, DICOM row/column spacing and signed frame spacing."""
    affine = np.asarray(affine, dtype=float)
    if affine.shape != (4, 4) or not np.isfinite(affine).all() or not np.allclose(affine[3], [0, 0, 0, 1]):
        raise FormatError('Invalid physical grid transform')
    axes = affine[:3, :3]
    lengths = np.linalg.norm(axes, axis=0)
    if np.any(lengths <= 0):
        raise FormatError('Zero grid spacing')
    directions = axes / lengths
    if not np.allclose(directions.T @ directions, np.eye(3), atol=1e-6, rtol=0):
        raise FormatError('DICOM grid cannot represent shear')
    normal = np.cross(directions[:, 0], directions[:, 1])
    zstep = float(normal @ axes[:, 2])
    return affine[:3, 3], np.r_[directions[:, 0], directions[:, 1]], lengths[[1, 0]], zstep


def dose_dataset(ctx: DicomContext, values_gy: np.ndarray, voxel_to_patient: np.ndarray,
                 plan_uid: str, *, deposition_medium: str | None = None) -> FileDataset:
    values = np.asarray(values_gy, dtype=float)
    if values.ndim != 3 or min(values.shape) <= 0 or not np.isfinite(values).all() or np.any(values < 0):
        raise FormatError('Dose must be a finite nonnegative 3-D physical array')
    origin, orientation, spacing, zstep = grid_geometry(voxel_to_patient)
    ds = ctx.new_dataset('RTDOSE')
    ds.SeriesNumber = 300
    ds.SeriesDescription = 'Converted iPlan physical dose'
    ds.ImageType = ['DERIVED', 'SECONDARY', 'DOSE']
    ds.DoseUnits = 'GY'
    ds.DoseType = 'PHYSICAL'
    if deposition_medium is not None:
        if deposition_medium != 'DOSE_TO_MEDIUM':
            raise FormatError('Unsupported dose deposition medium')
        # Newer standard tag; pydicom 2.4 does not know its keyword yet.
        ds.add_new((0x3004, 0x0084), 'CS', deposition_medium)
    ds.DoseSummationType = 'PLAN'
    ds.NumberOfFrames, ds.Rows, ds.Columns = values.shape
    ds.SamplesPerPixel = 1
    ds.PhotometricInterpretation = 'MONOCHROME2'
    ds.BitsAllocated = 32
    ds.BitsStored = 32
    ds.HighBit = 31
    ds.PixelRepresentation = 0
    ds.ImagePositionPatient = [dsnum(v) for v in origin]
    ds.ImageOrientationPatient = [dsnum(v) for v in orientation]
    ds.PixelSpacing = [dsnum(v) for v in spacing]
    ds.SliceThickness = dsnum(abs(zstep))
    ds.GridFrameOffsetVector = [dsnum(i * zstep) for i in range(values.shape[0])]
    ds.FrameIncrementPointer = [0x3004000C]
    scaling = float(dsnum(max(float(values.max()) / (2 ** 32 - 2), 1e-12)))
    ds.DoseGridScaling = dsnum(scaling)
    ds.PixelData = np.rint(values / scaling).astype('<u4').tobytes()
    reference = Dataset()
    reference.ReferencedSOPClassUID = RTPlanStorage
    reference.ReferencedSOPInstanceUID = plan_uid
    ds.ReferencedRTPlanSequence = Sequence([reference])
    ds.TissueHeterogeneityCorrection = ['IMAGE']
    return ds


def validate_references(datasets: Iterable[Dataset]) -> list[str]:
    """Check the generated set without logging UIDs or patient metadata."""
    datasets = list(datasets)
    errors = []
    ids = [str(d.get('SOPInstanceUID', '')) for d in datasets]
    if any(not uid for uid in ids):
        errors.append('Missing SOP instance identifier')
    if len(ids) != len(set(ids)):
        errors.append('Duplicate SOP instance identifiers')
    available = {str(d.get('SOPInstanceUID', '')): str(d.get('SOPClassUID', '')) for d in datasets}
    patients = {str(d.get('PatientID', '')) for d in datasets}
    studies = {str(d.get('StudyInstanceUID', '')) for d in datasets}
    if len(patients) > 1 or len(studies) > 1:
        errors.append('Mixed patient or study in one export')
    for ds in datasets:
        if hasattr(ds, 'file_meta') and ds.file_meta.get('MediaStorageSOPInstanceUID') != ds.get('SOPInstanceUID'):
            errors.append('File-meta instance identifier mismatch')
        if hasattr(ds, 'file_meta') and ds.file_meta.get('MediaStorageSOPClassUID') != ds.get('SOPClassUID'):
            errors.append('File-meta class identifier mismatch')
        def check_items(item: Dataset) -> None:
            if 'ReferencedSOPInstanceUID' in item:
                uid = str(item.ReferencedSOPInstanceUID)
                sop_class = str(item.get('ReferencedSOPClassUID', ''))
                is_study = sop_class == '1.2.840.10008.3.1.2.3.2' and uid in studies
                if not is_study and uid not in available:
                    errors.append('Unresolved SOP reference within export')
                elif not is_study and sop_class != available[uid]:
                    errors.append('SOP reference class mismatch')
            for element in item:
                if element.VR == 'SQ':
                    for nested in element.value:
                        check_items(nested)
        check_items(ds)
        if ds.get('Modality') == 'RTSTRUCT':
            definitions = list(ds.get('StructureSetROISequence', []))
            numbers = [int(item.ROINumber) for item in definitions]
            if len(numbers) != len(set(numbers)):
                errors.append('Duplicate ROI numbers')
            for key in ('ROIContourSequence', 'RTROIObservationsSequence'):
                if any(int(item.ReferencedROINumber) not in numbers for item in ds.get(key, [])):
                    errors.append('Unresolved ROI number')
            if any(str(item.get('ReferencedFrameOfReferenceUID', '')) != str(ds.get('FrameOfReferenceUID', '')) for item in definitions):
                errors.append('ROI frame differs from the structure set frame')
        if ds.get('Modality') in ('RTPLAN', 'RTDOSE'):
            sequence = 'ReferencedStructureSetSequence' if ds.Modality == 'RTPLAN' else 'ReferencedRTPlanSequence'
            for ref in ds.get(sequence, []):
                linked = next((d for d in datasets if d.get('SOPInstanceUID') == ref.ReferencedSOPInstanceUID), None)
                if linked is not None and ds.get('FrameOfReferenceUID') != linked.get('FrameOfReferenceUID'):
                    errors.append('RT objects use different planning frames')
    return sorted(set(errors))


def verify_readback(expected: Dataset, saved: Dataset) -> None:
    """Check every serialized data element, decoded pixels and file-meta links."""
    def equivalent(left, right):
        if set(left.keys()) != set(right.keys()):
            return False
        for element in left:
            other = right[element.tag]
            if element.VR != other.VR:
                return False
            if element.is_empty and other.is_empty:
                continue  # empty Type 2 numeric values read back as None
            if element.VR == 'SQ':
                if len(element.value) != len(other.value) or not all(equivalent(a, b) for a, b in zip(element.value, other.value)):
                    return False
            elif element != other:
                return False
        return True
    if not equivalent(expected, saved):
        raise FormatError('DICOM data elements changed during serialization')
    if (saved.file_meta.MediaStorageSOPClassUID != saved.SOPClassUID or
            saved.file_meta.MediaStorageSOPInstanceUID != saved.SOPInstanceUID):
        raise FormatError('DICOM readback file-meta mismatch')
    if 'PixelData' in expected and not np.array_equal(expected.pixel_array, saved.pixel_array):
        raise FormatError('Decoded DICOM pixels changed during serialization')


def write_dataset(ds: FileDataset, path: Path) -> None:
    path = Path(path)
    if path.exists():
        raise FormatError('Refusing to overwrite an existing output file')
    path.parent.mkdir(parents=True, exist_ok=True)
    ds.save_as(path, write_like_original=False)
