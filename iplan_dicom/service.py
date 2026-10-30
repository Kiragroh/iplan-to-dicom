"""Local orchestration. Reference exports are never inputs to conversion."""
from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from datetime import datetime
from typing import Callable
import json
import uuid

import numpy as np
import pydicom
from pydicom.sequence import Sequence
from pydicom.uid import generate_uid

from . import __version__
from .dicom import DicomContext, dose_dataset, validate_references, write_dataset, verify_readback
from .plan import discover_plans, objects, read_plan, referenced_object, resolve_directory, text, populate_rtplan_dataset
from .structures import decode_voxel_mask, mask_to_contours
from .writers import ROI, image_datasets, structure_dataset, sop_reference
from .xbrain import FormatError, read_payload, safe_join, scalar


def checked_output(case_dir: Path, output: Path) -> Path:
    case_dir, output = Path(case_dir).resolve(), Path(output).resolve()
    if not case_dir.is_dir():
        raise FormatError('Selected source case is not a directory')
    if output.exists():
        raise FormatError('Choose a new output directory; existing results are never overwritten')
    if output == case_dir or case_dir in output.parents:
        raise FormatError('Output must be outside the source case')
    return output


def acquisition_date(header) -> str:
    """Folder date from the image, never the conversion date."""
    for field in ('AcquisitionDate', 'AcquisitionDateTime', 'SeriesDate', 'StudyDate'):
        raw = str(header.get(field, '')).strip()
        if not raw:
            continue
        value = raw[:8] if field == 'AcquisitionDateTime' else raw
        try:
            if len(value) != 8:
                raise ValueError
            return datetime.strptime(value, '%Y%m%d').strftime('%Y-%m-%d')
        except ValueError:
            raise FormatError('Image acquisition date is invalid') from None
    return 'DateUnknown'


def inspect_case(case_dir: Path, *, include_names: bool = False) -> dict:
    case_dir = Path(case_dir)
    plans = []
    for index, folder in enumerate(discover_plans(case_dir)):
        entry = {'index': index, 'supported_plan': False}
        try:
            plan = read_plan(folder, case_dir)
            entry.update(supported_plan=True, beams=len(plan.beams), fractions=plan.fractions,
                         prescription_per_fraction_gy=plan.prescription_per_fraction,
                         dose_present=all(b.source.find('AlgoDoseData') is not None for b in plan.beams))
            if include_names:
                entry['name'] = text(plan.source, 'PlanName')
        except FormatError as exc:
            entry['reason'] = str(exc)
        plans.append(entry)
    return {'version': __version__, 'plans': plans}


def load_headers(case_dir: Path) -> dict:
    def read(path):
        try:
            return pydicom.dcmread(path, stop_before_pixels=True)
        except Exception:
            raise FormatError('An original DICOM header could not be read') from None
    with ThreadPoolExecutor(max_workers=6) as pool:
        datasets = list(pool.map(read, Path(case_dir).rglob('*.dcm')))
    headers = {}
    for ds in datasets:
        uid = str(ds.get('SOPInstanceUID', ''))
        if not uid:
            raise FormatError('Original image header has no SOP instance identifier')
        if uid in headers:
            raise FormatError('Duplicate original image identifiers require explicit reconciliation')
        headers[uid] = ds
    if not headers:
        raise FormatError('No original DICOM image headers are available')
    return headers


def selected_slice_set(plan, case_dir: Path):
    refs = [r for r in plan.source.findall('Reference') if r.get('contentType') == 'SliceSet']
    if len(refs) != 1:
        raise FormatError('Plan does not reference exactly one image SliceSet')
    directory = resolve_directory(case_dir, plan.source_directory, refs[0].get('directory', '.'))
    return referenced_object(plan.source, 'SliceSet', plan.source_directory, case_dir), directory


def load_rois(plan, case_dir: Path, image, *, headers: dict | None = None,
              image_cache: dict | None = None, native_graph=None) -> tuple[list[ROI], list[dict]]:
    from .images import load_slice_set, register_slice_set
    from .identity import NativeImageGraph
    if native_graph is None:
        native_graph = NativeImageGraph(plan, case_dir)
    if image_cache is None:
        _, primary_directory = selected_slice_set(plan, case_dir)
        image_cache = {primary_directory: image}
    result, metrics = [], []
    seen = set()
    coordinators = objects(plan.source_directory / 'rtobjectcoordinator.xbrain', 'RTObjectCoordinator')
    if len(coordinators) != 1:
        raise FormatError('Expected one structure coordinator for the selected plan')
    for item in coordinators[0].findall('RTObject'):
        refs = [r for r in item.findall('Reference') if r.get('contentType') == 'VoxelObject']
        if len(refs) != 1 or refs[0].get('contentUid') in seen:
            raise FormatError('Ambiguous or repeated ROI membership')
        seen.add(refs[0].get('contentUid'))
        directory = resolve_directory(case_dir, plan.source_directory, refs[0].get('directory', '.'))
        node = referenced_object(item, 'VoxelObject', plan.source_directory, case_dir)
        if directory not in image_cache:
            if headers is None:
                headers = load_headers(case_dir)
            owning_grids = objects(directory / 'sliceset.xbrain', 'SliceSet')
            if len(owning_grids) != 1:
                raise FormatError('ROI source grid is missing or ambiguous')
            own_image = load_slice_set(owning_grids[0], directory, case_dir, headers)
            own_image = native_graph.reconcile(own_image, image)
            image_cache[directory] = register_slice_set(own_image, own_image.source_uid, image.source_uid, plan.source_directory)
        roi_image = image_cache[directory]
        scales = [scalar(node, key) for key in ('XYScalingFactor', 'XYResolutionMultiplier') if node.find(key) is not None]
        if not scales or not np.allclose(scales, scales[0], atol=1e-10, rtol=0):
            raise FormatError('ROI scaling is missing or inconsistent')
        scale = scales[0]
        payload = read_payload(safe_join(directory, node.get('dataFileName', '')))
        mask = decode_voxel_mask(payload, base_width=roi_image.pixels.shape[2],
                                 base_height=roi_image.pixels.shape[1], xy_scale=scale)
        contours = mask_to_contours(mask, roi_image)
        color = tuple(int(node.get(k, '128')) for k in ('colorRed', 'colorGreen', 'colorBlue'))
        if any(v < 0 or v > 255 for v in color):
            raise FormatError('Invalid ROI display color')
        typ = item.get('defaultObjectType', '')
        interpreted = typ if typ in ('PTV', 'CTV', 'GTV', 'EXTERNAL', 'ORGAN') else ''
        result.append(ROI(text(node, 'Name'), color, contours, interpreted, source_grid=directory))
        metrics.append({'roi_index': len(result), 'mask_voxels': mask.voxel_count,
                        'xy_scale': scale, 'contours': len(contours)})
    if not result:
        raise FormatError('Selected plan has no reconstructible structures')
    return result, metrics


def convert_case(case_dir: Path, output: Path, *, plan_index: int,
                 include_mr: bool = True, progress: Callable[[str], None] | None = None) -> dict:
    from .images import load_slice_set, load_image_volume, regularize_image_volume, tissue_image_affine
    from .dose import read_stored_dose, attach_geometry, resample_to_tissue_grid
    from .identity import NativeImageGraph

    report_progress = progress or (lambda message: None)
    case_dir = Path(case_dir).resolve()
    output = checked_output(case_dir, output)
    folders = discover_plans(case_dir)
    if type(plan_index) is not int or not 0 <= plan_index < len(folders):
        raise FormatError('Select an existing plan index')
    report_progress('Plan und Bildreferenzen lesen')
    plan = read_plan(folders[plan_index], case_dir)
    headers = load_headers(case_dir)
    slices, slice_dir = selected_slice_set(plan, case_dir)
    image = load_slice_set(slices, slice_dir, case_dir, headers)
    if image.modality != 'CT':
        raise FormatError('Selected treatment plan must reference a verified CT stack')
    ctx = DicomContext.from_source(image.headers[0])
    original_ct = image_datasets(ctx, image, series_number=11)
    report_progress('Regelmäßiges CT-Raster für die Volumendarstellung erstellen')
    viewer_image = regularize_image_volume(image, tissue_image_affine(plan.source_directory)) if not image.uniform_grid else image
    ct = image_datasets(ctx, viewer_image, series_number=10) if viewer_image is not image else original_ct
    report_progress('Strukturmasken und Konturen rekonstruieren')
    image_cache = {slice_dir: image}
    native_graph = NativeImageGraph(plan, case_dir)
    rois, roi_metrics = load_rois(plan, case_dir, image, headers=headers, image_cache=image_cache,
                                native_graph=native_graph)
    registered_images = {slice_dir: original_ct}
    for directory, roi_image in image_cache.items():
        if directory == slice_dir:
            continue
        registered_images[directory] = image_datasets(ctx, roi_image, series_number=20 + len(registered_images))
    for roi in rois:
        roi.images = registered_images[roi.source_grid]
    rs = structure_dataset(ctx, rois, ct)
    rp = populate_rtplan_dataset(ctx.new_dataset('RTPLAN'), plan, image.raw_to_patient)
    rp.SeriesNumber = 200
    rp.RTPlanDate = rp.InstanceCreationDate
    rp.RTPlanTime = rp.InstanceCreationTime
    rp.ReferencedStructureSetSequence = Sequence([sop_reference(rs)])
    report_progress('Gespeicherte Felddosen zusammenführen')
    stored = attach_geometry(read_stored_dose(plan, case_dir), plan.source_directory)
    stored = resample_to_tissue_grid(stored, plan.source_directory, case_dir)
    rd = dose_dataset(ctx, stored.values_gy, stored.require_geometry(), str(rp.SOPInstanceUID),
                      deposition_medium=stored.dose_deposition_medium)
    # QuickLookRT's default scan is nonrecursive. Keep its primary CT and RT
    # objects together; source geometry and other modalities have own folders.
    datasets = [(f'CT_{i + 1:05d}.dcm', ds) for i, ds in enumerate(ct)]
    if viewer_image is not image:
        folder = 'CT_SourceGeometry_' + acquisition_date(image.headers[0])
        datasets.extend((f'{folder}/{i + 1:05d}.dcm', ds) for i, ds in enumerate(original_ct))
    for index, (directory, images) in enumerate(registered_images.items()):
        if directory != slice_dir:
            source_image = image_cache[directory]
            folder = f'{source_image.modality}_Registered_{acquisition_date(source_image.headers[0])}_{index:02d}'
            datasets.extend((f'{folder}/{i + 1:05d}.dcm', ds) for i, ds in enumerate(images))
    datasets.extend([('RTSTRUCT.dcm', rs), ('RTPLAN.dcm', rp), ('RTDOSE.dcm', rd)])
    warnings = list(plan.warnings)
    if not image.uniform_grid:
        warnings.append('Primary CT is resampled for regular-volume viewing. Exact native per-slice CT and contour references are preserved in a separate series.')
    selected_indices = [int(node.get('imageIndex')) for node in slices.findall('Slices/Slice')]
    skipped_intervals = sum(abs(b - a) > 1 for a, b in zip(selected_indices, selected_indices[1:]))
    if skipped_intervals:
        warnings.append('Native planning CT omits acquired slices. The viewing CT interpolates across these intervals; the complete original CT series is retained separately.')
    mr_count = 0
    original_counts = {'CT': 0, 'MR': 0}
    report_progress('Vollständige ursprüngliche Bildserien im Scanner-Bezugssystem erhalten')
    frames = {}
    for series, series_directory in native_graph.series.values():
        first = series.find('Image')
        modality = first.get('modality') if first is not None else None
        if modality not in ('CT', 'MR') or (modality == 'MR' and not include_mr):
            continue
        original = load_image_volume(series, series_directory, headers)
        original = native_graph.reconcile(original, image)
        source_frame = str(original.headers[0].get('FrameOfReferenceUID', ''))
        if not source_frame:
            raise FormatError('Original image frame identifier is missing')
        frame_uid = frames.setdefault(source_frame, generate_uid())
        original_counts[modality] += 1
        for index, ds in enumerate(image_datasets(ctx, original, series_number=50 + sum(original_counts.values()), frame_uid=frame_uid)):
            datasets.append((f'{modality}_Original_{acquisition_date(original.headers[0])}_{original_counts[modality]:02d}/{index + 1:05d}.dcm', ds))
    mr_count = original_counts['MR']
    warnings.append('Additional original CT/MR series retain their scanner frames and are limited to the selected native reference graph; fusion registrations are not exported as REG objects.')
    if native_graph.reconciled_series:
        warnings.append('Imported image identifiers were reconciled to the planning CT only after native reference and exact demographic checks.')
    errors = validate_references(ds for _, ds in datasets)
    if errors:
        raise FormatError('; '.join(errors))
    report = {'converter_version': __version__, 'status': 'converted_not_clinically_validated',
              'plan_index': plan_index, 'ct_slices': len(ct), 'source_ct_slices': len(original_ct), 'mr_series': mr_count,
              'registered_image_series': len(registered_images),
              'original_image_series': original_counts, 'ct_selection_gap_intervals': skipped_intervals,
              'identity_reconciled_image_grids': len(native_graph.reconciled_series),
              'rois': len(rois), 'beams': len(plan.beams), 'fractions': plan.fractions,
              'dose_shape': list(stored.values_gy.shape), 'dose_max_gy': float(stored.values_gy.max()),
              'dose_sampling': stored.sampling, 'dose_support_applied': stored.support_applied,
              'dose_deposition_medium': stored.dose_deposition_medium,
              'files': len(datasets), 'reference_errors': errors, 'roi_metrics': roi_metrics,
              'warnings': warnings}
    report_progress('DICOM-Dateien schreiben und zurücklesen')
    staging = output.parent / ('.' + output.name + '.incomplete-' + uuid.uuid4().hex[:10])
    staging.mkdir(parents=True, exist_ok=False)
    try:
        for relative, ds in datasets:
            destination = staging / relative
            write_dataset(ds, destination)
            saved = pydicom.dcmread(destination)
            verify_readback(ds, saved)
        (staging / 'conversion-report.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
        staging.rename(output)
    except Exception:
        (staging / 'INCOMPLETE.txt').write_text('Export did not finish. Do not use this incomplete set.\n', encoding='utf-8')
        raise
    report_progress('Konvertierung abgeschlossen')
    return report
