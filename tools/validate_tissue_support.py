"""Read-only numerical check of native tissue-buffer support, without identities."""
from pathlib import Path
import argparse
import json
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pydicom
from scipy.ndimage import map_coordinates

from iplan_dicom.plan import objects, discover_plans, read_plan, text, vector
from iplan_dicom.xbrain import FormatError, read_payload, scalar, resolve_inside
from iplan_dicom.dose import read_dose
from iplan_dicom.service import load_headers
from validate_geometry import compare_plan, safe_error


def inspect(source_root, reference_root, case_index):
    case = sorted(path for path in source_root.iterdir() if path.is_dir())[case_index - 1]
    headers = load_headers(case)
    ids = {str(ds.PatientID) for ds in headers.values() if ds.Modality == 'CT'}
    if len(ids) != 1:
        raise FormatError('CT identity pairing is ambiguous')
    reference_rd = reference_plan = None
    for directory in sorted(path for path in reference_root.iterdir() if path.is_dir()):
        for path in directory.rglob('*.dcm'):
            header = pydicom.dcmread(path, stop_before_pixels=True)
            if str(header.PatientID) not in ids:
                break
            if header.Modality == 'RTDOSE':
                reference_rd = pydicom.dcmread(path)
            if header.Modality == 'RTPLAN':
                reference_plan = header
        if reference_rd is not None and reference_plan is not None:
            break
    if reference_rd is None or reference_plan is None:
        raise FormatError('Native reference dose and plan are unavailable')
    candidates = []
    for index, directory in enumerate(discover_plans(case)):
        plan = read_plan(directory, case)
        if compare_plan(plan, reference_plan)['matches_reference']:
            candidates.append((index, directory))
    if not candidates:
        raise FormatError('No source plan matches native control points')
    plan_index, plan_dir = candidates[-1]
    tissue = objects(plan_dir / 'tissueprovider.xbrain', 'TissueProvider')
    if len(tissue) != 1:
        raise FormatError('Tissue-provider selection is ambiguous')
    tissue = tissue[0]
    dimensions = [int(tissue.get(key)) for key in ('tiProWidth', 'tiProHeight', 'tiProDepth')]
    payload = read_payload(resolve_inside(case, plan_dir, text(tissue, 'BinaryFileName')))
    if len(payload) != int(np.prod(dimensions)) * 2:
        raise FormatError('Tissue-buffer dimensions do not match payload')
    buffer = np.frombuffer(payload, dtype='<u2').reshape(dimensions[::-1])
    matrix = vector(tissue, 'TiProP3DSubsyst/Point_3D_Subsystem')[:16].reshape(4, 4)
    matrix[:3, 2] *= scalar(tissue, 'TiProSliceDistanceZ')
    matrix[:3, 3] += .5 * (matrix[:3, 0] + matrix[:3, 1])
    matrix = np.diag([-1., -1., 1., 1.]) @ matrix
    rd = reference_rd
    target = rd.pixel_array.astype(float) * float(rd.DoseGridScaling)
    z, y, x = np.indices(target.shape, dtype=float)
    iop, origin = np.asarray(rd.ImageOrientationPatient, float), np.asarray(rd.ImagePositionPatient, float)
    normal = np.cross(iop[:3], iop[3:])
    positions = (origin[:, None] + iop[:3, None] * (x.ravel() * float(rd.PixelSpacing[1]))
                 + iop[3:, None] * (y.ravel() * float(rd.PixelSpacing[0]))
                 + normal[:, None] * np.asarray(rd.GridFrameOffsetVector, float)[z.astype(int).ravel()])
    inverse = np.linalg.inv(matrix)
    coordinates = inverse[:3, :3] @ positions + inverse[:3, 3, None]
    tissue_values = map_coordinates(buffer, coordinates[::-1], order=0, mode='constant', cval=0)
    inside = tissue_values > 0
    dose = read_dose(plan_dir, case)
    inverse_dose = np.linalg.inv(dose.voxel_to_patient)
    dose_coordinates = inverse_dose[:3, :3] @ positions + inverse_dose[:3, 3, None]
    samples = map_coordinates(dose.values_gy, dose_coordinates[::-1], order=1, mode='nearest')
    flattened = target.ravel()
    positive = flattened > 0
    false_inside = np.flatnonzero(inside & ~positive & (samples > .01))
    false_outside = np.flatnonzero(~inside & positive)
    indices = np.floor(coordinates + .5).astype(int)
    valid = np.all((indices >= 0) & (indices < np.asarray(dimensions)[:, None]), axis=0)
    manual = np.zeros(len(samples), dtype=np.uint16)
    manual[valid] = buffer[indices[2, valid], indices[1, valid], indices[0, valid]]
    result = {'case_index': case_index, 'matching_plan_indices': [index for index, _ in candidates], 'selected_plan_index': plan_index,
              'source_tissue_shape': list(buffer.shape), 'reference_dose_shape': list(target.shape),
              'source_nonzero_tissue_voxels': int((buffer > 0).sum()), 'sampled_support_count': int(inside.sum()),
              'native_positive_count': int(positive.sum()), 'native_positive_outside_support': int(len(false_outside)),
              'native_zero_inside_support_with_dose_gt_0_01Gy': int(len(false_inside)),
              'manual_rounding_vs_scipy_value_difference_count': int((manual != tissue_values).sum()),
              'metrics': {}}
    for label, array in (('unclipped', samples), ('TP_nonzero_nearest_support', samples * inside)):
        error = array - flattened
        result['metrics'][label] = {'all_pixel_rmse_gy': float(np.sqrt(np.mean(error**2))), 'all_pixel_max_abs_gy': float(np.max(np.abs(error))),
                                     'positive_reference_rmse_gy': float(np.sqrt(np.mean(error[positive]**2))),
                                     'positive_reference_max_abs_gy': float(np.max(np.abs(error[positive])))}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', required=True, type=Path)
    parser.add_argument('--reference-root', required=True, type=Path)
    parser.add_argument('--case-indices', default='1,7')
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    source, reference, output = args.source_root.resolve(), args.reference_root.resolve(), args.output.resolve()
    if source in output.parents or reference in output.parents:
        raise FormatError('Report destination must be outside both source roots')
    report = {'read_only_sources': True, 'no_reference_values_used_for_reconstruction': True,
              'status': 'running', 'cases': []}
    for index in map(int, args.case_indices.split(',')):
        try:
            result = inspect(source, reference, index)
            result['status'] = 'compared'
        except Exception as exc:
            result = {'case_index': index, 'status': 'failed', 'error': safe_error(exc)}
        report['cases'].append(result)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')
        print(json.dumps({'case_index': index, 'status': result['status'], 'metrics': result.get('metrics'), 'error': result.get('error')}), flush=True)
    report['status'] = 'completed'
    output.write_text(json.dumps(report, indent=2, allow_nan=False), encoding='utf-8')


if __name__ == '__main__':
    try:
        main()
    except Exception as exc:
        print(json.dumps({'status': 'failed', 'error': safe_error(exc)}), flush=True)
        raise SystemExit(1)
