"""Read-only native-dose comparison; output contains anonymous metrics only.

Reference DICOM pixels are used exclusively here, never by dose reconstruction.
Run from the repository root with Python. Case/plan numbers are transient sorted
directory positions; no identities, paths, UIDs or source pixels are serialized.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
import pydicom
from scipy.ndimage import map_coordinates

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from iplan_dicom.dose import read_dose, load_tissue_support, sample_tissue_support, resample_to_tissue_grid
from iplan_dicom.plan import discover_plans
from iplan_dicom.xbrain import FormatError


def identities(case):
    # A source MR can have a different identifier than its CT. Keep both private
    # and use the reference patient's identity only to pair entire cases.
    found = set()
    for path in case.rglob("*.dcm"):
        ds = pydicom.dcmread(path, stop_before_pixels=True, specific_tags=["PatientID"])
        if ds.get("PatientID"):
            found.add(str(ds.PatientID))
    return found


def reference_doses(root):
    found = {}
    for folder in sorted(p for p in root.iterdir() if p.is_dir()):
        for path in folder.rglob("*"):
            if not path.is_file():
                continue
            try:
                ds = pydicom.dcmread(path, stop_before_pixels=True,
                                     specific_tags=["PatientID", "Modality"])
            except Exception:
                continue
            if ds.get("Modality") == "RTDOSE":
                found.setdefault(str(ds.PatientID), []).append(path)
    return found


def statistics(error, selection, bound):
    values = np.abs(error[selection])
    if not values.size:
        return {"count": 0, "passes_quantization_bound": None}
    signed = error[selection]
    return {"count": int(values.size), "rmse_gy": float(np.sqrt(np.mean(signed**2))),
            "mean_error_gy": float(np.mean(signed)), "max_abs_error_gy": float(values.max()),
            "p95_abs_error_gy": float(np.percentile(values, 95)),
            "exceeding_quantization_bound": int(np.count_nonzero(values > bound)),
            "passes_quantization_bound": bool(np.all(values <= bound))}


def compare(dose, reference, tissue_support=None):
    if (str(reference.DoseUnits), str(reference.DoseType), str(reference.DoseSummationType)) != ("GY", "PHYSICAL", "PLAN"):
        raise FormatError("Reference dose type is unsupported")
    values = reference.pixel_array.astype(np.float64) * float(reference.DoseGridScaling)
    # Native reference quantization permits half a stored integer step. The
    # independent float allowance is 8 float32 epsilons times maximum dose,
    # covering native intermediate float32 rounding; it is not tuned per case.
    floating = 8 * np.finfo(np.float32).eps * max(1., float(values.max()), float(dose.values_gy.max()))
    bound = float(reference.DoseGridScaling) / 2 + floating
    orientation = np.asarray(reference.ImageOrientationPatient, float)
    normal = np.cross(orientation[:3], orientation[3:])
    origin = np.asarray(reference.ImagePositionPatient, float)
    spacing = np.asarray(reference.PixelSpacing, float)
    offsets = np.asarray(reference.GridFrameOffsetVector, float)
    if abs(offsets[0]) > 1e-8:
        raise FormatError("Only relative reference dose frame offsets are supported")
    inverse = np.linalg.inv(dose.require_geometry())
    error = np.empty(values.shape, dtype=np.float64)
    inside = np.empty(values.shape, dtype=bool)
    supported_error = np.empty(values.shape, dtype=np.float64) if tissue_support is not None else None
    source_support = np.empty(values.shape, dtype=bool) if tissue_support is not None else None
    zz, yy, xx = values.shape
    row, col = np.indices((yy, xx))
    base = origin[:, None, None] + orientation[:3, None, None]*col*spacing[1] + orientation[3:, None, None]*row*spacing[0]
    for frame in range(zz):
        points = base + normal[:, None, None] * offsets[frame]
        coordinates = inverse[:3, :3] @ points.reshape(3, -1) + inverse[:3, 3, None]
        shape_xyz = np.array(dose.values_gy.shape[::-1])[:, None]
        # Outside the center lattice is reported separately; no acceptance is
        # inferred from clamping extrapolation beyond reconstructed grid support.
        support = np.all((coordinates >= 0) & (coordinates <= shape_xyz-1), axis=0)
        predicted = map_coordinates(dose.values_gy, coordinates[::-1], order=1, mode="nearest")
        error[frame] = predicted.reshape(yy, xx) - values[frame]
        inside[frame] = support.reshape(yy, xx)
        if source_support is not None:
            source_support[frame] = sample_tissue_support(tissue_support, points.reshape(3, -1)).reshape(yy, xx)
            cell_domain = np.all((coordinates >= -.5) & (coordinates < shape_xyz-.5), axis=0)
            source_support[frame] &= cell_domain.reshape(yy, xx)
            supported_error[frame] = np.where(source_support[frame], predicted.reshape(yy, xx), 0) - values[frame]
    result = {"source_shape_zyx": list(dose.values_gy.shape),
            "source_algorithm_versions": list(dose.algorithm_versions),
            "source_dose_deposition_medium": dose.dose_deposition_medium,
            "reference_software_version": str(reference.get("SoftwareVersions", "unspecified")),
            "reference_shape_zyx": list(values.shape),
            "source_max_gy": float(dose.values_gy.max()), "reference_max_gy": float(values.max()),
            "reference_quantum_gy": float(reference.DoseGridScaling),
            "floating_allowance_gy": float(floating), "quantization_bound_gy": float(bound),
            "all_inside": statistics(error, inside, bound),
            "positive_inside": statistics(error, inside & (values > 0), bound),
            "zero_reference_inside": statistics(error, inside & (values == 0), bound),
            "above_1gy_inside": statistics(error, inside & (values > 1), bound),
            "positive_outside_center_lattice": statistics(error, ~inside & (values > 0), bound),
            "zero_reference_outside_center_lattice_count": int(np.count_nonzero(~inside & (values == 0)))}
    if source_support is not None:
        result["source_support_at_reference_voxels"] = {
            "all_voxels": statistics(supported_error, np.ones(values.shape, dtype=bool), bound),
            "source_support_count": int(source_support.sum()),
            "reference_positive_outside_source_support": int(np.count_nonzero(~source_support & (values > 0))),
            "reference_zero_inside_source_support": int(np.count_nonzero(source_support & (values == 0)))}
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    refs = reference_doses(args.reference)
    report = {"validation_kind": "independent_native_dose_quantization", "case7_held_out_until_affine_fixed": True,
              "source_pixels_from_reference": False,
              "interpolation": "trilinear; clamp edge within calculation crop; center-lattice outside reported separately",
              "source_support_rule": "nearest source TP HU+1024 greater than zero, inside TP crop and native coarse half-cell domain [-0.5,n-0.5)",
              "default_output_grid": "fine source TissueProvider lattice; interpolate dose before support masking",
              "support_boundary_precision": "double coordinates; no fitted half-index epsilon", "cases": []}
    cases = sorted(p for p in args.source.iterdir() if p.is_dir())
    for ci, case in enumerate(cases, 1):
        matches = {p for pid in identities(case) for p in refs.get(pid, [])}
        record = {"case_index": ci, "plans": []}
        report["cases"].append(record)
        if len(matches) != 1:
            record["blocked"] = "Reference case identity is missing or ambiguous"
            continue
        reference = pydicom.dcmread(next(iter(matches)))
        for pi, plan in enumerate(discover_plans(case)):
            result = {"plan_index": pi}
            record["plans"].append(result)
            try:
                dose = read_dose(plan, case)
                support = load_tissue_support(plan, case)
                result.update(beam_count=dose.beam_count, fractions=dose.fractions,
                              comparison=compare(dose, reference, support))
                summary = result["comparison"]["above_1gy_inside"]
                print(json.dumps({"case": ci, "plan": pi, **summary}), flush=True)
            except FormatError as error:
                result["blocked"] = str(error)
                print(json.dumps({"case": ci, "plan": pi, "blocked": str(error)}), flush=True)
        candidates = [r for r in record["plans"] if "comparison" in r]
        if candidates:
            best = min(candidates, key=lambda r: r["comparison"]["positive_inside"].get("rmse_gy", float("inf")))
            record["closest_plan_index"] = best["plan_index"]
            record["matching_plan_index"] = best["plan_index"] if best["comparison"]["all_inside"]["passes_quantization_bound"] else None
            record["matching_positive_plan_index"] = best["plan_index"] if best["comparison"]["positive_inside"]["passes_quantization_bound"] else None
            masked = best["comparison"]["source_support_at_reference_voxels"]["all_voxels"]
            record["matching_supported_plan_index"] = best["plan_index"] if masked["passes_quantization_bound"] else None
            chosen = discover_plans(case)[best["plan_index"]]
            fine = resample_to_tissue_grid(read_dose(chosen, case), chosen, case)
            record["fine_tissue_export_grid_comparison"] = compare(fine, reference)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2)+"\n", encoding="utf-8")


if __name__ == "__main__":
    try:
        main()
    except Exception:
        print(json.dumps({'status': 'failed', 'error': 'Dose validation failed; source details omitted'}), flush=True)
        raise SystemExit(1)
