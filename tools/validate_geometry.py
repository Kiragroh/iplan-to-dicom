"""Read-only, identifier-free image/structure validation against native exports.

Source identities are used only in memory for pairing. The JSON deliberately
contains no paths, UIDs, names, dates of acquisition, or structure names.
The reference export is a validation input, never a reconstruction input.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np
import pydicom
from skimage.measure import points_in_poly

from iplan_dicom.images import RAW_TO_PATIENT, load_image_volume, load_slice_set
from iplan_dicom.plan import discover_plans, objects, read_plan
from iplan_dicom.service import load_headers, selected_slice_set, load_rois
from iplan_dicom.xbrain import FormatError


def safe_error(exc: Exception) -> dict:
    # FormatError messages are authored by the local package and omit identities.
    return {"type": type(exc).__name__, "reason": str(exc) if isinstance(exc, FormatError) else "Unexpected validation error; inspect locally without publishing source values"}


def reference_cases(root: Path) -> dict:
    result = {}
    for directory in sorted(path for path in root.iterdir() if path.is_dir()):
        files = list(directory.rglob("*.dcm"))
        if not files:
            raise FormatError("Reference case has no DICOM files")
        def read(path):
            try:
                return pydicom.dcmread(path, stop_before_pixels=True)
            except Exception:
                raise FormatError("Reference DICOM file cannot be read") from None
        with ThreadPoolExecutor(max_workers=6) as pool:
            datasets = list(pool.map(read, files))
        identities = {str(ds.get("PatientID", "")) for ds in datasets}
        if len(identities) != 1 or not next(iter(identities)):
            raise FormatError("Reference case patient identity is inconsistent")
        identity = next(iter(identities))
        if identity in result:
            raise FormatError("Reference patient pairing is ambiguous")
        result[identity] = datasets
    return result


def angle_error(a, b) -> float:
    return abs((float(a) - float(b) + 180) % 360 - 180)


def plane_sampling_summary(affines, width: int, height: int, image_indices) -> dict:
    centers = np.asarray([(matrix @ np.array([(width - 1) / 2, (height - 1) / 2, 0, 1]))[:3] for matrix in affines])
    distances = np.linalg.norm(np.diff(centers, axis=0), axis=1)
    jumps = np.diff(np.asarray(image_indices, int))
    return {"plane_center_distance_minimum_mm": float(distances.min()),
            "plane_center_distance_median_mm": float(np.median(distances)),
            "plane_center_distance_maximum_mm": float(distances.max()),
            "maximum_source_image_index_step": int(jumps.max()),
            "source_image_index_jumps_greater_than_one": int(np.sum(jumps > 1)),
            "omitted_internal_source_images": int(np.maximum(jumps - 1, 0).sum())}


def compare_plan(plan, reference) -> dict:
    """Identify reference-compatible plans using MU and all present CP geometry.

    Predeclared tolerances only cover serialized decimal rounding. The native
    reference lacks leaf positions, so this does not establish MLC equivalence.
    """
    groups = list(reference.FractionGroupSequence)
    if len(groups) != 1:
        raise FormatError("Reference plan has multiple fraction groups")
    beams = list(reference.BeamSequence)
    groups_ref = {int(item.ReferencedBeamNumber): float(item.BeamMeterset)
                  for item in groups[0].ReferencedBeamSequence}
    out = {"beam_count_equal": len(plan.beams) == len(beams),
           "fractions_equal": plan.fractions == int(groups[0].NumberOfFractionsPlanned),
           "control_point_counts_equal": True, "max_mu_difference": 0.0,
           "max_angle_difference_degrees": 0.0, "max_weight_difference": 0.0,
           "max_isocenter_difference_mm": 0.0}
    if not out["beam_count_equal"]:
        out["matches_reference"] = False
        return out
    for source_beam, beam in zip(plan.beams, beams):
        out["max_mu_difference"] = max(out["max_mu_difference"], abs(source_beam.mu - groups_ref[int(beam.BeamNumber)]))
        cps = list(beam.ControlPointSequence)
        if len(cps) != len(source_beam.points):
            out["control_point_counts_equal"] = False
            continue
        prior = {}
        for source_cp, cp in zip(source_beam.points, cps):
            for field in ("GantryAngle", "PatientSupportAngle", "BeamLimitingDeviceAngle", "IsocenterPosition"):
                if hasattr(cp, field):
                    prior[field] = getattr(cp, field)
            for field, value in (("GantryAngle", source_cp.gantry), ("PatientSupportAngle", source_cp.table),
                                 ("BeamLimitingDeviceAngle", source_cp.collimator)):
                if field not in prior:
                    raise FormatError("Reference plan control-point angle is missing")
                out["max_angle_difference_degrees"] = max(out["max_angle_difference_degrees"], angle_error(value, prior[field]))
            weight = float(cp.CumulativeMetersetWeight) / float(beam.FinalCumulativeMetersetWeight)
            out["max_weight_difference"] = max(out["max_weight_difference"], abs(source_cp.weight - weight))
            if "IsocenterPosition" not in prior:
                raise FormatError("Reference plan isocenter is missing")
            iso = (RAW_TO_PATIENT @ source_beam.isocenter)[:3]
            out["max_isocenter_difference_mm"] = max(out["max_isocenter_difference_mm"], float(np.max(np.abs(iso - np.asarray(prior["IsocenterPosition"], float)))))
    out["matches_reference"] = bool(out["beam_count_equal"] and out["fractions_equal"] and out["control_point_counts_equal"]
                                    and out["max_mu_difference"] <= 0.001
                                    and out["max_angle_difference_degrees"] <= 0.0001
                                    and out["max_weight_difference"] <= 0.00001
                                    and out["max_isocenter_difference_mm"] <= 0.001)
    return out


def polygon_statistics(planes: dict[int, list[np.ndarray]], basis: dict[int, np.ndarray], thickness: float) -> dict:
    """Area-weighted centroid and explicit volume proxy, with nested-ring XOR."""
    area_total = 0.0
    moment = np.zeros(3)
    lower, upper = np.full(3, np.inf), np.full(3, -np.inf)
    point_count = 0
    for index, rings in planes.items():
        axes = basis[index]
        ux = axes[:3, 0] / np.linalg.norm(axes[:3, 0])
        uy = axes[:3, 1] / np.linalg.norm(axes[:3, 1])
        normal = np.cross(ux, uy)
        projected = [np.column_stack((points @ ux, points @ uy)) for points in rings]
        for ring_index, (points, polygon) in enumerate(zip(rings, projected)):
            if len(points) < 3 or not np.isfinite(points).all():
                raise FormatError("Invalid contour coordinates")
            normal_positions = points @ normal
            if np.ptp(normal_positions) > 0.002:
                raise FormatError("Contour is not planar in its referenced image")
            x, y = polygon.T
            xn, yn = np.roll(x, -1), np.roll(y, -1)
            cross = x * yn - xn * y
            signed_area = cross.sum() / 2
            if abs(signed_area) < 1e-10:
                raise FormatError("Degenerate contour has zero area")
            nesting = sum(bool(points_in_poly(polygon[:1], other)[0]) for other_index, other in enumerate(projected) if other_index != ring_index)
            signed_weight = abs(signed_area) * (-1 if nesting % 2 else 1)
            centroid2d = np.array([((x + xn) * cross).sum(), ((y + yn) * cross).sum()]) / (6 * signed_area)
            centroid = centroid2d[0] * ux + centroid2d[1] * uy + normal_positions.mean() * normal
            area_total += signed_weight
            moment += signed_weight * centroid
            lower = np.minimum(lower, points.min(axis=0))
            upper = np.maximum(upper, points.max(axis=0))
            point_count += len(points)
    if area_total <= 0:
        raise FormatError("Contour set has no positive XOR area")
    return {"centroid": moment / area_total, "lower": lower, "upper": upper,
            "volume_proxy": area_total * thickness, "points": point_count, "planes": len(planes)}


def structure_comparisons(rois, metrics, image_cache, reference_rs, reference_ct) -> dict:
    names = {}
    for entry in reference_rs.StructureSetROISequence:
        key = str(entry.ROIName).strip().casefold()
        if key in names:
            raise FormatError("Reference structure names are ambiguous")
        names[key] = int(entry.ROINumber)
    contour_by_number = {int(item.ReferencedROINumber): item for item in reference_rs.ROIContourSequence}
    ct_by_uid = {str(ds.SOPInstanceUID): ds for ds in reference_ct}
    ct_positions = sorted(float(ds.ImagePositionPatient[2]) for ds in reference_ct)
    reference_thickness = float(np.median(np.diff(ct_positions)))
    rows, matched_names = [], set()
    for index, (roi, metric) in enumerate(zip(rois, metrics), 1):
        entry = {"roi_index": index, **metric}
        image = image_cache[roi.source_grid]
        source_thickness = float(np.median(np.linalg.norm(np.diff(np.array([image.slice_affine(i)[:3, 3] for i in range(len(image.pixels))]), axis=0), axis=1)))
        entry.update(source_grid_modality=image.modality, source_grid_planes=len(image.pixels))
        key = roi.name.strip().casefold()
        if key not in names:
            entry["status"] = "source_name_unmatched"
            rows.append(entry)
            continue
        matched_names.add(key)
        reference = contour_by_number.get(names[key])
        if reference is None or not getattr(reference, "ContourSequence", None):
            entry["status"] = "reference_has_no_contours"
            rows.append(entry)
            continue
        if not roi.contours:
            entry["status"] = "source_has_no_contours"
            rows.append(entry)
            continue
        source_planes, source_basis = {}, {}
        for plane_index, points in roi.contours:
            source_planes.setdefault(plane_index, []).append(points)
            source_basis[plane_index] = image.slice_affine(plane_index)
        ref_planes, ref_basis, plane_ids = {}, {}, {}
        for contour in reference.ContourSequence:
            if str(contour.ContourGeometricType) not in ("CLOSED_PLANAR", "CLOSEDPLANAR_XOR"):
                raise FormatError("Unsupported reference contour type")
            points = np.asarray(contour.ContourData, float).reshape(-1, 3)
            links = list(getattr(contour, "ContourImageSequence", []))
            if len(links) != 1 or str(links[0].ReferencedSOPInstanceUID) not in ct_by_uid:
                raise FormatError("Reference contour image link is unresolved")
            uid = str(links[0].ReferencedSOPInstanceUID)
            plane_index = plane_ids.setdefault(uid, len(plane_ids))
            header = ct_by_uid[uid]
            iop = np.asarray(header.ImageOrientationPatient, float)
            axes = np.eye(4)
            axes[:3, 0], axes[:3, 1] = iop[:3], iop[3:]
            ref_planes.setdefault(plane_index, []).append(points)
            ref_basis[plane_index] = axes
        source_stats = polygon_statistics(source_planes, source_basis, source_thickness)
        ref_stats = polygon_statistics(ref_planes, ref_basis, reference_thickness)
        entry.update(status="compared", centroid_distance_mm=float(np.linalg.norm(source_stats["centroid"] - ref_stats["centroid"])),
                     bounding_extrema_max_difference_mm=float(max(np.max(np.abs(source_stats["lower"] - ref_stats["lower"])), np.max(np.abs(source_stats["upper"] - ref_stats["upper"])))),
                     source_volume_proxy_mm3=source_stats["volume_proxy"], reference_volume_proxy_mm3=ref_stats["volume_proxy"],
                     volume_proxy_relative_difference_percent=100 * (source_stats["volume_proxy"] / ref_stats["volume_proxy"] - 1),
                     source_contour_planes=source_stats["planes"], reference_contour_planes=ref_stats["planes"])
        rows.append(entry)
    compared = [row for row in rows if row["status"] == "compared"]
    return {"source_roi_count": len(rois), "reference_roi_count": len(names), "compared_count": len(compared),
            "unmatched_reference_count": len(set(names) - matched_names),
            "maximum_centroid_distance_mm": max((r["centroid_distance_mm"] for r in compared), default=None),
            "maximum_bounding_extrema_difference_mm": max((r["bounding_extrema_max_difference_mm"] for r in compared), default=None),
            "maximum_absolute_volume_proxy_difference_percent": max((abs(r["volume_proxy_relative_difference_percent"]) for r in compared), default=None),
            "rois": rows}


def validate_case(case: Path, reference: list, index: int, progress) -> dict:
    result = {"case_index": index, "held_out_until_this_validation": index == 7, "status": "running", "stage": "headers"}
    try:
        headers = load_headers(case)
        plan_refs = [ds for ds in reference if ds.Modality == "RTPLAN"]
        structures = [ds for ds in reference if ds.Modality == "RTSTRUCT"]
        ref_ct = [ds for ds in reference if ds.Modality == "CT"]
        if len(plan_refs) != 1 or len(structures) != 1 or len(ref_ct) < 2:
            raise FormatError("Reference must contain one plan, one structure set and a CT series")
        result["stage"] = "plan_matching"
        candidates, supported = [], {}
        for plan_index, directory in enumerate(discover_plans(case)):
            row = {"plan_index": plan_index}
            try:
                plan = read_plan(directory, case)
                row.update(compare_plan(plan, plan_refs[0]))
                supported[plan_index] = (plan, (directory / "rtplan.xbrain").stat().st_mtime_ns)
            except Exception as exc:
                row.update(matches_reference=False, error=safe_error(exc))
            candidates.append(row)
        result["plan_candidates"] = candidates
        matching = [row["plan_index"] for row in candidates if row.get("matches_reference")]
        result["matching_plan_indices"] = matching
        result["plan_identity_ambiguous"] = len(matching) > 1
        if not matching:
            raise FormatError("No source plan matches reference MU and control-point geometry")
        selected = max(matching, key=lambda candidate: (supported[candidate][1], candidate))
        result["selected_plan_index"] = selected
        result["selection_rule"] = "latest_source_pointer_mtime_then_highest_index_among_geometry_matches; ambiguity remains explicit"
        plan = supported[selected][0]
        result["stage"] = "selected_CT"
        progress(index, "selected_CT")
        slices, directory = selected_slice_set(plan, case)
        image = load_slice_set(slices, directory, case, headers)
        axes = np.array([image.slice_affine(i)[:3, :2] for i in range(len(image.pixels))])
        orthogonality = np.abs(np.sum(axes[:, :, 0] * axes[:, :, 1], axis=1))
        result["image"] = {"ct_shape": list(image.pixels.shape), "ct_pixel_count": int(image.pixels.size),
                           "ct_sample_minimum": int(image.pixels.min()), "ct_sample_maximum": int(image.pixels.max()),
                           "ct_dtype": str(image.pixels.dtype), "individual_plane_affines": len(image.slice_affines),
                           "uniform_grid": image.uniform_grid, "maximum_plane_axis_dot_product": float(orthogonality.max()),
                           "source_HU_slope": image.rescale_slope, "source_HU_intercept": image.rescale_intercept,
                           "reference_ct_slices": len(ref_ct)}
        result["image"].update(plane_sampling_summary(image.slice_affines, image.pixels.shape[2], image.pixels.shape[1],
                                                       [int(node.get('imageIndex')) for node in slices.findall('Slices/Slice')]))
        result["stage"] = "all_structure_masks"
        progress(index, "all_structure_masks")
        image_cache = {directory: image}
        rois, metrics = load_rois(plan, case, image, headers=headers, image_cache=image_cache)
        result["stage"] = "structure_reference_comparison"
        result["structures"] = structure_comparisons(rois, metrics, image_cache, structures[0], ref_ct)
        result["registered_ROI_grids"] = [{"grid_index": i + 1, "modality": grid.modality, "shape": list(grid.pixels.shape),
                                            "uniform_grid": grid.uniform_grid, "primary_CT_grid": grid is image}
                                           for i, grid in enumerate(image_cache.values())]
        del rois, image, image_cache
        result["stage"] = "original_MR"
        progress(index, "original_MR")
        result["original_MR"] = []
        for source in sorted(case.rglob("imageseries.xbrain")):
            for series in objects(source, "ImageSeries"):
                first = series.find("Image")
                if first is None or first.get("modality") != "MR":
                    continue
                mr = load_image_volume(series, source.parent, headers)
                result["original_MR"].append({"series_index": len(result["original_MR"]) + 1, "shape": list(mr.pixels.shape),
                                             "pixel_count": int(mr.pixels.size), "dtype": str(mr.pixels.dtype),
                                             "sample_minimum": int(mr.pixels.min()), "sample_maximum": int(mr.pixels.max()),
                                             "header_geometry_checked": True, "fusion_applied": False})
                del mr
        result["status"] = "completed_with_explicit_comparison_limits"
        result["stage"] = "complete"
    except Exception as exc:
        result["status"] = "failed"
        result["error"] = safe_error(exc)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-root", required=True, type=Path)
    parser.add_argument("--reference-root", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--case-indices", help="Optional comma-separated 1-based indices; default all")
    args = parser.parse_args()
    source, reference_root, output = args.source_root.resolve(), args.reference_root.resolve(), args.output.resolve()
    if source == output or source in output.parents or reference_root == output or reference_root in output.parents:
        raise FormatError("Validation report must remain outside both source roots")
    cases = sorted(path for path in source.iterdir() if path.is_dir())
    selected = set(map(int, args.case_indices.split(","))) if args.case_indices else set(range(1, len(cases) + 1))
    if not selected <= set(range(1, len(cases) + 1)):
        raise FormatError("Requested validation case index does not exist")
    report = {"schema_version": 1, "validation_started_utc": datetime.now(timezone.utc).isoformat(),
              "status": "running", "source_case_count": len(cases), "selected_case_indices": sorted(selected),
              "read_only_sources": True, "reference_used_only_for_validation": True,
              "plan_match_tolerances": {"MU": .001, "angle_degrees": .0001, "meterset_weight": .00001, "isocenter_mm": .001},
              "comparison_limits": ["Centroids are polygon-area-weighted, not a clinical acceptance criterion.",
                                    "Volume proxies use median plane spacing; registered nonparallel planes and native resampling differ.",
                                    "No anatomy-dependent tolerance is fitted from these results.",
                                    "Native reference lacks MLC positions; plan matching cannot establish leaf equivalence.",
                                    "ROI-owning MR grids use the selected plan's native fusion; optional original MR remains in its scanner frame.",
                                    "Unsupported source objects and comparison failures are reported explicitly."], "cases": []}
    def save():
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(report, indent=2, allow_nan=False), encoding="utf-8")
    save()
    references = reference_cases(reference_root)
    def progress(index, stage):
        print(json.dumps({"case_index": index, "stage": stage}), flush=True)
    for index, case in enumerate(cases, 1):
        if index not in selected:
            continue
        progress(index, "pairing")
        try:
            header_map = load_headers(case)
            identities = {str(ds.get("PatientID", "")) for ds in header_map.values() if ds.get("Modality") == "CT"}
            if len(identities) != 1 or next(iter(identities)) not in references:
                raise FormatError("Source-to-reference case pairing is unresolved")
            reference = references[next(iter(identities))]
            del header_map
            result = validate_case(case, reference, index, progress)
        except Exception as exc:
            result = {"case_index": index, "status": "failed", "stage": "pairing", "error": safe_error(exc)}
        report["cases"].append(result)
        save()
        progress(index, result["status"])
    report["status"] = "completed_with_failures" if any(c["status"] == "failed" for c in report["cases"]) else "completed"
    report["validation_finished_utc"] = datetime.now(timezone.utc).isoformat()
    save()
    return 1 if report["status"] == "completed_with_failures" else 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(json.dumps({"status": "failed", "error": safe_error(exc)}), flush=True)
        raise SystemExit(1)
