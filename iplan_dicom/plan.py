"""Read native conformal beams/arcs and preserve them in research RTPLANs.

Native references are resolved only inside the selected case. No reference
DICOM is used to construct a plan. Source machine leaf geometry is required.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path, PureWindowsPath
import xml.etree.ElementTree as ET

import numpy as np
from pydicom.dataset import Dataset
from pydicom.sequence import Sequence

from .xbrain import FormatError, numbers, read_document, scalar


def content(root: ET.Element) -> ET.Element:
    versions = root.findall("./Content/ContentVersion")
    if not versions:
        raise FormatError("Object has no content schema")
    # These are alternative serialization schemas, not historical revisions.
    try:
        supported = [n for n in versions if int(n.attrib["value"]) in (1, 2, 3, 4)]
        if not supported:
            raise FormatError("Unsupported content schema version")
        return max(supported, key=lambda n: int(n.attrib["value"]))
    except (KeyError, ValueError):
        raise FormatError("Invalid content schema version") from None


def text(node: ET.Element, path: str) -> str:
    child = node.find(path)
    if child is None:
        raise FormatError("Required text field is missing")
    for kind in ("_ASCII_", "_UNICODE_"):
        value = child.find(".//" + kind)
        if value is not None and value.text is not None:
            return value.text
    if child.text and child.text.strip():
        return child.text.strip()
    raise FormatError("Required text field is empty")


def vector(node: ET.Element, path: str) -> np.ndarray:
    return np.asarray(numbers(text(node, path), counted=True), dtype=float)


def resolve_directory(case_dir: Path, containing_dir: Path, relative: str) -> Path:
    win = PureWindowsPath(relative or ".")
    if win.drive or win.root:
        raise FormatError("Absolute object reference is unsupported")
    resolved = containing_dir.joinpath(*win.parts).resolve()
    try:
        resolved.relative_to(case_dir.resolve())
    except ValueError:
        raise FormatError("Object reference escapes the selected case") from None
    return resolved


def objects(path: Path, tag: str) -> list[ET.Element]:
    if not path.exists():
        return []
    root = read_document(path)
    container = root.find("Content")
    if container is not None and container.get("count") == "0":
        return []
    return list(content(root).iter(tag))


def referenced_object(node: ET.Element, tag: str, folder: Path, case_dir: Path) -> ET.Element:
    refs = [r for r in node.findall("Reference") if r.get("contentType") == tag]
    if len(refs) != 1:
        raise FormatError("Object must have one required typed reference")
    ref = refs[0]
    directory = resolve_directory(case_dir, folder, ref.get("directory", "."))
    found = [x for x in objects(directory / (tag.lower() + ".xbrain"), tag)
             if x.get("uId") == ref.get("contentUid")]
    if len(found) != 1:
        raise FormatError("Required object reference is unresolved or ambiguous")
    return found[0]


def discover_plans(case_dir: Path) -> list[Path]:
    return [p.parent for p in sorted(Path(case_dir).rglob("rtplan.xbrain"))
            if objects(p, "RTPlan")]


def treatment_elements(plan_dir: Path) -> list[tuple[Path, ET.Element]]:
    """Inventory every current treatment collection, including unknown types."""
    result = []
    expected = {"confarcelement.xbrain": "ConfArcElement", "confbeamelement.xbrain": "ConfBeamElement"}
    for path in sorted(plan_dir.glob("*element.xbrain")):
        root = read_document(path)
        container = root.find("Content")
        if container is None:
            raise FormatError("Treatment collection has no content container")
        if not container.findall("ContentVersion"):
            if container.get("count") == "0" and not list(container):
                continue
            raise FormatError("Populated treatment collection has no content schema")
        schema = content(root)
        elements = [item for child in schema for item in (list(child) if child.tag == "Multiple" else [child])]
        if not elements:
            continue
        wanted = expected.get(path.name.lower())
        if wanted is None or any(item.tag != wanted for item in elements):
            raise FormatError("Unsupported treatment element collection or type")
        if container.get("count") == "0":
            raise FormatError("Treatment collection content contradicts its count")
        result.extend((path, item) for item in elements)
    return result


def local_component(element: ET.Element, components: list[ET.Element],
                    plan_dir: Path, case_dir: Path) -> ET.Element:
    refs = [r for r in element.findall("Reference") if r.get("contentType") == "LINACIsocComponent"]
    if len(refs) != 1:
        raise FormatError("Beam must have exactly one component reference")
    directory = resolve_directory(case_dir, plan_dir, refs[0].get("directory", "."))
    if directory != plan_dir.resolve():
        raise FormatError("Beam component reference targets another plan directory")
    linked = referenced_object(element, "LINACIsocComponent", plan_dir, case_dir)
    matches = [c for c in components if c.get("uId") == linked.get("uId")]
    if len(matches) != 1 or ET.tostring(matches[0]) != ET.tostring(linked):
        raise FormatError("Beam component reference differs from selected local component")
    return matches[0]


def blocking_devices(profile: ET.Element, definitions: dict[tuple[str, str], ET.Element]) -> tuple[ET.Element, ET.Element]:
    """Resolve every declared blocking device; the supported set is jaws + MLC."""
    devices: dict[str, ET.Element] = {}
    for reference in profile.findall("BlockingDeviceId"):
        identifier = text(reference, ".")
        matches = [obj for (tag, key), obj in definitions.items() if key == identifier]
        if len(matches) != 1:
            raise FormatError("Blocking device reference is unresolved or ambiguous")
        device = matches[0]
        if device.tag not in ("JawsParameters", "MultileafParameters"):
            raise FormatError("Unsupported blocking device or accessory")
        if device.tag in devices:
            raise FormatError("Duplicate blocking device type")
        devices[device.tag] = device
    if set(devices) != {"JawsParameters", "MultileafParameters"}:
        raise FormatError("Exactly one jaw and one MLC device are required")
    jaws = devices["JawsParameters"]
    if (jaws.get("availableX") != "true" or jaws.get("availableY") != "true" or
            jaws.get("symmetricX") != "false" or jaws.get("symmetricY") != "false" or
            scalar(jaws, "Orientation") != 0):
        raise FormatError("Jaw availability, symmetry or orientation has not been verified")
    return jaws, devices["MultileafParameters"]


@dataclass
class Machine:
    name: str
    manufacturer: str
    model: str
    energy: float
    sad: float
    leaf_boundaries: np.ndarray


@dataclass
class ControlPoint:
    gantry: float
    table: float
    collimator: float
    weight: float
    jaws: np.ndarray
    leaves: np.ndarray


@dataclass
class Beam:
    number: int
    mu: float
    source_single_dose: float
    dose_rate: float
    direction: str
    machine: Machine
    isocenter: np.ndarray
    points: list[ControlPoint]
    source: ET.Element = field(repr=False)
    source_path: Path = field(repr=False)
    name: str = field(default="", repr=False)


@dataclass
class Plan:
    beams: list[Beam]
    fractions: int
    prescription_per_fraction: float
    source: ET.Element = field(repr=False)
    source_directory: Path = field(repr=False)
    warnings: tuple[str, ...] = (
        "MLC positions are preserved from the native machine and field records; "
        "the available manufacturer export omits MLC positions and cannot validate them.",
        "Native SingleDose is a planning weight, not an independently verified DICOM BeamDose.",
    )


def read_machine(component: ET.Element, plan_dir: Path) -> Machine:
    definitions: dict[tuple[str, str], ET.Element] = {}
    # Each plan stores its own machine snapshot. The same identifier can have
    # revised definitions in another plan, so a case-wide index is incorrect.
    for path in sorted(plan_dir.rglob("*.xmach")):
        obj = read_document(path)
        key = (obj.tag, obj.get("identifier", ""))
        previous = definitions.get(key)
        if previous is not None and ET.tostring(previous) != ET.tostring(obj):
            raise FormatError("Conflicting machine definitions")
        definitions[key] = obj

    def lookup(tag: str, identifier: str) -> ET.Element:
        obj = definitions.get((tag, identifier))
        if obj is None:
            raise FormatError("Required machine definition is unavailable")
        return obj

    profile = lookup("MachineProfileFileable", component.get("machineConfiguration", ""))
    linac = lookup("LINACParameters", profile.get("linacId", ""))
    if linac.get("convention") != "Generic.IEC1217":
        raise FormatError("Machine angle convention has not been verified")
    _, mlc = blocking_devices(profile, definitions)
    geom = lookup("MultileafGeometry", mlc.get("multileafGeometry", ""))
    if geom.get("isocentric") != "true" or scalar(mlc, "Orientation") != 0:
        raise FormatError("Non-isocentric or rotated MLC geometry is unsupported")
    if mlc.get("swapBanksOnExport") != "false":
        raise FormatError("Machine requires an unverified MLC bank transformation")
    begins, ends = vector(geom, "LeafBeginY"), vector(geom, "LeafEndY")
    n = int(geom.get("numberLeafPairs", "0"))
    if len(begins) != n or len(ends) != n or n < 1:
        raise FormatError("MLC geometry count mismatch")
    if np.any(ends <= begins) or not np.allclose(ends[:-1], begins[1:], atol=1e-9, rtol=0):
        raise FormatError("MLC leaf boundaries are not contiguous and ascending")
    sad = scalar(linac, "SID")
    if not np.isclose(scalar(geom, "DefinitionSID"), sad, atol=1e-9):
        raise FormatError("MLC definition plane differs from source-axis distance")
    name = linac.get("exportId", "")
    if not name or len(name) > 16:
        raise FormatError("Machine export identifier is missing or too long")
    return Machine(name, text(linac, "Manufacturer"), text(linac, "Model"),
                   scalar(profile, "NominalEnergy"), sad, np.r_[begins, ends[-1]])


def read_points(element: ET.Element, leaf_count: int) -> tuple[list[ControlPoint], str]:
    states = element.findall("ConformalField")
    is_arc = element.tag == "ConfArcElement"
    if (is_arc and len(states) < 2) or (not is_arc and len(states) != 1):
        raise FormatError("Unsupported conformal field control-point count")
    if is_arc and element.get("arcType") not in ("dynamic", "conformal"):
        raise FormatError("Unknown conformal arc type")
    raw_direction = element.get("desiredGantryRotationDirection")
    if is_arc and raw_direction not in ("clockwise", "counterclockwise"):
        raise FormatError("Unknown arc direction")
    reverse = is_arc and raw_direction == "counterclockwise"
    points = []
    for state in states:
        jaws = np.asarray([scalar(state, "JawsState/" + k) for k in ("X1", "X2", "Y1", "Y2")])
        left, right = vector(state, "MultileafState/X1"), vector(state, "MultileafState/X2")
        if len(left) != leaf_count or len(right) != leaf_count:
            raise FormatError("Field MLC leaf count differs from the machine")
        if np.any(left > right) or jaws[0] > jaws[1] or jaws[2] > jaws[3]:
            raise FormatError("Reversed leaf or jaw bank positions")
        native_weight = scalar(state, "CumMeterSetWeight")
        if not is_arc and native_weight != -1:
            raise FormatError("Static field requires the verified native meterset sentinel")
        weight = native_weight if is_arc else 0.0
        points.append(ControlPoint(scalar(state, "GantryAngle") % 360,
                                   scalar(state, "TableAngle") % 360,
                                   scalar(state, "CollimatorAngle") % 360,
                                   1 - weight if reverse else weight, jaws, np.r_[left, right]))
    if reverse:
        points.reverse()
    if any(not np.isclose(point.table, points[0].table) or
           not np.isclose(point.collimator, points[0].collimator) for point in points[1:]):
        raise FormatError("Moving couch or collimator has not been verified")
    if is_arc and element.get("arcType") == "conformal":
        if any(not np.array_equal(point.leaves, points[0].leaves) for point in points[1:]):
            raise FormatError("Conformal arc has unexpectedly varying leaf positions")
    if not is_arc:
        first = points[0]
        points.append(ControlPoint(first.gantry, first.table, first.collimator, 1.0,
                                   first.jaws.copy(), first.leaves.copy()))
    weights = np.asarray([p.weight for p in points])
    if not np.isclose(weights[0], 0) or not np.isclose(weights[-1], 1) or np.any(np.diff(weights) < 0):
        raise FormatError("Cumulative meterset weights are invalid")
    return points, ("CC" if reverse else "CW") if is_arc else "NONE"


def read_plan(plan_dir: Path, case_dir: Path) -> Plan:
    plan_dir, case_dir = Path(plan_dir), Path(case_dir)
    plans = objects(plan_dir / "rtplan.xbrain", "RTPlan")
    if len(plans) != 1:
        raise FormatError("Selected folder must contain one RTPlan")
    plan = plans[0]
    if plan.get("patientOrientation") != "Supine" or plan.get("headFeetOrientation") != "HeadFirst":
        raise FormatError("Patient orientation has not been verified")
    components = sorted(objects(plan_dir / "linacisoccomponent.xbrain", "LINACIsocComponent"),
                        key=lambda x: int(x.get("componentSortId", "-1")))
    sources = treatment_elements(plan_dir)
    assigned = {id(element): local_component(element, components, plan_dir, case_dir)
                for _, element in sources}
    beams = []
    used = set()
    for component in components:
        if component.get("useDynamicTracking", "false") != "false":
            raise FormatError("Dynamic tracking has not been verified")
        if any(scalar(component, key) != 0 for key in
               ("TableVertical", "TableLongitudinal", "TableLateral")):
            raise FormatError("Nonzero native table translations have not been verified")
        linked_plan = referenced_object(component, "RTPlan", plan_dir, case_dir)
        if linked_plan.get("uId") != plan.get("uId"):
            raise FormatError("Component belongs to another plan")
        machine = read_machine(component, plan_dir)
        iso_obj = referenced_object(component, "Isocenter", plan_dir, case_dir)
        iso = vector(iso_obj, "NamedPoint/Point/Point_3D")
        if iso.shape != (4,) or not np.isclose(iso[3], 1):
            raise FormatError("Invalid homogeneous isocenter point")
        selected = [(p, x) for p, x in sources if assigned[id(x)] is component]
        selected.sort(key=lambda item: int(item[1].get("sortId", "-1")))
        for path, element in selected:
            if scalar(element, "TableEccentricAngle") != 0:
                raise FormatError("Nonzero table eccentric rotation has not been verified")
            if id(element) in used:
                raise FormatError("Beam is associated with multiple components")
            used.add(id(element))
            points, direction = read_points(element, len(machine.leaf_boundaries) - 1)
            mu, rate = scalar(element, "SingleMU"), scalar(component, "DoseRate") * 60
            if mu <= 0 or rate <= 0:
                raise FormatError("Beam meterset and dose rate must be positive")
            beams.append(Beam(len(beams) + 1, mu, scalar(element, "SingleDose"), rate,
                              direction, machine, iso, points, element, path, text(element, "Name")))
    if not beams or len(used) != len(sources):
        raise FormatError("Plan has no beams or unresolved beam-component references")
    coords = objects(plan_dir / "rtobjectcoordinator.xbrain", "RTObjectCoordinator")
    fractions = {(int(r.get("numberOfFractions", "0")), scalar(r, "DosePerFraction"))
                 for co in coords for r in co.findall("RTObject") if r.get("hasFractionation") == "true"}
    if len(fractions) != 1:
        raise FormatError("Fractionation is missing or differs across targets")
    count, dose = fractions.pop()
    if count < 1 or dose <= 0:
        raise FormatError("Invalid fractionation")
    return Plan(beams, count, dose, plan, plan_dir)


def populate_rtplan_dataset(ds: Dataset, plan: Plan, raw_to_patient: np.ndarray) -> Dataset:
    matrix = np.asarray(raw_to_patient, dtype=float)
    if matrix.shape != (4, 4) or not np.all(np.isfinite(matrix)) or not np.allclose(matrix[3], [0, 0, 0, 1]):
        raise FormatError("Invalid isocenter coordinate transform")
    ds.RTPlanLabel = "IPLAN_RESEARCH"
    ds.RTPlanName = "Native iPlan conversion"
    ds.RTPlanDescription = "Research reconstruction; not clinically released"
    ds.RTPlanGeometry = "PATIENT"
    ds.ApprovalStatus = "UNAPPROVED"
    ds.BeamSequence = Sequence()
    group = Dataset()
    group.FractionGroupNumber = 1
    group.NumberOfFractionsPlanned = plan.fractions
    group.NumberOfBeams = len(plan.beams)
    group.NumberOfBrachyApplicationSetups = 0
    group.ReferencedBeamSequence = Sequence()
    ds.FractionGroupSequence = Sequence([group])
    setup = Dataset()
    setup.PatientSetupNumber = 1
    setup.PatientPosition = "HFS"
    ds.PatientSetupSequence = Sequence([setup])
    for beam in plan.beams:
        b = Dataset()
        b.BeamNumber = beam.number
        b.BeamName = beam.name[:64] if beam.name else f"B{beam.number:02d}"
        b.BeamType = "STATIC" if beam.direction == "NONE" else "DYNAMIC"
        b.TreatmentDeliveryType = "TREATMENT"
        b.RadiationType = "PHOTON"
        b.PrimaryDosimeterUnit = "MU"
        b.TreatmentMachineName = beam.machine.name
        b.Manufacturer = beam.machine.manufacturer
        b.ManufacturerModelName = beam.machine.model
        b.SourceAxisDistance = beam.machine.sad
        b.NumberOfWedges = b.NumberOfCompensators = b.NumberOfBoli = b.NumberOfBlocks = 0
        b.ReferencedPatientSetupNumber = 1
        b.FinalCumulativeMetersetWeight = 1.0
        b.NumberOfControlPoints = len(beam.points)
        b.BeamLimitingDeviceSequence = Sequence()
        for typ, count in (("ASYMX", 1), ("ASYMY", 1), ("MLCX", len(beam.machine.leaf_boundaries) - 1)):
            device = Dataset()
            device.RTBeamLimitingDeviceType = typ
            device.NumberOfLeafJawPairs = count
            if typ == "MLCX":
                device.LeafPositionBoundaries = beam.machine.leaf_boundaries.tolist()
            b.BeamLimitingDeviceSequence.append(device)
        b.ControlPointSequence = Sequence()
        for index, point in enumerate(beam.points):
            cp = Dataset()
            cp.ControlPointIndex = index
            cp.CumulativeMetersetWeight = point.weight
            cp.NominalBeamEnergy = beam.machine.energy
            cp.DoseRateSet = beam.dose_rate
            cp.GantryAngle = point.gantry
            cp.GantryRotationDirection = beam.direction if index < len(beam.points) - 1 else "NONE"
            cp.PatientSupportAngle = point.table
            cp.PatientSupportRotationDirection = "NONE"
            cp.TableTopEccentricAngle = 0.0
            cp.TableTopEccentricRotationDirection = "NONE"
            # Native zeros are also used as absent positioning defaults. Do not
            # turn them into a claim about an actual treatment table position.
            cp.TableTopVerticalPosition = ""
            cp.TableTopLongitudinalPosition = ""
            cp.TableTopLateralPosition = ""
            cp.BeamLimitingDeviceAngle = point.collimator
            cp.BeamLimitingDeviceRotationDirection = "NONE"
            cp.IsocenterPosition = (matrix @ beam.isocenter)[:3].tolist()
            cp.BeamLimitingDevicePositionSequence = Sequence()
            for typ, positions in (("ASYMX", point.jaws[:2]), ("ASYMY", point.jaws[2:]), ("MLCX", point.leaves)):
                device = Dataset()
                device.RTBeamLimitingDeviceType = typ
                device.LeafJawPositions = positions.tolist()
                cp.BeamLimitingDevicePositionSequence.append(device)
            b.ControlPointSequence.append(cp)
        ds.BeamSequence.append(b)
        ref = Dataset()
        ref.ReferencedBeamNumber = beam.number
        ref.BeamMeterset = beam.mu
        group.ReferencedBeamSequence.append(ref)
    from .dicom import format_decimal_values
    format_decimal_values(ds)
    return ds
