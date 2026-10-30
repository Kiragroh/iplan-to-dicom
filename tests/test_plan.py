from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np
import pytest
from pydicom.dataset import Dataset

from iplan_dicom.plan import (Beam, Machine, Plan, content, populate_rtplan_dataset,
                              read_points, resolve_directory, blocking_devices,
                              local_component, treatment_elements)
from iplan_dicom.xbrain import FormatError


def value(parent, tag, values):
    node = ET.SubElement(parent, tag)
    ET.SubElement(ET.SubElement(node, "_DOUBLE_"), "_ASCII_").text = str(values)


def synthetic_arc(direction="counterclockwise"):
    arc = ET.Element("ConfArcElement", arcType="dynamic", desiredGantryRotationDirection=direction)
    for index, (angle, weight) in enumerate([(350, 0), (360, 0.2), (370, 1)]):
        point = ET.SubElement(arc, "ConformalField")
        for tag, number in [("GantryAngle", angle), ("TableAngle", -10),
                            ("CollimatorAngle", 0), ("CumMeterSetWeight", weight)]:
            value(point, tag, number)
        jaws = ET.SubElement(point, "JawsState")
        for tag, number in zip(["X1", "X2", "Y1", "Y2"], [-20, 20, -10, 10]):
            value(jaws, tag, number)
        leaves = ET.SubElement(point, "MultileafState")
        value(leaves, "X1", f"2,{-2-index},-1")
        value(leaves, "X2", "2,1,2")
    return arc


def test_counterclockwise_reverses_geometry_and_complements_weights():
    points, direction = read_points(synthetic_arc(), 2)
    assert direction == "CC"
    assert [p.gantry for p in points] == [10, 0, 350]
    assert [p.weight for p in points] == [0, 0.8, 1]
    assert [p.leaves[0] for p in points] == [-4, -3, -2]
    assert all(p.table == 350 for p in points)


def test_clockwise_preserves_geometry_and_nonuniform_weights():
    points, direction = read_points(synthetic_arc("clockwise"), 2)
    assert direction == "CW"
    assert [p.gantry for p in points] == [350, 0, 10]
    assert [p.weight for p in points] == [0, 0.2, 1]


def test_static_native_minus_one_weight_becomes_two_dicom_points():
    arc = synthetic_arc()
    arc.tag = "ConfBeamElement"
    for point in list(arc)[1:]:
        arc.remove(point)
    arc.find("ConformalField/CumMeterSetWeight/_DOUBLE_/_ASCII_").text = "-1"
    points, direction = read_points(arc, 2)
    assert direction == "NONE"
    assert [p.weight for p in points] == [0, 1]
    assert np.array_equal(points[0].leaves, points[1].leaves)


@pytest.mark.parametrize("weight", [None, "0", "1", "0.5", "nan"])
def test_static_missing_or_arbitrary_meterset_is_rejected(weight):
    element = synthetic_arc()
    element.tag = "ConfBeamElement"
    for point in list(element)[1:]:
        element.remove(point)
    point = element[0]
    if weight is None:
        point.remove(point.find("CumMeterSetWeight"))
    else:
        point.find("CumMeterSetWeight/_DOUBLE_/_ASCII_").text = weight
    with pytest.raises(FormatError):
        read_points(element, 2)


def write_collection(path, elements, count=None):
    root = ET.Element("File")
    container = ET.SubElement(root, "Content", count=str(len(elements) if count is None else count))
    if elements:
        multiple = ET.SubElement(ET.SubElement(container, "ContentVersion", value="4"), "Multiple")
        multiple.extend(elements)
    path.write_bytes(ET.tostring(root))


@pytest.mark.parametrize("directory", ["../other", "../../escape"])
def test_component_uid_cannot_override_reference_directory(tmp_path, directory):
    plan = tmp_path / "plan"
    plan.mkdir()
    component = ET.Element("LINACIsocComponent", uId="synthetic-component")
    element = ET.Element("ConfArcElement")
    ET.SubElement(element, "Reference", contentType="LINACIsocComponent",
                  contentUid="synthetic-component", directory=directory)
    with pytest.raises(FormatError):
        local_component(element, [component], plan, tmp_path)


@pytest.mark.parametrize("count", [0, 2])
def test_component_reference_requires_exactly_one_occurrence(tmp_path, count):
    component = ET.Element("LINACIsocComponent", uId="synthetic-component")
    element = ET.Element("ConfArcElement")
    for _ in range(count):
        ET.SubElement(element, "Reference", contentType="LINACIsocComponent",
                      contentUid="synthetic-component", directory=".")
    with pytest.raises(FormatError, match="exactly one"):
        local_component(element, [component], tmp_path, tmp_path)


def test_local_component_reference_is_resolved_and_checked(tmp_path):
    component = ET.Element("LINACIsocComponent", uId="synthetic-component")
    write_collection(tmp_path / "linacisoccomponent.xbrain", [component])
    element = ET.Element("ConfArcElement")
    ET.SubElement(element, "Reference", contentType="LINACIsocComponent",
                  contentUid="synthetic-component", directory=".")
    assert local_component(element, [component], tmp_path, tmp_path) is component
    component.set("changed", "true")
    with pytest.raises(FormatError, match="differs"):
        local_component(element, [component], tmp_path, tmp_path)


@pytest.mark.parametrize("filename,tag", [("futureelement.xbrain", "FutureElement"),
                                          ("confarcelement.xbrain", "UnexpectedElement")])
def test_inventory_rejects_every_unknown_populated_treatment_collection(tmp_path, filename, tag):
    write_collection(tmp_path / filename, [ET.Element(tag)])
    with pytest.raises(FormatError, match="Unsupported treatment"):
        treatment_elements(tmp_path)


def test_inventory_ignores_empty_unknown_collection_but_not_contradictory_count(tmp_path):
    path = tmp_path / "futureelement.xbrain"
    write_collection(path, [])
    assert treatment_elements(tmp_path) == []
    write_collection(path, [ET.Element("FutureElement")], count=0)
    with pytest.raises(FormatError):
        treatment_elements(tmp_path)


def device_fixture():
    profile = ET.Element("MachineProfileFileable")
    for identifier in ("jaws", "mlc"):
        ET.SubElement(ET.SubElement(profile, "BlockingDeviceId"), "_ASCII_").text = identifier
    jaws = ET.Element("JawsParameters", availableX="true", availableY="true",
                      symmetricX="false", symmetricY="false")
    value(jaws, "Orientation", 0)
    mlc = ET.Element("MultileafParameters")
    return profile, {("JawsParameters", "jaws"): jaws, ("MultileafParameters", "mlc"): mlc}


def test_all_blocking_device_references_must_resolve():
    profile, definitions = device_fixture()
    ET.SubElement(ET.SubElement(profile, "BlockingDeviceId"), "_ASCII_").text = "missing"
    with pytest.raises(FormatError, match="unresolved"):
        blocking_devices(profile, definitions)


def test_additional_resolved_accessory_is_rejected():
    profile, definitions = device_fixture()
    ET.SubElement(ET.SubElement(profile, "BlockingDeviceId"), "_ASCII_").text = "wedge"
    definitions[("WedgeParameters", "wedge")] = ET.Element("WedgeParameters")
    with pytest.raises(FormatError, match="accessory"):
        blocking_devices(profile, definitions)


@pytest.mark.parametrize("attribute,value_text", [("availableX", "false"), ("availableY", "false"),
                                                   ("symmetricX", "true"), ("symmetricY", "true")])
def test_unverified_jaw_configuration_is_rejected(attribute, value_text):
    profile, definitions = device_fixture()
    definitions[("JawsParameters", "jaws")].set(attribute, value_text)
    with pytest.raises(FormatError, match="Jaw availability"):
        blocking_devices(profile, definitions)


def test_supported_jaws_and_mlc_resolve_without_implicit_extra_devices():
    profile, definitions = device_fixture()
    jaws, mlc = blocking_devices(profile, definitions)
    assert jaws is definitions[("JawsParameters", "jaws")]
    assert mlc is definitions[("MultileafParameters", "mlc")]


def test_mlc_count_mismatch_is_not_silently_dropped():
    with pytest.raises(FormatError, match="leaf count"):
        read_points(synthetic_arc(), 3)


def test_moving_couch_does_not_receive_a_false_none_rotation_direction():
    arc = synthetic_arc()
    arc.findall("ConformalField")[1].find("TableAngle/_DOUBLE_/_ASCII_").text = "20"
    with pytest.raises(FormatError, match="Moving couch"):
        read_points(arc, 2)


def test_fixed_aperture_conformal_arc_preserves_gantry_motion():
    arc = synthetic_arc()
    arc.set("arcType", "conformal")
    for point in arc:
        point.find("MultileafState/X1/_DOUBLE_/_ASCII_").text = "2,-2,-1"
    points, direction = read_points(arc, 2)
    assert direction == "CC"
    assert [p.gantry for p in points] == [10, 0, 350]
    assert all(np.array_equal(p.leaves, points[0].leaves) for p in points)


def test_alternative_content_versions_are_not_concatenated():
    root = ET.fromstring('<File><Content><ContentVersion value="4"><Multiple><ConfArcElement/></Multiple></ContentVersion>'
                         '<ContentVersion value="3"><Multiple><ConfArcElement/></Multiple></ContentVersion></Content></File>')
    assert len(list(content(root).iter("ConfArcElement"))) == 1
    assert content(root).get("value") == "4"


def test_parent_reference_inside_case_is_permitted_but_escape_is_rejected(tmp_path):
    plan = tmp_path / "study" / "plan"
    assert resolve_directory(tmp_path, plan, "..\\shared") == tmp_path / "study" / "shared"
    with pytest.raises(FormatError, match="escapes"):
        resolve_directory(tmp_path, plan, "..\\..\\..\\other")


def test_rtplan_preserves_mlc_and_mu_and_transforms_isocenter():
    arc = synthetic_arc()
    points, direction = read_points(arc, 2)
    machine = Machine("SYNTHETIC", "TEST", "TEST", 6.0, 1000.0, np.array([-5, 0, 5]))
    beam = Beam(1, 23.5, 0.2, 300.0, direction, machine, np.array([3, 5, 7, 1]),
                points, arc, Path("synthetic.xbrain"))
    model = Plan([beam], 5, 2.0, ET.Element("RTPlan"), Path("synthetic"))
    ds = populate_rtplan_dataset(Dataset(), model, np.diag([-1, -1, 1, 1]))
    emitted = ds.BeamSequence[0]
    assert emitted.ControlPointSequence[0].IsocenterPosition == [-3, -5, 7]
    assert ds.FractionGroupSequence[0].ReferencedBeamSequence[0].BeamMeterset == 23.5
    assert emitted.BeamLimitingDeviceSequence[-1].LeafPositionBoundaries == [-5, 0, 5]
    assert emitted.ControlPointSequence[0].BeamLimitingDevicePositionSequence[-1].LeafJawPositions == [-4, -1, 1, 2]
    assert emitted.ControlPointSequence[-1].GantryRotationDirection == "NONE"
    assert "BeamDose" not in ds.FractionGroupSequence[0].ReferencedBeamSequence[0]
