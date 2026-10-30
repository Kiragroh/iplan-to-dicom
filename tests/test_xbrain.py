import struct
import zlib
from pathlib import Path

import pytest

from iplan_dicom.xbrain import FormatError, decode_blzip, read_document, safe_join, numbers, scalar, xml_from_bytes, resolve_inside


def pack(data):
    return b"BLZIP" + struct.pack(">I", len(data)) + zlib.compress(data)


def test_binary_bytes_survive_decompression():
    binary = bytes(range(256)) * 2
    assert decode_blzip(pack(binary)) == binary


@pytest.mark.parametrize("corrupt", [b"BLZIP", b"BLZIP" + struct.pack(">I", 99) + zlib.compress(b"abc"), pack(b"abc")[:-2], pack(b"abc") + b"junk"])
def test_invalid_stream_rejected(corrupt):
    with pytest.raises(FormatError):
        decode_blzip(corrupt)


def test_revision_pointer_selects_named_revision_not_latest(tmp_path):
    pointer = b'<?xml version="1.0"?><xbrain:Revision name="Object" revision="1"/>'
    tmp_path.joinpath("Object.xbrain").write_bytes(pack(bytes(b ^ 255 for b in pointer)))
    xml = b'<xbrain:File xmlns:xbrain="http://brainlab.com/XBrain/1.5" revision="1"><xbrain:Content><xbrain:ContentVersion value="2"><xbrain:Plan/></xbrain:ContentVersion></xbrain:Content></xbrain:File>'
    tmp_path.joinpath("Object.xbrain,1").write_bytes(pack(bytes(b ^ 255 for b in xml)))
    tmp_path.joinpath("Object.xbrain,2").write_bytes(pack(b"do not select this"))
    doc = read_document(tmp_path / "Object.xbrain")
    assert doc.tag == "File"
    assert doc.find("./Content/ContentVersion/Plan") is not None


def test_missing_revision_is_explicit(tmp_path):
    tmp_path.joinpath("Object.xbrain").write_bytes(b'<xbrain:Revision name="Object" revision="2"/>')
    with pytest.raises(FormatError, match="revision"):
        read_document(tmp_path / "Object.xbrain")


@pytest.mark.parametrize("name", ["../outside", "..\\outside", "C:\\outside", "//server/share", "/outside"])
def test_references_cannot_escape_case(tmp_path, name):
    with pytest.raises(FormatError):
        safe_join(tmp_path, name)


def test_ascii_counted_double_vector():
    assert numbers("3 1.5 -2 4", counted=True) == [1.5, -2.0, 4.0]
    with pytest.raises(FormatError):
        numbers("4 1.5 -2 4", counted=True)


def test_native_comma_vector_and_nested_scalar():
    assert numbers("3,1.5,-2,4", counted=True) == [1.5, -2, 4]
    node = xml_from_bytes(b'<File><Dose><_DOUBLE_>\n<_ASCII_>2.5</_ASCII_>\n</_DOUBLE_></Dose></File>')
    assert scalar(node, "Dose") == 2.5


def test_legacy_machine_namespace_fragment():
    node = xml_from_bytes(b'<rtmachine:MachineProfileFileable><rtmachine:LINACParameters/></rtmachine:MachineProfileFileable>')
    assert node.find('LINACParameters') is not None


def test_parent_reference_allowed_only_inside_case(tmp_path):
    child = tmp_path / 'plan'
    child.mkdir()
    assert resolve_inside(tmp_path, child, '..\\images') == tmp_path / 'images'
    with pytest.raises(FormatError):
        resolve_inside(tmp_path, child, '..\\..\\outside')
