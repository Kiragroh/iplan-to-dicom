"""BLZIP decoding and explicit XBrain revision resolution.

Errors deliberately omit source paths and XML values: these may identify patients.
"""

from __future__ import annotations

import math
from pathlib import Path, PureWindowsPath
import re
import xml.etree.ElementTree as ET
import zlib

INVERT = bytes.maketrans(bytes(range(256)), bytes(v ^ 255 for v in range(256)))
MAX_BYTES = 1024 * 1024 * 1024


class FormatError(ValueError):
    """Unsupported or inconsistent source data; never substitute a guessed value."""


def decode_blzip(data: bytes, *, max_bytes: int = MAX_BYTES) -> bytes:
    """Return the uncompressed bytes, without changing their encoding."""
    if not data.startswith(b"BLZIP"):
        if len(data) > max_bytes:
            raise FormatError("Source exceeds size limit")
        return data
    if len(data) < 11:
        raise FormatError("Truncated BLZIP header")
    length = int.from_bytes(data[5:9], "big")
    if length > max_bytes:
        raise FormatError("BLZIP length exceeds size limit")
    stream = zlib.decompressobj()
    try:
        decoded = stream.decompress(data[9:], length + 1)
    except zlib.error:
        raise FormatError("Invalid BLZIP compressed stream") from None
    if len(decoded) != length or not stream.eof or stream.unused_data or stream.unconsumed_tail:
        raise FormatError("BLZIP length or stream integrity mismatch")
    return decoded


def read_payload(path: Path, *, invert: bool = True) -> bytes:
    """Read an XBrain payload. Binary readers explicitly choose inversion."""
    try:
        with Path(path).open("rb") as source:
            data = source.read(MAX_BYTES + 10)
    except OSError:
        raise FormatError("Source file is unavailable") from None
    decoded = decode_blzip(data)
    return decoded.translate(INVERT) if invert else decoded


def xml_from_bytes(data: bytes) -> ET.Element:
    """Parse normal or inverted XML, including the standalone Revision fragment."""
    data = decode_blzip(data)
    if not data.lstrip().startswith(b"<"):
        data = data.translate(INVERT)
    if not data.lstrip().startswith(b"<"):
        raise FormatError("Payload is binary, not an XML document")
    if b"<!DOCTYPE" in data.upper() or b"<!ENTITY" in data.upper():
        raise FormatError("XML declarations with external or expanded entities are unsupported")
    # Native revision pointers omit the namespace declaration, unlike File documents.
    if b"<xbrain:Revision" in data and b"xmlns:xbrain=" not in data:
        data = data.replace(b"<xbrain:Revision", b'<xbrain:Revision xmlns:xbrain="http://brainlab.com/XBrain/1.5"', 1)
    for prefix in (b"rtmachine", b"settings", b"xvmc"):
        if b"<" + prefix + b":" in data and b"xmlns:" + prefix + b"=" not in data:
            # These legacy fragments use a prefix without its declaration.
            # A parsing-only namespace preserves local tag names; it is not
            # presented as an authoritative vendor schema URI.
            root_match = re.search(rb"<[A-Za-z_][A-Za-z_0-9:.-]*", data)
            if root_match:
                end = root_match.end()
                data = data[:end] + b' xmlns:' + prefix + b'="urn:iplan-legacy:' + prefix + b'"' + data[end:]
    try:
        root = ET.fromstring(data)
    except ET.ParseError:
        raise FormatError("Invalid XBrain XML document") from None
    for node in root.iter():
        node.tag = node.tag.rsplit("}", 1)[-1]
    return root


def read_document(path: Path, *, revision: int | None = None) -> ET.Element:
    """Resolve a native Revision pointer; never infer a revision by file recency."""
    path = Path(path)
    if revision is not None:
        if revision < 0:
            raise FormatError("Invalid revision")
        path = path.with_name(re.sub(r",\d+$", "", path.name) + f",{revision}")
    try:
        root = xml_from_bytes(path.read_bytes())
    except OSError:
        raise FormatError("Requested source revision is unavailable") from None
    if root.tag != "Revision":
        return root
    rev = root.get("revision", "")
    if not rev.isdigit() or "," in path.suffix:
        raise FormatError("Invalid or recursive revision pointer")
    # The pointer's name is a historical absolute database path. Only the local
    # sibling is authoritative after a case has been copied to another location.
    sibling = path.with_name(path.name + "," + rev)
    try:
        resolved = xml_from_bytes(sibling.read_bytes())
    except OSError:
        raise FormatError("Referenced revision is unavailable") from None
    if resolved.tag == "Revision":
        raise FormatError("Recursive revision pointer")
    if resolved.get("revision") not in (None, rev):
        raise FormatError("Resolved revision does not match pointer")
    return resolved


def safe_join(root: Path, relative: str) -> Path:
    """Resolve case-relative references, rejecting drive/UNC/parent escapes."""
    windows = PureWindowsPath(relative)
    if windows.drive or windows.root or ".." in windows.parts:
        raise FormatError("Reference escapes the selected case")
    candidate = Path(root).joinpath(*windows.parts).resolve()
    try:
        candidate.relative_to(Path(root).resolve())
    except ValueError:
        raise FormatError("Reference escapes the selected case") from None
    return candidate


def resolve_inside(case_root: Path, base: Path, relative: str) -> Path:
    """Resolve native parent-relative links but never leave the selected case."""
    parts = PureWindowsPath(relative)
    if parts.drive or parts.root:
        raise FormatError("Absolute source reference is unsupported")
    candidate = Path(base).joinpath(*parts.parts).resolve()
    try:
        candidate.relative_to(Path(case_root).resolve())
    except ValueError:
        raise FormatError("Reference escapes the selected case") from None
    return candidate


def numbers(text: str | None, *, counted: bool = False) -> list[float]:
    try:
        values = [float(x) for x in re.split(r"[,\s]+", (text or "").strip()) if x]
    except ValueError:
        raise FormatError("Invalid numeric field") from None
    if not values or not all(math.isfinite(v) for v in values):
        raise FormatError("Missing or non-finite numeric field")
    if counted:
        if values[0] != len(values) - 1:
            raise FormatError("Numeric vector count mismatch")
        values = values[1:]
    return values


def scalar(node: ET.Element, path: str) -> float:
    child = node.find(path)
    if child is None:
        raise FormatError("Required numeric field is missing")
    value = child.find("_DOUBLE_")
    if value is None:
        value = child.find("_INTEGER_")
    if value is None:
        value = child
    ascii_value = value.find("_ASCII_")
    vals = numbers(ascii_value.text if ascii_value is not None else value.text)
    if len(vals) != 1:
        raise FormatError("Scalar field contains multiple values")
    return vals[0]
