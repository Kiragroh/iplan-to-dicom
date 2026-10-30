"""Folder and generic ZIP input, with bounded private temporary extraction.

This supports ZIPs that contain the existing native iPlan case layout. The
other vendor archive layouts have not been supplied or verified; unsupported layouts
fail explicitly instead of inferring a patient or changing relative references.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path, PureWindowsPath
import stat
import tempfile
import unicodedata
import zipfile

from .xbrain import FormatError


@dataclass(frozen=True)
class ZipLimits:
    max_files: int = 100_000
    max_expanded_bytes: int = 20 * 1024**3
    max_member_bytes: int = 2 * 1024**3


class SourceInput:
    """Opened source; GUI retains this object while its extracted files are used."""

    def __init__(self, original: Path, root: Path, temporary=None):
        self.original = original
        self.root = root
        self._temporary = temporary
        self._owned_root = Path(temporary.name).resolve() if temporary is not None else None

    @property
    def is_zip(self):
        return self._temporary is not None

    def close(self):
        if self._temporary is not None:
            # Delete only the exact temporary directory this adapter created.
            if Path(self._temporary.name).resolve() != self._owned_root:
                raise FormatError('Temporary source cleanup target changed')
            self._temporary.cleanup()
            self._temporary = None

    def check_output(self, output: Path):
        if self._owned_root is not None and Path(output).resolve().is_relative_to(self._owned_root):
            raise FormatError('Output must be outside the temporary ZIP extraction')

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.close()


def _members(archive: zipfile.ZipFile, limits: ZipLimits):
    entries = archive.infolist()
    if len(entries) > limits.max_files:
        raise FormatError('ZIP exceeds the entry-count limit')
    total, explicit, nodes, result = 0, set(), {}, []
    reserved = {'CON', 'PRN', 'AUX', 'NUL', 'CONIN$', 'CONOUT$', 'CLOCK$',
                *(f'COM{i}' for i in '123456789¹²³'), *(f'LPT{i}' for i in '123456789¹²³')}
    for entry in entries:
        raw = entry.orig_filename.replace('\\', '/')
        windows = PureWindowsPath(raw)
        if windows.drive or windows.root or raw.startswith('/') or '\x00' in raw:
            raise FormatError('ZIP contains an absolute or invalid member path')
        directory = entry.is_dir() or raw.endswith('/')
        if directory and entry.file_size:
            raise FormatError('ZIP directory entry contains unexpected data')
        parts = tuple(raw[:-1].split('/') if directory else raw.split('/'))
        if (not parts or any(not p or p in ('.', '..') or ':' in p or p.rstrip(' .') != p or
                             any(ord(c) < 32 or c in '<>"|?*' for c in p) or
                             p.split('.')[0].upper() in reserved for p in parts)):
            raise FormatError('ZIP contains an unsafe member path')
        mode = entry.external_attr >> 16
        kind = stat.S_IFMT(mode)
        if (kind not in (0, stat.S_IFREG, stat.S_IFDIR) or (kind == stat.S_IFDIR and not directory)
                or entry.external_attr & 0x400):
            raise FormatError('ZIP links and special files are unsupported')
        if entry.flag_bits & 1:
            raise FormatError('Encrypted ZIP archives are unsupported')
        canonical = tuple(unicodedata.normalize('NFC', p).casefold() for p in parts)
        if canonical in explicit:
            raise FormatError('ZIP has duplicate or case-colliding member paths')
        explicit.add(canonical)
        for length in range(1, len(parts) + 1):
            key, original = canonical[:length], parts[:length]
            is_dir = length < len(parts) or directory
            previous = nodes.get(key)
            if previous is not None and previous != (original, is_dir):
                raise FormatError('ZIP has conflicting file, directory or case-colliding paths')
            nodes[key] = (original, is_dir)
        if entry.file_size < 0 or entry.file_size > limits.max_member_bytes:
            raise FormatError('ZIP member exceeds the expanded-size limit')
        total += entry.file_size
        if total > limits.max_expanded_bytes:
            raise FormatError('ZIP exceeds the total expanded-size limit')
        result.append((entry, parts, directory))
    return result


def _native_root(extracted: Path) -> Path:
    cases = sorted(p.parent for p in extracted.rglob('*') if p.is_file() and p.name.casefold() == 'patient.xbrain')
    if not cases:
        raise FormatError('ZIP has no supported native iPlan case layout. Other vendor archive layouts are not yet verified.')
    if len(cases) == 1:
        return cases[0]
    parents = {p.parent for p in cases}
    if len(parents) != 1:
        raise FormatError('ZIP cases are not one unambiguous immediate-child collection')
    return parents.pop()


def prepare_source(path: Path, *, limits: ZipLimits = ZipLimits()) -> SourceInput:
    """Open a directory unchanged, or safely extract one generic native-layout ZIP."""
    path = Path(path).resolve()
    if path.is_dir():
        return SourceInput(path, path)
    if not path.is_file() or path.suffix.lower() != '.zip':
        raise FormatError('Select a source directory or ZIP archive')
    temporary = tempfile.TemporaryDirectory(prefix='iplan-source-')
    extraction = Path(temporary.name).resolve()
    source = SourceInput(path, extraction, temporary)
    try:
        with zipfile.ZipFile(path) as archive:
            members = _members(archive, limits)
            written = 0
            for entry, parts, directory in members:
                destination = extraction.joinpath(*parts)
                if not destination.resolve().is_relative_to(extraction):
                    raise FormatError('ZIP member escapes the temporary extraction directory')
                if directory:
                    destination.mkdir(parents=True, exist_ok=True)
                    continue
                destination.parent.mkdir(parents=True, exist_ok=True)
                member_size = 0
                with archive.open(entry) as input_stream, destination.open('xb') as output_stream:
                    while block := input_stream.read(1024 * 1024):
                        written += len(block)
                        member_size += len(block)
                        if (written > limits.max_expanded_bytes or member_size > limits.max_member_bytes
                                or member_size > entry.file_size):
                            raise FormatError('ZIP expanded data exceed the declared or configured limits')
                        output_stream.write(block)
                if member_size != entry.file_size:
                    raise FormatError('ZIP member expanded length differs from its declaration')
        source.root = _native_root(extraction)
        return source
    except FormatError:
        source.close()
        raise
    except Exception:
        source.close()
        raise FormatError('ZIP extraction failed; source details omitted') from None
