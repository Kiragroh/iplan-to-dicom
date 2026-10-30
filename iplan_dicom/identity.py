"""Identity reconciliation restricted to the selected plan's resolved image graph.

Native Patient/Study metadata do not contain typed ownership edges. Ownership
here means a selected RTPlan SliceSet or a SliceSet reachable through its
verified Fusion graph, followed by the exact ImageSeries and original SOPs.
This is not permission to merge arbitrary same-name images from a case folder.
"""
from __future__ import annotations

from copy import deepcopy
from dataclasses import replace
from pathlib import Path
import unicodedata
import xml.etree.ElementTree as ET

from .images import ImageVolume, registration_matrix, validate_source_headers
from .plan import objects, referenced_object, resolve_directory
from .xbrain import FormatError


def _demographics(header) -> tuple[str, str, str]:
    # Keep PN component/group separators; do not use fuzzy name matching.
    name = unicodedata.normalize('NFC', str(header.get('PatientName', ''))).casefold()
    name = '='.join('^'.join(' '.join(part.split()) for part in group.split('^'))
                    for group in name.split('='))
    dob = str(header.get('PatientBirthDate', '')).strip()
    sex = str(header.get('PatientSex', '')).strip().upper()
    if not name.replace('^', '').replace('=', '').strip() or not dob or not sex:
        raise FormatError('Source image identity review required')
    return name, dob, sex


class NativeImageGraph:
    """Resolve source ownership and registration before considering identity."""

    def __init__(self, plan, case_dir: Path):
        self.series = {}
        self.grids = {}
        self.reconciled_series = set()
        case_dir = Path(case_dir).resolve()
        refs = [r for r in plan.source.findall('Reference') if r.get('contentType') == 'SliceSet']
        if len(refs) != 1:
            raise FormatError('Plan image ownership is ambiguous')
        self.primary_uid = refs[0].get('contentUid', '')
        if not self.primary_uid:
            raise FormatError('Plan image ownership is missing')
        all_refs = list(refs)
        edges = {}
        for fusion in objects(plan.source_directory / 'fusion.xbrain', 'Fusion'):
            linked = fusion.findall('Reference')
            if len(linked) != 2 or any(r.get('contentType') != 'SliceSet' or not r.get('contentUid') for r in linked):
                raise FormatError('Native image ownership references are invalid')
            a, b = (r.get('contentUid') for r in linked)
            edges.setdefault(a, set()).add(b)
            edges.setdefault(b, set()).add(a)
            all_refs.extend(linked)
        reachable, pending = set(), [self.primary_uid]
        while pending:
            uid = pending.pop()
            if uid in reachable:
                continue
            reachable.add(uid)
            pending.extend(edges.get(uid, ()))
        for uid in reachable:
            if uid != self.primary_uid:
                registration_matrix(uid, self.primary_uid, plan.source_directory)
        for reference in all_refs:
            uid = reference.get('contentUid')
            if uid not in reachable:
                continue
            owner = ET.Element('Owner')
            owner.append(deepcopy(reference))
            grid = referenced_object(owner, 'SliceSet', plan.source_directory, case_dir)
            directory = resolve_directory(case_dir, plan.source_directory, reference.get('directory', '.'))
            existing = self.grids.get(uid)
            if existing is not None:
                if existing[1] != directory or ET.tostring(existing[0]) != ET.tostring(grid):
                    raise FormatError('Native image grid ownership is ambiguous')
                continue
            slices = grid.find('Slices')
            if slices is None:
                raise FormatError('Native image grid has no source series')
            series = referenced_object(slices, 'ImageSeries', directory, case_dir)
            series_ref = [r for r in slices.findall('Reference') if r.get('contentType') == 'ImageSeries'][0]
            series_dir = resolve_directory(case_dir, directory, series_ref.get('directory', '.'))
            series_uid = series.get('uId', '')
            if not series_uid:
                raise FormatError('Native image series ownership is missing')
            previous = self.series.get(series_uid)
            if previous is not None and (previous[1] != series_dir or ET.tostring(previous[0]) != ET.tostring(series)):
                raise FormatError('Native image series ownership is ambiguous')
            self.series[series_uid] = (series, series_dir)
            self.grids[uid] = (grid, directory, series_uid)

    def _check_image(self, image: ImageVolume) -> None:
        validate_source_headers(image.headers, image.modality)
        if image.source_uid in self.grids:
            grid, _, series_uid = self.grids[image.source_uid]
            series = self.series[series_uid][0]
            nodes = sorted(series.findall('Image'), key=lambda n: int(n.get('number', '-1')))
            try:
                selected = [nodes[int(s.attrib['imageIndex'])] for s in grid.findall('Slices/Slice')]
            except (KeyError, ValueError, IndexError):
                raise FormatError('Native image ownership indices are invalid') from None
        elif image.source_uid in self.series:
            selected = sorted(self.series[image.source_uid][0].findall('Image'), key=lambda n: int(n.get('number', '-1')))
        else:
            raise FormatError('Image is not owned by the selected native plan graph')
        expected = [n.get('dicomInstanceUID', '') for n in selected]
        actual = [str(h.get('SOPInstanceUID', '')) for h in image.headers]
        if not expected or not all(expected) or expected != actual:
            raise FormatError('Native image ownership differs from original SOP references')

    def reconcile(self, image: ImageVolume, canonical: ImageVolume) -> ImageVolume:
        """Copy canonical patient identity only after exact graph and demographic checks."""
        if canonical.source_uid != self.primary_uid:
            raise FormatError('Canonical image is not the selected planning grid')
        self._check_image(canonical)
        self._check_image(image)
        canonical_id = str(canonical.headers[0].PatientID)
        identities_differ = any(str(h.PatientID) != canonical_id for h in image.headers)
        # All original planes must agree, including the canonical CT planes.
        expected = _demographics(canonical.headers[0])
        if any(_demographics(h) != expected for h in (*canonical.headers, *image.headers)):
            raise FormatError('Source image identity review required')
        if not identities_differ:
            return image
        copied = []
        for original in image.headers:
            header = deepcopy(original)
            for key in ('PatientID', 'PatientName', 'PatientBirthDate', 'PatientSex'):
                setattr(header, key, deepcopy(getattr(canonical.headers[0], key)))
            copied.append(header)
        self.reconciled_series.add(image.source_uid)
        return replace(image, headers=tuple(copied))
