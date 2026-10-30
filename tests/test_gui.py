"""Controller-state regressions without opening a desktop window."""
from pathlib import Path
from types import SimpleNamespace
import zipfile

from iplan_dicom.gui import ConverterWindow


class Box(dict):
    def __init__(self, index=0):
        super().__init__()
        self.index, self.value = index, 'old plan'

    def current(self, index=None):
        if index is not None:
            self.index = index
        return self.index

    def set(self, value):
        self.value = value
        if not value:
            self.index = -1


def test_case_scan_clears_old_plan_before_any_async_result():
    window = ConverterWindow.__new__(ConverterWindow)
    window.case_box = Box(1)
    window.plan_box = Box(0)
    window.cases = [Path('synthetic-first'), Path('synthetic-second')]
    window.plans = [{'index': 4, 'supported_plan': True, 'dose_present': True}]
    window.busy = False
    window.append = lambda message: None
    invoked = []
    window.run = lambda operation, kind: invoked.append(kind)
    window.scan()
    assert window.plans == []
    assert window.plan_box['values'] == []
    assert window.plan_box.current() == -1
    assert invoked == ['plans']


def test_gui_retains_zip_until_source_change_or_idle_close(tmp_path):
    from iplan_dicom.source import prepare_source
    opened = []
    for index in range(2):
        path = tmp_path / f'synthetic{index}.zip'
        with zipfile.ZipFile(path, 'w') as archive:
            archive.writestr('case/patient.xbrain', 'synthetic')
        opened.append(prepare_source(path))
    old, new = opened
    old_root, new_root = old.root, new.root
    window = ConverterWindow.__new__(ConverterWindow)
    window.source_session = old
    window.source = SimpleNamespace(set=lambda value: None)
    window.case_box, window.plan_box = Box(), Box()
    window.scan = lambda: None
    window.append = lambda message: None
    window.plan_changed = lambda: None
    destroyed = []
    window.root = SimpleNamespace(destroy=lambda: destroyed.append(True))
    window.apply_source(new, [new.root])
    assert not old_root.exists()
    assert new_root.exists()
    window.busy = True
    window.close()
    assert new_root.exists() and not destroyed
    window.busy = False
    window.close()
    assert not new_root.exists() and destroyed
