"""Native Tkinter interface. No shell launchers or browser/native bridge."""
from __future__ import annotations

from pathlib import Path
import queue
import threading
import tkinter as tk
from tkinter import filedialog, ttk

from .service import inspect_case
from .batch import convert_plan, convert_batch, discover_cases
from .source import prepare_source
from .xbrain import FormatError


class ConverterWindow:
    def __init__(self, root):
        self.root = root
        self.events = queue.Queue()
        self.cases = []
        self.plans = []
        self.busy = False
        self.source_session = None
        root.protocol('WM_DELETE_WINDOW', self.close)
        root.title('BrainScan · iPlan → DICOM')
        root.geometry('950x740')
        root.minsize(760, 520)
        self.source = tk.StringVar()
        self.output = tk.StringVar()
        self.include_mr = tk.BooleanVar(value=True)
        self.resume = tk.BooleanVar(value=True)
        self.status = tk.StringVar(value='Quellordner wählen. Die Quelldaten werden ausschließlich gelesen.')
        body = ttk.Frame(root, padding=22)
        body.pack(fill='both', expand=True)
        body.columnconfigure(1, weight=1)
        ttk.Label(body, text='iPlan → DICOM', font=('Segoe UI', 22, 'bold')).grid(row=0, column=0, columnspan=3, sticky='w')
        ttk.Label(body, text='Bilder · Strukturen · gespeicherte Dosis · Bestrahlungsplan').grid(row=1, column=0, columnspan=3, sticky='w', pady=(0, 20))
        ttk.Label(body, text='Quelle').grid(row=2, column=0, sticky='w', padx=(0, 12))
        ttk.Entry(body, textvariable=self.source, state='readonly').grid(row=2, column=1, sticky='ew')
        source_actions = ttk.Frame(body)
        source_actions.grid(row=2, column=2, padx=(10, 0))
        self.source_button = ttk.Button(source_actions, text='Ordner', command=self.choose_source)
        self.source_button.pack(side='left')
        self.zip_button = ttk.Button(source_actions, text='ZIP', command=self.choose_zip)
        self.zip_button.pack(side='left', padx=(5, 0))
        ttk.Label(body, text='Fall').grid(row=3, column=0, sticky='w', pady=12)
        self.case_box = ttk.Combobox(body, state='readonly')
        self.case_box.grid(row=3, column=1, sticky='ew')
        self.case_box.bind('<<ComboboxSelected>>', lambda event: self.scan())
        self.scan_button = ttk.Button(body, text='Pläne lesen', command=self.scan)
        self.scan_button.grid(row=3, column=2, padx=(10, 0))
        ttk.Label(body, text='Plan').grid(row=4, column=0, sticky='w')
        self.plan_box = ttk.Combobox(body, state='readonly')
        self.plan_box.grid(row=4, column=1, columnspan=2, sticky='ew')
        self.plan_box.bind('<<ComboboxSelected>>', lambda event: self.plan_changed())
        ttk.Label(body, text='Exportwurzel').grid(row=5, column=0, sticky='w', pady=12)
        ttk.Entry(body, textvariable=self.output).grid(row=5, column=1, sticky='ew')
        self.output_button = ttk.Button(body, text='Ziel wählen', command=self.choose_output)
        self.output_button.grid(row=5, column=2, padx=(10, 0))
        ttk.Checkbutton(body, text='MR-Serien zusätzlich in ihren ursprünglichen Bezugssystemen erhalten', variable=self.include_mr).grid(row=6, column=0, columnspan=3, sticky='w')
        ttk.Checkbutton(body, text='Bereits vollständig geprüfte Exporte überspringen (Fortsetzen)', variable=self.resume).grid(row=7, column=0, columnspan=3, sticky='w')
        actions = ttk.Frame(body)
        actions.grid(row=8, column=0, columnspan=3, sticky='w', pady=(18, 12))
        self.convert_button = ttk.Button(actions, text='Ausgewählten Plan exportieren', command=self.convert, state='disabled')
        self.convert_button.pack(side='left')
        self.batch_button = ttk.Button(actions, text='Alle Fälle und Pläne exportieren', command=self.convert_all, state='disabled')
        self.batch_button.pack(side='left', padx=14)
        self.bar = ttk.Progressbar(body, mode='indeterminate')
        self.bar.grid(row=9, column=0, columnspan=3, sticky='ew')
        ttk.Label(body, textvariable=self.status, wraplength=860).grid(row=10, column=0, columnspan=3, sticky='w', pady=10)
        self.log = tk.Text(body, height=11, wrap='word', state='disabled', font=('Segoe UI', 10))
        self.log.grid(row=11, column=0, columnspan=3, sticky='nsew')
        body.rowconfigure(11, weight=1)
        ttk.Label(body, text='Entwicklungsversion · Technische Konvertierung ist keine klinische Freigabe.', foreground='#755d2b').grid(row=12, column=0, columnspan=3, sticky='w', pady=(12, 0))
        root.after(100, self.poll)

    def append(self, message):
        self.log.configure(state='normal')
        self.log.insert('end', message + '\n')
        self.log.see('end')
        self.log.configure(state='disabled')
        self.status.set(message)

    def choose_source(self):
        folder = filedialog.askdirectory(title='iPlan-Fallordner oder Ordner mit mehreren Fällen wählen')
        if folder:
            self.open_source(Path(folder))

    def choose_zip(self):
        archive = filedialog.askopenfilename(title='ZIP mit nativer iPlan-Fallstruktur wählen',
                                            filetypes=[('ZIP-Archive', '*.zip')])
        if archive:
            self.open_source(Path(archive))

    def open_source(self, path):
        if self.busy:
            return
        self.append('Quelle prüfen; ZIP-Inhalte werden vorübergehend lokal entpackt …')
        def prepare():
            opened = prepare_source(path)
            try:
                return opened, discover_cases(opened.root)
            except Exception:
                opened.close()
                raise
        self.run(prepare, 'source_ready')

    def apply_source(self, opened, cases):
        try:
            if self.source_session is not None:
                self.source_session.close()
        except Exception:
            opened.close()
            self.cases, self.plans = [], []
            self.append('Temporäre Quelle konnte nicht bereinigt werden. Bitte erneut versuchen.')
            return
        self.source_session = opened
        self.source.set(str(opened.original))
        self.cases = cases
        self.case_box['values'] = [p.name for p in self.cases]
        self.plan_box.set('')
        self.plans = []
        if self.cases:
            self.case_box.current(0)
            self.scan()
        else:
            self.case_box.set('')
            self.append('Keine iPlan-Fälle mit patient.xbrain gefunden.')
            self.plan_changed()

    def close(self):
        if self.busy:
            self.append('Bitte die laufende Verarbeitung abwarten; die temporäre Quelle wird noch benötigt.')
            return
        try:
            if self.source_session is not None:
                self.source_session.close()
        except Exception:
            self.append('Temporäre Quelle konnte nicht bereinigt werden. Bitte erneut schließen.')
            return
        self.root.destroy()

    def choose_output(self):
        parent = filedialog.askdirectory(title='Übergeordneten Ausgabeordner wählen')
        if parent:
            self.output.set(parent)

    def set_busy(self, value):
        self.busy = value
        for button in (self.source_button, self.zip_button, self.scan_button, self.output_button):
            button.configure(state='disabled' if value else 'normal')
        self.case_box.configure(state='disabled' if value else 'readonly')
        self.plan_box.configure(state='disabled' if value else 'readonly')
        self.convert_button.configure(state='disabled')
        self.batch_button.configure(state='disabled')
        self.bar.start(12) if value else self.bar.stop()
        if not value:
            self.plan_changed()

    def run(self, operation, kind):
        if self.busy:
            return
        self.set_busy(True)
        def worker():
            try:
                self.events.put((kind, operation()))
            except FormatError as exc:
                self.events.put(('error', str(exc)))
            except Exception as exc:
                self.events.put(('error', f'Verarbeitung gestoppt ({type(exc).__name__}). Keine vollständige Ausgabe bestätigt.'))
        threading.Thread(target=worker, daemon=True).start()

    def scan(self):
        index = self.case_box.current()
        if index < 0 or self.busy:
            return
        case = self.cases[index]
        # A failed scan must not leave the previous case's plan selectable.
        self.plans = []
        self.plan_box['values'] = []
        self.plan_box.set('')
        self.append('Planstände und gespeicherte Daten prüfen …')
        self.run(lambda: inspect_case(case, include_names=True), 'plans')

    def plan_changed(self):
        index = self.plan_box.current()
        valid = index >= 0 and index < len(self.plans) and self.plans[index]['supported_plan'] and self.plans[index].get('dose_present', False)
        self.convert_button.configure(state='normal' if valid and not self.busy else 'disabled')
        self.batch_button.configure(state='normal' if self.cases and not self.busy else 'disabled')

    def convert(self):
        case_index, plan_index = self.case_box.current(), self.plan_box.current()
        if case_index < 0 or plan_index < 0 or not self.output.get().strip():
            self.append('Bitte Fall, Plan und eine Exportwurzel auswählen.')
            return
        case, output = self.cases[case_index], Path(self.output.get())
        try:
            if self.source_session is not None:
                self.source_session.check_output(output)
        except FormatError as exc:
            self.append(str(exc))
            return
        native_index = self.plans[plan_index]['index']
        include_mr = self.include_mr.get()
        resume = self.resume.get()
        self.run(lambda: convert_plan(case, output, plan_index=native_index, include_mr=include_mr, resume=resume,
                                     progress=lambda message: self.events.put(('progress', message))), 'done')

    def convert_all(self):
        if not self.source.get().strip() or not self.output.get().strip():
            self.append('Bitte Quelle und Exportwurzel auswählen.')
            return
        source = self.source_session.root if self.source_session is not None else Path(self.source.get())
        output = Path(self.output.get())
        try:
            if self.source_session is not None:
                self.source_session.check_output(output)
        except FormatError as exc:
            self.append(str(exc))
            return
        include_mr, resume = self.include_mr.get(), self.resume.get()
        self.append('Batch-Export: Alle aktuellen Planvarianten werden einzeln geprüft. Fehler stoppen nur den betroffenen Plan.')
        self.run(lambda: convert_batch(source, output, include_mr=include_mr, resume=resume,
                                      progress=lambda message: self.events.put(('progress', message))), 'batch_done')

    def poll(self):
        try:
            while True:
                kind, value = self.events.get_nowait()
                if kind == 'progress':
                    self.append(value)
                    continue
                self.set_busy(False)
                if kind == 'source_ready':
                    self.apply_source(*value)
                elif kind == 'plans':
                    self.plans = value['plans']
                    labels = []
                    for plan in self.plans:
                        if plan['supported_plan']:
                            labels.append(f"{plan['index']}: {plan.get('name', 'Plan')} · {plan['beams']} Felder · {plan['fractions']} Fraktionen")
                        else:
                            labels.append(f"{plan['index']}: Nicht unterstützt – {plan['reason']}")
                    self.plan_box['values'] = labels
                    self.plan_box.set('Planstand ausdrücklich auswählen')
                    self.append(f'{len(labels)} Planstände gefunden. Bitte den gewünschten Stand auswählen.')
                elif kind == 'done':
                    if value['status'] == 'resumed_verified':
                        self.append('Vorhandener Export vollständig anhand Quelldaten und Datei-Prüfsummen bestätigt.')
                    else:
                        self.append(f"Geschrieben und zurückgelesen: {value['ct_slices']} CT-Schichten, {value['rois']} Strukturen, {value['beams']} Felder, RTDOSE.")
                    for warning in value.get('warnings', []):
                        self.append(warning)
                elif kind == 'batch_done':
                    self.append(f"Batch beendet: {value['converted']} konvertiert, {value['resumed_verified']} geprüft übersprungen, {value['failed']} fehlgeschlagen. Manifest: {value['manifest']}")
                    for entry in value['entries']:
                        if entry['status'] == 'failed':
                            self.append(f"Fall {entry['case_index']}, Plan {entry.get('plan_index', '?')}: {entry.get('reason', 'Verarbeitung fehlgeschlagen')}")
                else:
                    self.append('Gestoppt: ' + value)
                self.plan_changed()
        except queue.Empty:
            pass
        self.root.after(100, self.poll)


def main():
    root = tk.Tk()
    ConverterWindow(root)
    root.mainloop()
