#!/usr/bin/env python
"""
LLAMAS Exposure Time Calculator -- PyQt6 GUI.

A single-window front end to observe.observe_spectrum(), for users who prefer a GUI to the Jupyter
notebook. It wraps the SAME calculation engine (spectrograph + observe), so results match the notebook.

Run from anywhere:
    python etc_gui.py

Controls: input spectrum (defaults to the bundled SN1a_R20mag.fits), exposure time, airmass, seeing,
source type (point / extended surface brightness), aperture, PSF, throughput mode, and which channels.
The embedded plot shows counts and SNR vs wavelength; "Save results" writes a PNG + a CSV.
"""
import os
import sys
import io
import csv
import contextlib

import numpy as np
from astropy.table import Table

# run relative to this file so the bare `import spectrograph/observe`, the *.def files, and COATINGS/ resolve
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
os.chdir(SCRIPT_DIR)
sys.path.insert(0, SCRIPT_DIR)
os.environ.setdefault('COATINGS_PATH', './COATINGS/')

from PyQt6 import QtWidgets, QtCore                                      # noqa: E402
from matplotlib.figure import Figure                                    # noqa: E402
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg, NavigationToolbar2QT   # noqa: E402

import spectrograph as spec                                             # noqa: E402
import observe                                                          # noqa: E402

CHANNELS = [('blue', 'LLAMAS_BLUE', 'llamas_blue.def', 'tab:blue'),
            ('green', 'LLAMAS_GREEN', 'llamas_green.def', 'tab:green'),
            ('red', 'LLAMAS_RED', 'llamas_red.def', 'tab:red')]
DEFAULT_SPECTRUM = os.path.join(SCRIPT_DIR, 'SN1a_R20mag.fits')


class ETCWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle('LLAMAS Exposure Time Calculator')
        self._models = {}                 # cache: (channel, mode) -> Spectrograph
        self._results = {}                # channel -> (waves, counts, noise)
        self.wave_nm = None
        self.flux = None

        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        layout = QtWidgets.QHBoxLayout(central)
        layout.addLayout(self._build_controls(), 0)
        layout.addLayout(self._build_plot(), 1)

        self._load_spectrum(DEFAULT_SPECTRUM)
        self._on_source_changed()
        self.compute()

    # ---------------- UI ----------------
    def _build_controls(self):
        form = QtWidgets.QFormLayout()

        # input spectrum path + browse
        self.path_edit = QtWidgets.QLineEdit(DEFAULT_SPECTRUM)
        browse = QtWidgets.QPushButton('Browse...')
        browse.clicked.connect(self._browse)
        row = QtWidgets.QHBoxLayout()
        row.addWidget(self.path_edit); row.addWidget(browse)
        w = QtWidgets.QWidget(); w.setLayout(row)
        form.addRow('Input spectrum:', w)

        self.texp = QtWidgets.QDoubleSpinBox()
        self.texp.setRange(1, 1e6); self.texp.setValue(1200); self.texp.setSuffix(' s')
        form.addRow('Exposure time:', self.texp)

        self.airmass = QtWidgets.QDoubleSpinBox()
        self.airmass.setRange(1.0, 4.0); self.airmass.setSingleStep(0.05); self.airmass.setValue(1.2)
        form.addRow('Airmass (sec z):', self.airmass)

        self.source = QtWidgets.QComboBox(); self.source.addItems(['point', 'extended'])
        self.source.currentTextChanged.connect(self._on_source_changed)
        form.addRow('Source type:', self.source)

        self.seeing = QtWidgets.QDoubleSpinBox()
        self.seeing.setRange(0.0, 3.0); self.seeing.setSingleStep(0.1); self.seeing.setValue(0.8)
        self.seeing.setSuffix(' "')
        form.addRow('Seeing FWHM:', self.seeing)

        self.aperture = QtWidgets.QComboBox()
        self.aperture.addItems(['optimal', 'single', '2', '3', '5', '7', '9'])
        form.addRow('Aperture (fibres):', self.aperture)

        self.psf = QtWidgets.QComboBox(); self.psf.addItems(['moffat', 'gaussian'])
        form.addRow('Seeing PSF:', self.psf)

        self.nbin = QtWidgets.QSpinBox(); self.nbin.setRange(1, 37); self.nbin.setValue(1)
        form.addRow('Bin fibres (extended):', self.nbin)

        self.mode = QtWidgets.QComboBox(); self.mode.addItems(['measured', 'theoretical'])
        form.addRow('Throughput model:', self.mode)

        self.ch_boxes = {}
        chrow = QtWidgets.QHBoxLayout()
        for name, _, _, _ in CHANNELS:
            cb = QtWidgets.QCheckBox(name); cb.setChecked(True)
            self.ch_boxes[name] = cb; chrow.addWidget(cb)
        cw = QtWidgets.QWidget(); cw.setLayout(chrow)
        form.addRow('Channels:', cw)

        compute = QtWidgets.QPushButton('Compute')
        compute.clicked.connect(self.compute)
        save = QtWidgets.QPushButton('Save results...')
        save.clicked.connect(self._save)
        form.addRow(compute); form.addRow(save)

        self.status = QtWidgets.QPlainTextEdit(); self.status.setReadOnly(True)
        self.status.setMaximumWidth(360); self.status.setFixedHeight(150)
        form.addRow('Log:', self.status)
        return form

    def _build_plot(self):
        col = QtWidgets.QVBoxLayout()
        self.fig = Figure(figsize=(7, 6))
        self.canvas = FigureCanvasQTAgg(self.fig)
        self.ax_counts = self.fig.add_subplot(2, 1, 1)
        self.ax_snr = self.fig.add_subplot(2, 1, 2, sharex=self.ax_counts)
        col.addWidget(NavigationToolbar2QT(self.canvas, self))
        col.addWidget(self.canvas)
        return col

    # ------------- behaviour -------------
    def _on_source_changed(self):
        point = self.source.currentText() == 'point'
        for w in (self.seeing, self.aperture, self.psf):
            w.setEnabled(point)
        self.nbin.setEnabled(not point)

    def _browse(self):
        path, _ = QtWidgets.QFileDialog.getOpenFileName(
            self, 'Select input spectrum', SCRIPT_DIR, 'Spectra (*.fits *.dat *.txt *.csv);;All files (*)')
        if path:
            self.path_edit.setText(path)
            self._load_spectrum(path)

    def _load_spectrum(self, path):
        try:
            t = Table.read(path)
            cols = t.colnames
            wcol = next((c for c in cols if 'wave' in c.lower()), cols[0])
            fcol = next((c for c in cols if 'flux' in c.lower() or 'sb' in c.lower()), cols[1])
            self.wave_nm = np.asarray(t[wcol], float)
            self.flux = np.asarray(t[fcol], float)
            self._log('loaded %s: %s vs %s (%d points)' % (os.path.basename(path), wcol, fcol, len(self.wave_nm)))
        except Exception as e:
            self._log('ERROR loading spectrum: %s' % e)
            self.wave_nm = self.flux = None

    def _model(self, channel, name, deffile, mode):
        key = (channel, mode)
        if key not in self._models:
            m = spec.Spectrograph(name)
            with contextlib.redirect_stdout(io.StringIO()):     # hush the optical-surface printout
                m.build_model(deffile, throughput_mode=mode)
            self._models[key] = m
        return self._models[key]

    def compute(self):
        if self.wave_nm is None:
            self._log('no input spectrum loaded'); return
        texp = self.texp.value(); airmass = self.airmass.value()
        source = self.source.currentText(); mode = self.mode.currentText()
        ap_txt = self.aperture.currentText()
        aperture = ap_txt if ap_txt in ('optimal', 'single') else int(ap_txt)
        self._results = {}
        self._log('--- compute (%s, %s throughput) ---' % (source, mode))
        for name, specname, deffile, _ in CHANNELS:
            if not self.ch_boxes[name].isChecked():
                continue
            model = self._model(name, specname, deffile, mode)
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                if source == 'point':
                    counts, noise = observe.observe_spectrum(
                        model, texp, self.wave_nm, self.flux, airmass=airmass,
                        source='point', seeing=self.seeing.value(), aperture=aperture,
                        psf=self.psf.currentText())
                else:
                    counts, noise = observe.observe_spectrum(
                        model, texp, self.wave_nm, self.flux, airmass=airmass,
                        source='extended', nbin=self.nbin.value())
            self._results[name] = (model.waves, counts, noise)
            for line in buf.getvalue().strip().splitlines():
                self._log('[%s] %s' % (name, line.strip()))
        self._replot()

    def _replot(self):
        self.ax_counts.clear(); self.ax_snr.clear()
        for name, _, _, color in CHANNELS:
            if name not in self._results:
                continue
            waves, counts, noise = self._results[name]
            with np.errstate(invalid='ignore', divide='ignore'):
                snr = counts / noise
            self.ax_counts.plot(waves, counts, color=color, lw=0.8, label=name)
            self.ax_snr.plot(waves, snr, color=color, lw=0.8, label=name)
        self.ax_counts.set_ylabel('counts [e-]'); self.ax_counts.grid(alpha=0.3)
        self.ax_snr.set_ylabel('SNR / pixel'); self.ax_snr.set_xlabel('wavelength [nm]')
        self.ax_snr.grid(alpha=0.3)
        if self._results:
            self.ax_counts.legend(fontsize=8)
        self.fig.tight_layout()
        self.canvas.draw()

    def _save(self):
        if not self._results:
            self._log('nothing to save -- compute first'); return
        base, _ = QtWidgets.QFileDialog.getSaveFileName(self, 'Save results (basename)', SCRIPT_DIR)
        if not base:
            return
        base = os.path.splitext(base)[0]
        self.fig.savefig(base + '.png', dpi=130)
        with open(base + '.csv', 'w', newline='') as fh:
            wr = csv.writer(fh)
            wr.writerow(['channel', 'wave_nm', 'counts_e', 'noise_e', 'snr'])
            for name, (waves, counts, noise) in self._results.items():
                with np.errstate(invalid='ignore', divide='ignore'):
                    snr = counts / noise
                for i in range(len(waves)):
                    wr.writerow([name, '%.4f' % waves[i], '%.6g' % counts[i],
                                 '%.6g' % noise[i], '%.4f' % snr[i]])
        self._log('saved %s.png and %s.csv' % (base, base))

    def _log(self, msg):
        self.status.appendPlainText(msg)


def main():
    app = QtWidgets.QApplication(sys.argv)
    win = ETCWindow()
    win.resize(1100, 680)
    win.show()
    sys.exit(app.exec())


if __name__ == '__main__':
    main()
