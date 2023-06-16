import sys
import os
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

import PyQt5
import PyQt5.Qt
from PyQt5 import QtCore, QtGui
from PyQt5.QtWidgets import QMainWindow, QApplication
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT

import measurement
import config
import beamdynamics
import devices
import plot_results
import logbook
import workers

if __name__ == '__main__' and (not os.path.isfile('./gui.py') or os.path.getmtime('./gui.ui') > os.path.getmtime('./gui.py')):
    cmd = 'bash ./ui2py.sh'
    print(cmd)
    os.system(cmd)

from gui import Ui_MainWindow


class Main(QMainWindow):

    data_dir = './data/'

    def __init__(self):
        QMainWindow.__init__(self)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        #For debug purposes
        #self.verbose = True

        #Connect buttons
        self.ui.ResetButton.clicked.connect(self.initialize)
        self.ui.RestoreCorrectorsButton.clicked.connect(self.restore_correctors)
        self.ui.CalcCorrAnglesButton.clicked.connect(self.calc_corr_angles)
        self.ui.CalcAPhiButton.clicked.connect(self.calc_A_phi)
        self.ui.StartMeasurementButton.clicked.connect(self.do_measurement)
        self.ui.AbortMeasurementButton.clicked.connect(self.abort_measurement)
        self.ui.LogbookButton.clicked.connect(self.do_logbook)
        self.ui.SaveAnalysisButton.clicked.connect(self.save_result)

        #Initialize some GUI widgets
        self.ui.BeamlineSelect.addItems(config.beamlines)
        self.ui.StatusLabel.setText('Tool not initialized')

        # Disable certain buttons on startup
        self.disable_buttons = [
            self.ui.StartMeasurementButton,
            self.ui.RestoreCorrectorsButton,
            self.ui.CalcCorrAnglesButton,
            self.ui.CalcAPhiButton,
            ]
        for button in self.disable_buttons:
            button.setEnabled(False)
        self.ui.AbortMeasurementButton.setEnabled(False)
        self.ui.LogbookButton.setEnabled(False)
        self.ui.SaveAnalysisButton.setEnabled(False)

        # Meas_lock must be false
        self.meas_lock = False

        # Plot tabs
        self.orbit_fig = None
        self.performance_fig = None

        # Always show first tab on startup
        self.ui.tabWidget.setCurrentIndex(0)

    def initialize(self):
        """
        Initializes self.mi, self.tg and sets init corrector values.
        """
        self.beamline = self.ui.BeamlineSelect.currentText()
        self.plane = self.ui.PlaneSelect.currentText()
        self.dry_run = self.ui.DryRunCheck.isChecked()

        self.mi = devices.XFEL_interface(self.dry_run, self.beamline, self.plane)
        self.correctors = config.corrector_names[self.beamline][self.plane]
        self.init_values = [self.mi.read_corrector(x) for x in self.correctors]

        energy_eV = self.mi.read_beam_energy_eV()

        eps_norm = self.ui.emit_norm_nm.value()*1e-9

        beta0, alpha0, mu0 = config.corrector_beta_alpha_mu[self.correctors[0]]
        beta1, alpha1, mu1 = config.corrector_beta_alpha_mu[self.correctors[1]]
        r12, r22 = beamdynamics.calc_r(beta0, alpha0, mu0, beta1, alpha1, mu1)
        self.tg = beamdynamics.TwissGymnastics(beta1, alpha1, eps_norm, energy_eV, r12, r22)

        self.ui.StatusLabel.setText('Tool initalized. DRY_RUN=%i' % self.dry_run)
        for button in self.disable_buttons:
            button.setEnabled(True)
        self.ui.LogbookButton.setEnabled(False)
        self.ui.SaveAnalysisButton.setEnabled(False)

        info = [
                'Information on tool status',
                'Correctors: init val (mrad); β (m); α',
                ]
        for n, (beta, alpha) in enumerate([
                (beta0, alpha0),
                (beta1, alpha1),
                ]):
            corr_info = 'c%i %s: %+.4f; %.1f; %+.2f' % (n, self.correctors[n], self.init_values[n]*1e3, beta, alpha)
            info.append(corr_info)
        self.ui.InformationLabel.setText('\n'.join(info))
        self.ui.InformationLabel.setFont(QtGui.QFont('Monospace'))

        self.ui.CorrAngleResult.setText('Result:')
        self.ui.A_Phi_Result.setText('Result:')
        self.clear_plot_tabs()

    def restore_correctors(self):
        for corr, init_val in zip(self.correctors, self.init_values):
            self.mi.write_corrector(corr, init_val)

    def calc_corr_angles(self):
        A = self.ui.A_calc.value()
        phi = self.ui.Phi_calc.value()/180*np.pi
        x, xp = self.tg.Aphi_to_trajoffset(A, phi)
        c0, c1 = self.tg.trajoffset_to_corr(x, xp)
        self.ui.CorrAngleResult.setText('A=%.3f, ψ=%.1f deg --> Δc0=%.5f mrad, Δc1=%.5f mrad' % (A, phi*180/np.pi, c0*1e3, c1*1e3))

    def calc_A_phi(self):
        c0 = self.ui.C0_calc.value()/1e3
        c1 = self.ui.C1_calc.value()/1e3
        x, xp = self.tg.corr_to_trajoffset(c0, c1)
        A, phi = self.tg.trajoffset_to_Aphi(x, xp)
        self.ui.A_Phi_Result.setText('Δc0=%.5f mrad, Δc1=%.5f mrad --> A=%.3f, ψ=%.1f deg' % (c0*1e3, c1*1e3, A, phi*180/np.pi))

    def post_measurement(self):
        self.result_dict = self.meas_worker.outp
        self.new_figures()
        plot_results.plot_Aphi_scan(self.result_dict, plot_handles=self.performance_plot_handles)
        plot_results.plot_orbit(self.result_dict, plot_handles=self.orbit_plot_handles)
        self.ui.tabWidget.setCurrentIndex(1)
        self.ui.LogbookButton.setEnabled(True)
        self.ui.SaveAnalysisButton.setEnabled(True)

    def measurement_progress(self, val):
        self.ui.progressBar.setValue(val)

    def do_measurement(self):
        if self.mi.read_orbit() == 1:
            raise ValueError('Orbit feedback is active!')
        args = (
                self.dry_run,
                self.tg,
                self.mi,
                self.correctors,
                self.init_values,
                self.ui.A_min.value(),
                self.ui.A_max.value(),
                self.ui.A_steps.value(),
                self.ui.Phi_steps.value(),
                self.ui.SettleTime.value(),
                self.ui.MeasurementTime.value(),
                )
        start_funcs = (self.lock_meas, )
        finish_funcs = (self.unlock_meas, self.restore_correctors, self.post_measurement, )
        progress_funcs = (self.measurement_progress, )
        if self.meas_lock:
            raise RuntimeError('Cannot start new analysis while lock is active')
        self.meas_thread, self.meas_worker = workers.threaded_func(self, measurement.MeasWorker, args, {}, start_funcs, finish_funcs, progress_funcs)

    def abort_measurement(self):
        self.meas_worker.abort = True

    def lock_meas(self):
        self.ui.ResetButton.setEnabled(False)
        self.ui.StartMeasurementButton.setEnabled(False)
        self.ui.AbortMeasurementButton.setEnabled(True)
        self.meas_lock = True

    def unlock_meas(self):
        self.ui.ResetButton.setEnabled(True)
        self.ui.StartMeasurementButton.setEnabled(True)
        self.ui.AbortMeasurementButton.setEnabled(False)
        self.meas_lock = False

    def clear_plot_tabs(self):
        for fig, layout in [
                (self.orbit_fig, self.ui.OrbitLayout),
                (self.performance_fig, self.ui.PerformanceLayout),
                ]:
            if fig is not None:
                fig.clf()
                plt.close(fig)
            if layout is not None:
                for i in reversed(range(layout.count())):
                    layout.itemAt(i).widget().deleteLater()

    def new_figures(self):
        self.clear_plot_tabs()
        correctors = self.result_dict['input']['correctors']
        minA = self.result_dict['data']['A'].min()
        plane = self.result_dict['input']['plane']
        self.performance_plot_handles = plot_results.performance_figure(correctors[1], plane)
        canvas = FigureCanvasQTAgg(self.performance_plot_handles[0])
        toolbar = NavigationToolbar2QT(canvas, self)
        self.ui.PerformanceLayout.addWidget(canvas)
        self.ui.PerformanceLayout.addWidget(toolbar)

        self.orbit_plot_handles = plot_results.orbit_figure(correctors, minA, plane)
        canvas = FigureCanvasQTAgg(self.orbit_plot_handles[0])
        toolbar = NavigationToolbar2QT(canvas, self)
        self.ui.OrbitLayout.addWidget(canvas)
        self.ui.OrbitLayout.addWidget(toolbar)

    def do_logbook(self):
        index = self.ui.tabWidget.currentIndex()
        self.ui.tabWidget.setCurrentIndex(1)
        logbook.log_screen(self.save_result, self)
        self.ui.tabWidget.setCurrentIndex(index)

    def save_result(self):
        ok, comment = logbook.dialog(self)
        data = self.result_dict
        data2 = {}
        for key1, subdict1 in data.items():
            for key2, val in subdict1.items():
                data2['%s/%s' % (key1, key2)] = val

        Path(self.data_dir).mkdir(parents=True, exist_ok=True)
        filename = self.data_dir + time.strftime("%Y%m%d-%H_%M_%S") + "_snail_scan.npz"
        np.savez(filename, **data2)
        print('Saved analysis under %s' % os.path.abspath(filename))
        return filename, comment, ok

if __name__ == "__main__":
    # for pdb to work
    PyQt5.QtCore.pyqtRemoveInputHook()

    #make pyqt threadsafe
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_X11InitThreads)

    app = QApplication(sys.argv)
    app.setWindowIcon(QtGui.QIcon('./snails.png'))
    MainWindow = Main()
    MainWindow.show()
    sys.exit(app.exec_())

