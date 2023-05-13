import sys
import os
import time
from pathlib import Path

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

import PyQt5
import PyQt5.Qt
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import QThread
from PyQt5.QtWidgets import QMainWindow, QApplication, QInputDialog

#import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT

import measurement
import config
import beamdynamics
import devices
import plot_results
import logbook

if __name__ == '__main__' and (not os.path.isfile('./gui.py') or os.path.getmtime('./gui.ui') > os.path.getmtime('./gui.py')):
    cmd = 'bash ./ui2py.sh'
    print(cmd)
    os.system(cmd)

from gui import Ui_MainWindow

matplotlib.use('Qt5Agg')

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

        #Initialize some GUI widgets
        self.ui.BeamlineSelect.addItems(config.beamlines)
        self.ui.StatusLabel.setText('Tool not initialized')

        # Disable certain buttons on startup
        self.disable_buttons = [
            self.ui.StartMeasurementButton,
            self.ui.RestoreCorrectorsButton,
            self.ui.CalcCorrAnglesButton,
            self.ui.CalcAPhiButton,
            self.ui.LogbookButton,
            ]
        for button in self.disable_buttons:
            button.setEnabled(False)
        self.ui.AbortMeasurementButton.setEnabled(False)

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

        info = [
                'Information on tool status',
                'Correctors: init val (mrad); β (m); α',
                ]
        for n, (beta, alpha) in enumerate([
                (beta0, alpha0),
                (beta1, alpha1),
                ]):
            corr_info = 'c%i %s: %+.4f; %.1f; %+.2f' % (n, self.correctors[0], self.init_values[n]*1e3, beta, alpha)
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
        self.result_dict = self.func_worker.outp
        self.restore_correctors()
        self.new_figures()
        plot_results.plot_Aphi_scan(self.result_dict, plot_handles=self.performance_plot_handles)
        self.ui.tabWidget.setCurrentIndex(1)

    def measurement_progress(self, val):
        self.ui.progressBar.setValue(val)

    def threaded_func(self, worker, args, kwargs, post_func):
        # Written using example under https://realpython.com/python-pyqt-qthread/
        if self.meas_lock:
            raise RuntimeError('Cannot start new analysis while lock is active')

        self.func_thread = QThread(parent=self)
        self.func_worker = worker(*args, **kwargs)
        self.func_worker.moveToThread(self.func_thread)

        self.func_thread.started.connect(self.lock_meas)
        self.func_thread.started.connect(self.func_worker.run)
        self.func_thread.finished.connect(self.unlock_meas)

        self.func_worker.progress.connect(self.measurement_progress)

        self.func_worker.finished.connect(post_func)
        self.func_worker.finished.connect(self.func_thread.quit)
        self.func_worker.finished.connect(self.func_thread.wait)
        self.func_worker.finished.connect(self.func_thread.deleteLater)

        self.func_thread.start()

    def do_measurement(self):
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
        self.threaded_func(measurement.MeasWorker, args, {}, self.post_measurement)

    def abort_measurement(self):
        self.func_worker.abort = True

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
        rec_point = self.correctors[1]
        self.performance_plot_handles = plot_results.performance_figure(rec_point)
        canvas = FigureCanvasQTAgg(self.performance_plot_handles[0])
        toolbar = NavigationToolbar2QT(canvas, self)
        self.ui.PerformanceLayout.addWidget(canvas)
        self.ui.PerformanceLayout.addWidget(toolbar)

    def logbook(self, widget, text=""):
        screenshot = self.get_screenshot(widget)

        res = logbook.send_to_desy_elog(author='Dr. Snail', title='OrbitSnailScan', severity='INFO', text=text, elog='xfellog', image=screenshot)
        if not res:
            print('error during eLogBook sending')

    def get_screenshot(self, window_widget):
        screenshot_tmp = QtCore.QByteArray()
        screeshot_buffer = QtCore.QBuffer(screenshot_tmp)
        screeshot_buffer.open(QtCore.QIODevice.WriteOnly)
        widget = QtWidgets.QWidget.grab(window_widget)
        widget.save(screeshot_buffer, "png")
        return screenshot_tmp.toBase64().data().decode()

    def log_screen(self, widget, auto_comment=""):
        dlg = QInputDialog(self)
        dlg.setInputMode(QInputDialog.TextInput)
        dlg.setLabelText("Comment :")
        dlg.resize(400, 100)
        ok = dlg.exec_()
        comment = dlg.textValue()

        filename = self.save_result()
        text = 'Data is saved in %s' % filename
        if ok:
            text = comment + "\n" +"\n" + text
        text = auto_comment + text
        self.logbook(widget, text=text)

    def do_logbook(self):
        self.ui.tabWidget.setCurrentIndex(1)
        self.log_screen(self)

    def save_result(self):
        data = self.result_dict
        data2 = {}
        for key1, subdict1 in data.items():
            for key2, val in subdict1.items():
                data2['%s/%s' % (key1, key2)] = val

        Path(self.data_dir).mkdir(parents=True, exist_ok=True)
        filename = self.data_dir + time.strftime("%Y%m%d-%H_%M_%S") + "_snail_scan.npz"
        np.savez(filename, **data2)
        return filename

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

