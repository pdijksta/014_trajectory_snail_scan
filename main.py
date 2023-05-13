import sys
import os
import time

import numpy as np
import matplotlib

import PyQt5
import PyQt5.Qt
from PyQt5 import QtCore, QtGui
from PyQt5.QtWidgets import QMainWindow, QApplication

#import matplotlib.pyplot as plt
#from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT

import config
import beamdynamics
import devices

if __name__ == '__main__' and (not os.path.isfile('./gui.py') or os.path.getmtime('./gui.ui') > os.path.getmtime('./gui.py')):
    cmd = 'bash ./ui2py.sh'
    print(cmd)
    os.system(cmd)

from gui import Ui_MainWindow

matplotlib.use('Qt5Agg')

class Main(QMainWindow):
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
        self.ui.StartMeasurementButton.clicked.connect(self.measurement)

        #Initialize GUI
        self.ui.BeamlineSelect.addItems(config.beamlines)
        self.ui.StatusLabel.setText('Tool not initialized')

        self.disable_buttons = [
            self.ui.StartMeasurementButton,
            self.ui.RestoreCorrectorsButton,
            self.ui.CalcCorrAnglesButton,
            self.ui.CalcAPhiButton
            ]
        for button in self.disable_buttons:
            button.setEnabled(False)

        # Always show first tab on startup
        self.ui.tabWidget.setCurrentIndex(0)

    def initialize(self):
        """
        Initializes self.mi, self.tg and sets init corrector values.
        """
        self.beamline = self.ui.BeamlineSelect.currentText()
        self.plane = self.ui.PlaneSelect.currentText()
        self.dry_run = self.ui.DryRunCheck.isChecked()

        self.mi = devices.XFEL_interface(self.dry_run, self.beamline)
        self.correctors = config.corrector_names[self.beamline][self.plane]
        self.init_values = [self.mi.read_corrector(x) for x in self.correctors]

        energy_eV = self.mi.read_beam_energy_eV()

        eps_norm = self.ui.emit_norm_nm.value()*1e-9

        eps_geo = beamdynamics.eps_norm_to_geo(eps_norm, energy_eV)
        beta0, alpha0, mu0 = config.corrector_beta_alpha_mu[self.correctors[0]]
        beta1, alpha1, mu1 = config.corrector_beta_alpha_mu[self.correctors[1]]
        r12, r22 = beamdynamics.calc_r(beta0, alpha0, mu0, beta1, alpha1, mu1)
        self.tg = beamdynamics.TwissGymnastics(beta1, alpha1, eps_geo, r12, r22)

        self.ui.StatusLabel.setText('Tool initalized. DRY_RUN=%i' % self.dry_run)
        for button in self.disable_buttons:
            button.setEnabled(True)

        info = [
                'Information on tool status',
                'Correctors: init value (mrad); β (m); α',
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

    def measurement(self):
        a_min = self.ui.A_min.value()
        a_max = self.ui.A_max.value()
        a_steps = self.ui.A_steps.value()
        phi_steps = self.ui.Phi_steps.value()
        settle_time = self.ui.SettleTime.value()
        measurement_time = self.ui.MeasurementTime.value()

        A_range, phi_range = self.tg.gen_Aphi_range(a_min, a_max, a_steps, phi_steps)
        delta_corr_arr = self.tg.Aphi_range_to_corr(A_range, phi_range)

        pulse_ene_mean = np.full([len(A_range), len(phi_range)], np.nan, dtype=float)
        pulse_ene_std = pulse_ene_mean.copy()

        first = True

        n_meas = len(A_range)*len(phi_range)

        ctr = 0
        for n_A, A in enumerate(A_range):
            for n_phi, phi in enumerate(phi_range):

                this_pulse_energy = []
                this_orbit = []

                delta_corr = delta_corr_arr[n_A, n_phi]
                for corr, init, delta in zip(self.correctors, self.init_values, delta_corr):
                    self.mi.write_corrector(corr, init+delta)

                if not self.dry_run:
                    time.sleep(settle_time)

                time_begin = time_now = time.time()
                time_end = time_begin + measurement_time
                while time_now < time_end:
                    orbit_vals = self.mi.read_orbit()
                    pulse_energy_vals = self.mi.read_pulse_energy()
                    if first:
                        orbit_mean = np.full([len(A_range), len(phi_range), len(orbit_vals)], np.nan, dtype=float)
                        orbit_std = orbit_mean.copy()
                        first = False

                    this_orbit.append(orbit_vals)
                    this_pulse_energy.append(pulse_energy_vals)

                    time_prev = time_now
                    time_now = time.time()
                    delta = time_prev + 0.1 - time_now
                    if delta > 0:
                        time.sleep(delta)
                        time_now += delta

                this_pulse_energy = np.array(this_pulse_energy)
                this_orbit = np.array(this_orbit)

                pulse_ene_mean[n_A, n_phi] = np.mean(this_pulse_energy)
                pulse_ene_std[n_A, n_phi] = np.std(this_pulse_energy)
                orbit_mean[n_A, n_phi] = np.mean(this_orbit, axis=0)
                orbit_std[n_A, n_phi] = np.std(this_orbit, axis=0)

                ctr += 1
                print('Measurement %i out of %i done.' % (ctr, n_meas))

        result_dict = {
                'input': {
                    'A_min_max_steps': [a_min, a_max, a_steps],
                    'phi_steps': phi_steps,
                    'settle_time': settle_time,
                    'measurement_time': measurement_time,
                    'correctors': self.correctors,
                    'init_corrector_vals': self.init_values,
                    'beamline': self.beamline,
                    'plane': self.plane,
                    },
                'data': {
                    'pulse_ene_mean': pulse_ene_mean,
                    'pulse_ene_std': pulse_ene_std,
                    'orbit_mean': orbit_mean,
                    'orbit_std': orbit_std,
                    'A': A_range,
                    'phi': phi_range,
                    'delta_corr_angles': delta_corr_arr,
                    },
                }
        return result_dict


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

