import sys
import os

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

if __name__ == '__main__' and os.path.getmtime('./gui.ui') > os.path.getmtime('./gui.py'):
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
        beamline = self.ui.BeamlineSelect.currentText()
        plane = self.ui.PlaneSelect.currentText()
        dry_run = self.ui.DryRunCheck.isChecked()

        self.mi = devices.XFEL_interface(dry_run, beamline)
        self.correctors = config.corrector_names[beamline][plane]
        self.init_values = [self.mi.read_corrector(x) for x in self.correctors]

        energy_eV = self.mi.read_beam_energy_eV()

        eps_norm = self.ui.emit_norm_nm.value()*1e-9

        eps_geo = beamdynamics.eps_norm_to_geo(eps_norm, energy_eV)
        beta0, alpha0, mu0 = config.corrector_beta_alpha_mu[self.correctors[0]]
        beta1, alpha1, mu1 = config.corrector_beta_alpha_mu[self.correctors[1]]
        r12, r22 = beamdynamics.calc_r(beta0, alpha0, mu0, beta1, alpha1, mu1)
        self.tg = beamdynamics.TwissGymnastics(beta1, alpha1, eps_geo, r12, r22)

        self.ui.StatusLabel.setText('Tool initalized. DRY_RUN=%i' % self.mi.dry_run)
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

        A_range, phi_range = self.tg.gen_Aphi_range(a_min, a_max, a_steps, phi_steps)
        delta_corr_arr = self.tg.Aphi_range_to_corr(A_range, phi_range)


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

