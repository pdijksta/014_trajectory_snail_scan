import sys
import os
import time

import PyQt5
import PyQt5.Qt
from PyQt5 import QtCore, QtGui
from PyQt5.QtWidgets import QMainWindow, QApplication

#import matplotlib.pyplot as plt
#from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT
import matplotlib

from gui import Ui_MainWindow
import config
import beamdynamics
import devices


matplotlib.use('Qt5Agg')


class Main(QMainWindow):
    def __init__(self):
        QMainWindow.__init__(self)
        self.ui = Ui_MainWindow()
        self.ui.setupUi(self)

        #Connect buttons
        self.ui.ResetButton.clicked.connect(self.initialize)
        self.ui.RestoreCorrectorsButton.clicked.connect(self.restore_correctors)

        #Initialize GUI
        self.ui.BeamlineSelect.addItems(config.beamlines)
        self.ui.StatusLabel.setText('Tool not initialized')

        self.disable_buttons = [
            self.ui.StartMeasurementButton,
            self.ui.RestoreCorrectorsButton,
            ]
        for button in self.disable_buttons:
            button.setEnabled(False)

    def initialize(self):
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
                'Correctors (init values)',
                '%s %.3f mrad' % (self.correctors[0], self.init_values[0]*1e3),
                '%s %.3f mrad' % (self.correctors[1], self.init_values[1]*1e3),
                ]
        self.ui.InformationLabel.setText('\n'.join(info))

    def restore_correctors(self):
        for corr, init_val in zip(self.correctors, self.init_values):
            self.mi.write_corrector(corr, init_val)


if __name__ == "__main__":

    if os.path.getmtime('./gui.ui') > os.path.getmtime('./gui.py'):
        cmd = 'bash ./ui2py.sh'
        print(cmd)
        os.system(cmd)
        time.sleep(1)

    # for pdb to work
    PyQt5.QtCore.pyqtRemoveInputHook()

    #make pyqt threadsafe
    QtCore.QCoreApplication.setAttribute(QtCore.Qt.AA_X11InitThreads)

    app = QApplication(sys.argv)
    app.setWindowIcon(QtGui.QIcon('./snails.png'))
    MainWindow = Main()
    MainWindow.show()
    sys.exit(app.exec_())

