import time
import numpy as np

from PyQt5.QtCore import QObject, pyqtSignal

class WorkerBase(QObject):
    """
    Wraps a function call such that it can be used by QThread.
    """
    finished = pyqtSignal()
    progress = pyqtSignal(int)

    def __init__(self, *args, **kwargs):
        QObject.__init__(self)
        self.args = args
        self.kwargs = kwargs
        self.outp = None
        self.error = False

    def run(self):
        try:
            self.outp = self.func(*self.args, **self.kwargs)
        except Exception as e:
            self.error = True
            print('%s failed with following error message:\n%s' % (self.function, e))
        finally:
            self.finished.emit()

class MeasWorker(WorkerBase):

    def func(self, parent, a_min, a_max, a_steps, phi_steps, settle_time, measurement_time):

        A_range, phi_range = parent.tg.gen_Aphi_range(a_min, a_max, a_steps, phi_steps)
        delta_corr_arr = parent.tg.Aphi_range_to_corr(A_range, phi_range)

        pulse_ene_mean = np.full([len(A_range), len(phi_range)], np.nan, dtype=float)
        pulse_ene_std = pulse_ene_mean.copy()

        n_meas = len(A_range)*len(phi_range)

        ctr = 0
        for n_A, A in enumerate(A_range):
            for n_phi, phi in enumerate(phi_range):

                this_pulse_energy = []
                this_orbit = []

                delta_corr = delta_corr_arr[n_A, n_phi]
                for corr, init, delta in zip(parent.correctors, parent.init_values, delta_corr):
                    parent.mi.write_corrector(corr, init+delta)

                if not parent.dry_run:
                    time.sleep(settle_time)

                time_begin = time_now = time.time()
                time_end = time_begin + measurement_time
                while time_now < time_end:
                    orbit_vals = parent.mi.read_orbit()
                    pulse_energy_vals = parent.mi.read_pulse_energy()
                    if ctr == 0:
                        orbit_mean = np.full([len(A_range), len(phi_range), len(orbit_vals)], np.nan, dtype=float)
                        orbit_std = orbit_mean.copy()

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
                self.progress.emit((ctr*100)//n_meas)

        result_dict = {
                'input': {
                    'A_min_max_steps': [a_min, a_max, a_steps],
                    'phi_steps': phi_steps,
                    'settle_time': settle_time,
                    'measurement_time': measurement_time,
                    'correctors': parent.correctors,
                    'init_corrector_vals': parent.init_values,
                    'beamline': parent.beamline,
                    'plane': parent.plane,
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


