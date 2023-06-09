import time
import numpy as np

from workers import WorkerBase

def read(parent, mi, measurement_time):
    this_pulse_energy = []
    this_orbit = []
    bpm_names = []
    time_begin = time_now = time.time()
    time_end = time_begin + measurement_time
    while time_now < time_end:
        if parent.abort:
            break
        bpm_names, orbit_vals = mi.read_orbit()
        pulse_energy_vals = mi.read_pulse_energy()

        this_orbit.append(orbit_vals)
        this_pulse_energy.append(pulse_energy_vals)

        time_prev = time_now
        time_now = time.time()
        delta = time_prev + 0.1 - time_now
        if parent.abort:
            break
        if delta > 0:
            time.sleep(delta)
            time_now += delta
    return bpm_names, np.array(this_orbit), np.array(this_pulse_energy)


class MeasWorker(WorkerBase):

    def func(self, dry_run, tg, mi, correctors, init_values, a_min, a_max, a_steps, phi_steps, settle_time, measurement_time):

        A_range, phi_range = tg.gen_Aphi_range(a_min, a_max, a_steps, phi_steps)
        delta_corr_arr, xxp_arr = tg.Aphi_range_to_corr(A_range, phi_range)

        pulse_ene_mean = np.full([len(A_range), len(phi_range)], np.nan, dtype=float)
        pulse_ene_std = pulse_ene_mean.copy()
        bpm_names, init_orbit, init_pulse_energy = read(self, mi, measurement_time)

        orbit_mean = np.full([len(A_range), len(phi_range), init_orbit.shape[1]], np.nan, dtype=float)
        orbit_std = orbit_mean.copy()
        ctr = 0
        n_meas = len(A_range)*len(phi_range)
        for n_A, A in enumerate(A_range):
            if self.abort:
                break
            for n_phi, phi in enumerate(phi_range):
                if self.abort:
                    break

                delta_corr = delta_corr_arr[n_A, n_phi]
                for corr, init, delta in zip(correctors, init_values, delta_corr):
                    mi.write_corrector(corr, init+delta)

                time_end = time.time() + settle_time
                while time.time() < time_end:
                    time.sleep(0.1)
                    if self.abort or dry_run:
                        break

                _, this_orbit, this_pulse_energy = read(self, mi, measurement_time)
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
                    'correctors': correctors,
                    'init_corrector_vals': init_values,
                    'beamline': mi.beamline,
                    'plane': mi.plane,
                    'emit_geo': tg.emit_geo,
                    'emit_norm': tg.emit_norm,
                    'energy_eV': tg.energy_eV,
                    'beta': tg.beta,
                    'alpha': tg.alpha,
                    'R12': tg.R12,
                    'R22': tg.R22,
                    },
                'data': {
                    'pulse_ene_mean': pulse_ene_mean,
                    'pulse_ene_std': pulse_ene_std,
                    'orbit_mean': orbit_mean,
                    'orbit_std': orbit_std,
                    'A': A_range,
                    'phi': phi_range,
                    'delta_corr_angles': delta_corr_arr,
                    'delta_xxp': xxp_arr,
                    'bpm_names': bpm_names,
                    'init_orbit': np.mean(init_orbit, axis=0),
                    'init_pulse_energy': np.mean(init_pulse_energy),
                    },
                }
        if self.abort:
            print('Measurement aborted')
        else:
            print('Measurement complete')
        return result_dict

