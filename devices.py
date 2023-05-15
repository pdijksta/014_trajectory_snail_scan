import numpy as np
import config

try:
    from pydoocs import read, write
    hasdoocs = True
except ImportError:
    print('Cannot import pydoocs')
    hasdoocs = False

class XFEL_interface:
    def __init__(self, dry_run, beamline, plane):
        if not hasdoocs:
            dry_run = True
        self.dry_run = dry_run
        self.beamline = beamline
        self.plane = plane
        self.orbit_fb_channel = config.orbit_fb_channels[self.beamline]
        self.energy_ch = config.energy_ch(self.beamline)
        self.bpm_channels = config.bpm_chs(self.beamline, self.plane)
        self.pulse_energy_ch = config.fast_xgm_ch[self.beamline]

    def read_ch(self, ch):
        if self.dry_run and ch.endswith('ACTIVATE_FB'):
            return 0

        if hasdoocs:
            return read(ch)['data']

        if ch.endswith('KICK_MRAD.SP'):
            return (np.random.rand()-0.5)*1e-2

        if 'BEAM_ENERGY_MEASUREMENT' in ch:
            return 16385.

        if 'BPM' in ch:
            outp = []
            for ctr in range(20):
                outp.append([0, np.random.rand(), np.random.rand(), 0, 'BPME.%04i.SA%i' % (ctr, int(self.beamline[-1]))])
            return outp

        if 'RAW.TRAIN' in ch:
            return 500 + np.random.rand()*50

        raise ValueError(ch)

    def write_ch(self, ch, val):
        if self.dry_run:
            print('I would write %.3e to %s' % (val, ch))
        else:
            write(ch, val)

    def read_corrector(self, device_name):
        ch = 'XFEL.MAGNETS/MAGNET.ML/%s/KICK_MRAD.SP' % device_name
        return self.read_ch(ch)/1e3

    def write_corrector(self, device_name, val):
        ch = 'XFEL.MAGNETS/MAGNET.ML/%s/KICK_MRAD.SP' % device_name
        self.write_ch(ch, val*1e3)

    def read_fb_status(self):
        return self.read_ch(self.orbit_fb_channel)

    def read_beam_energy_eV(self):
        return self.read_ch(self.energy_ch)*1e6

    def read_orbit(self):
        bpm_data = self.read_ch(self.bpm_channels)
        bpm_names = [x[-1] for x in bpm_data]
        bpm_vals = np.array([x[1] for x in bpm_data])
        return bpm_names, bpm_vals*1e-3

    def read_pulse_energy(self):
        return self.read_ch(self.pulse_energy_ch)/1e6

