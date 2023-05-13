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
        self.bpm_channels = config.bpm_channels[self.beamline]
        self.pulse_energy_ch = config.fast_xgm_ch[self.beamline]

    def read_ch(self, ch):
        if not self.dry_run:
            return read(ch)['data']

        if ch.endswith('KICK_MRAD.SP'):
            return (np.random.rand()-0.5)*1e-2

        if ch.endswith('ACTIVATE_FB'):
            return 0

        if 'BEAM_ENERGY_MEASUREMENT' in ch:
            return 16385.

        if 'BPM' in ch:
            return (np.random.rand(20)-0.5)*1e-6

        if 'INTENSITY.RAW' in ch:
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
        # Also need to obtain bpm names and positions
        return self.read_ch(self.bpm_channels)

    def read_pulse_energy(self):
        return self.read_ch(self.pulse_energy_ch)


