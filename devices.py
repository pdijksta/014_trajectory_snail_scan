import numpy as np

import config

try:
    from pydoocs import read, write
    hasdoocs = True
except ImportError:
    print('Cannot import pydoocs')
    hasdoocs = False

class XFEL_interface:
    def __init__(self, dry_run, beamline):
        if not hasdoocs:
            dry_run = True
        self.dry_run = dry_run
        self.beamline = beamline

    def read_ch(self, ch):
        if not self.dry_run:
            return read(ch)['data']

        if ch.endswith('KICK_MRAD.SP'):
            return (np.random.rand()-0.5)*1e-2

        if ch.endswith('ACTIVATE_FB'):
            return 0

        if 'BEAM_ENERGY_MEASUREMENT' in ch:
            return 16385.

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
        ch = config.orbit_fb_channels[self.beamline]
        return self.read_ch(ch)

    def read_beam_energy_eV(self):
        ch = 'XFEL.DIAG/BEAM_ENERGY_MEASUREMENT/T4/ENERGY.SA%i' % int(self.beamline[-1])
        return self.read_ch(ch)*1e6





