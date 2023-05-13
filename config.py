
beamlines = ['SASE1', 'SASE2', 'SASE3']

corrector_names = {
        'SASE1': {
            'X': ['CFX.2133.T2', 'CFX.2162.T2'],
            'Y': ['CFY.2146.T2', 'CFY.2177.T2'],
            },
        'SASE2': {
            'X': ['CFX.2154.T1', 'CEX.2196.T1'],
            'Y': ['CFY.2168.T1', 'CNY.2196.T1'],
            },
        'SASE3': {
            'X': ['CFX.2762.T4', 'CEX.2795.T4'],
            'Y': ['CFY.2774.T4', 'CNY.2799.T4'],
            },
        }

# From "Optics calculation server" on 2023-05-12 afternoon
corrector_beta_alpha_mu = {
        #SASE1
        'CFX.2133.T2': [50.43, 2.33, 123.239],
        'CFY.2146.T2': [48.55, 2.06, 127.69035],
        'CFX.2162.T2': [48.66, 2.76, 124.4975],
        'CFY.2177.T2': [57.87, 2.13, 129.10465],

        #SASE2
        'CFX.2154.T1': [60.48, 1.77, 124.05925],
        'CFY.2168.T1': [44.09, 1.59, 129.14872],
        'CEX.2196.T1': [30.87, -0.83, 125.2908],
        'CNY.2196.T1': [28.50, 0.79, 130.17615],

        #SASE3
        'CFX.2762.T4': [60.29, 2.65, 148.56375],
        'CFY.2774.T4': [41.92, 2.16, 152.232],
        'CEX.2795.T4': [8.98, -0.54, 150.13044],
        'CNY.2799.T4': [10.99, 0.96, 153.7058],
        }

fast_xgm_ch = {
        'SASE1': 'XFEL.FEL/XGM/XGM.2643.T9/INTENSITY.SA1.RAW.TRAIN',
        'SASE2': 'XFEL.FEL/XGM/XGM.2595.T6/INTENSITY.RAW.TRAIN',
        'SASE3': 'XFEL.FEL/XGM/XGM.3130.T10/INTENSITY.SA3.RAW.TRAIN',
        }

orbit_fb_channels = {
        'SASE1': 'XFEL.FEEDBACK/ORBIT.SA1/ORBITFEEDBACK.ACTIVATE_FB',
        'SASE2': 'XFEL.FEEDBACK/ORBIT.SA2/ORBITFEEDBACK.ACTIVATE_FB',
        'SASE3': 'XFEL.FEEDBACK/ORBIT.SA3/ORBITFEEDBACK.ACTIVATE_FB',
        }

bpm_chs = lambda beamline, plane: 'XFEL.DIAG/BPM/*.SA%i/%s.ALL' % (int(beamline[-1]), plane)

energy_ch = lambda beamline: 'XFEL.DIAG/BEAM_ENERGY_MEASUREMENT/CL/ENERGY.SA%i' % int(beamline[-1])

