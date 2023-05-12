import numpy as np
import beamdynamics

eps_norm = 500e-9
energy_eV = 16.350e9
eps_geo = beamdynamics.eps_norm_to_geo(eps_norm, energy_eV)

r12 = 10
r22 = 2.5
beta = 15
alpha = -0.75

tg = beamdynamics.TwissGymnastics(beta, alpha, eps_geo, r12, r22)

a_arr = np.array([0.5, 1, 1.5])
phi_arr = np.linspace(-np.pi, np.pi, 11)[:10]

corr_angles = tg.Aphi_range_to_corrangles(a_arr, phi_arr)

