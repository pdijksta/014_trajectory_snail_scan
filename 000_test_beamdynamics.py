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

corr_angles, _ = tg.Aphi_range_to_corr(a_arr, phi_arr)


corrs0, corrs1 = np.random.rand(100), np.random.rand(100)

x, xp = tg.corr_to_trajoffset(corrs0, corrs1)
A, phi = tg.trajoffset_to_Aphi(x, xp)

x2, xp2, = tg.Aphi_to_trajoffset(A, phi)
corrs0_2, corrs1_2 = tg.trajoffset_to_corr(x2, xp2)

assert np.all(np.abs(corrs0 - corrs0_2) < 1e-10)
assert np.all(np.abs(corrs1 - corrs1_2) < 1e-10)


