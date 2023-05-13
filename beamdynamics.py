import numpy as np
from scipy.constants import m_e, c, e

m_e_eV = m_e*c**2/e

def eps_norm_to_geo(emit_norm, energy_eV):
    return emit_norm/(energy_eV/m_e_eV)

def calc_r(beta0, alpha0, mu0, beta1, alpha1, mu1):
    delta_mu = mu1 - mu0
    r12 = np.sqrt(beta1*beta0) * np.sin(delta_mu)
    r22 = np.sqrt(beta0/beta1) * (np.cos(delta_mu) - alpha1*np.sin(delta_mu))
    return r12, r22

class TwissGymnastics:
    def __init__(self, beta, alpha, emit_norm, energy_eV, R12, R22):
        self.beta = beta
        self.alpha = alpha
        self.gamma = (1+self.alpha**2)/self.beta
        self.emit_norm = emit_norm
        self.energy_eV = energy_eV
        self.emit_geo = eps_norm_to_geo(emit_norm, energy_eV)
        self.R12 = R12
        self.R22 = R22
        assert R12 != 0

    def A(self, x, xp):
        return 1/self.emit_geo * (self.gamma*x**2 + 2*self.alpha*x*xp + self.beta*xp**2)

    def rescaled_xxp(self, x, xp):
        x2 = x/np.sqrt(self.beta)
        xp2 = (self.alpha*x + self.beta*xp)/np.sqrt(self.beta)
        return x2, xp2

    def trajoffset_to_Aphi(self, x, xp):
        A = self.A(x, xp)
        x2, xp2 = self.rescaled_xxp(x, xp)
        phi = np.arctan2(xp2, x2)
        return A, phi

    def Aphi_to_trajoffset(self, A, phi):
        gradient = (np.tan(phi) - self.alpha)/self.beta
        x = np.sqrt(self.emit_geo*A/(self.gamma + 2*self.alpha*gradient + self.beta*gradient**2))
        if type(x) == np.ndarray:
            x[np.abs(phi) > np.pi/2] *= -1
        elif abs(x) > np.pi/2:
            x *= -1
        xp = gradient*x
        return x, xp

    def trajoffset_to_corr(self, x, xp):
        angle0 = x/self.R12
        angle1 = -self.R22/self.R12 * x + xp
        return angle0, angle1

    def corr_to_trajoffset(self, angle0, angle1):
        x = self.R12*angle0
        xp = self.R22*angle0 + angle1
        return x, xp

    def gen_Aphi_range(self, A_min, A_max, A_points, phi_points):
        if A_points == 1:
            A_arr = np.array([A_min])
        else:
            A_arr = np.linspace(A_min, A_max, A_points)
        phi_arr = np.linspace(-np.pi, np.pi, phi_points+1)[:-1]
        return A_arr, phi_arr

    def Aphi_range_to_corr(self, A_arr, phi_arr):
        corr_angles = np.zeros([A_arr.size, phi_arr.size, 2])
        xxp_arr = corr_angles.copy()
        for n_A, A in enumerate(A_arr):
            x_arr, xp_arr = self.Aphi_to_trajoffset(A, phi_arr)
            angle0, angle1 = self.trajoffset_to_corr(x_arr, xp_arr)
            corr_angles[n_A,:,0] = angle0
            corr_angles[n_A,:,1] = angle1
            xxp_arr[n_A,:,0] = x_arr
            xxp_arr[n_A,:,1] = xp_arr
        return corr_angles, xxp_arr

