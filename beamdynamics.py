import numpy as np
from scipy.constants import m_e, c, e

def eps_norm_to_geo(eps_norm, energy_eV):
    m_e_eV = m_e*c**2/e
    return eps_norm/(energy_eV/m_e_eV)

class TwissGymnastics:
    def __init__(self, beta, alpha, eps_geo, R12, R22):
        self.beta = beta
        self.alpha = alpha
        self.gamma = (1+self.alpha**2)/self.beta
        self.eps_geo = eps_geo
        self.R12 = R12
        self.R22 = R22
        assert R12 != 0

    def A(self, x, xp):
        return 1/self.eps_geo * (self.gamma*x**2 + 2*self.alpha*x*xp + self.beta*xp**2)

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
        x = np.sqrt(self.eps_geo*A/(self.gamma + 2*self.alpha*gradient + self.beta*gradient**2))
        if type(x) == np.ndarray:
            x[np.abs(phi) > np.pi/2] *= -1
        elif abs(x) > np.pi/2:
            x *= -1
        xp = gradient*x
        return x, xp

    def trajoffset_to_corrangle(self, x, xp):
        angle1 = x/self.R12
        angle2 = -self.R22/self.R12 * x + xp
        return angle1, angle2

    def corrangle_to_trajoffset(self, angle1, angle2):
        x = self.R12*angle1
        xp = self.R22*angle1 + angle2
        return x, xp

    def gen_Aphi_range(self, A_min, A_max, A_points, phi_points):
        if A_points == 1:
            A_arr = np.array([A_min])
        else:
            A_arr = np.linspace(A_min, A_max, A_points)
        phi_arr = np.linspace(-np.pi, np.pi, phi_points+1)[:-1]
        return A_arr, phi_arr

    def Aphi_range_to_corrangles(self, A_arr, phi_arr):
        corr_angles = np.zeros([A_arr.size, phi_arr.size, 2])
        for n_A, A in enumerate(A_arr):
            x_arr, xp_arr = self.Aphi_to_trajoffset(A, phi_arr)
            angle1, angle2 = self.trajoffset_to_corrangle(x_arr, xp_arr)
            corr_angles[n_A,:,0] = angle1
            corr_angles[n_A,:,1] = angle2
        return corr_angles

