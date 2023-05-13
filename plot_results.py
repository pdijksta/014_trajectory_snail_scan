import numpy as np
import matplotlib.pyplot as plt

def phase_space_ellipse(emittance, beta, alpha, n_points):

    # - parametrize an ellipse using emittance, beta, alpha
    # - evaluate this parametrization for n_points

    gamma = (1 + alpha**2)/beta

    if gamma == beta:
        tan_2phi = 0
    else:
        tan_2phi = 2*alpha / (gamma - beta)
    phi = np.arctan(tan_2phi)/2

    # Find out length of ellipse axes
    x_sq = emittance/gamma
    y_sq = emittance/beta
    sin_sq, cos_sq = np.sin(phi)**2, np.cos(phi)**2
    b4 = (-sin_sq * x_sq + cos_sq * y_sq) / (cos_sq * x_sq - sin_sq * y_sq) * emittance**2
    b = b4**(1/4)
    a = emittance/b

    # Straight ellipse
    t = np.linspace(0, 2*np.pi, n_points)
    x = a * np.cos(t)
    y = b * np.sin(t)

    # Rotate it
    x2 = np.cos(phi) * x - np.sin(phi) * y
    y2 = np.cos(phi) * y + np.sin(phi) * x

    return x2, y2


def plot_Aphi_scan(result_dict, plot_handles=None):
    if plot_handles is None:
        fig = plt.figure(figsize=(12,10))
        sp_ellipse = plt.subplot(2, 2, 1)
        sp_ellipse.set_title('Phase space at %s' % result_dict['input']['correctors'][1])
        sp_ellipse.set_xlabel('$\Delta x$ ($\mu$m)')
        sp_ellipse.set_ylabel('$\Delta x\'$ ($\mu$rad)')

        sp_ellipse_norm = plt.subplot(2, 2, 2)
        sp_ellipse_norm.set_title('Norm. phase space at %s' % result_dict['input']['correctors'][1])
        sp_ellipse_norm.set_xlabel(r'$\Delta x/\sqrt{\epsilon\beta}$')
        sp_ellipse_norm.set_ylabel(r'$(\alpha\Delta x +\Delta x\'\beta)/\sqrt{\epsilon\beta}$')

    else:
        fig, (sp_ellipse, ) = plot_handles

    emit_geo = result_dict['input']['emit_geo']
    beta = result_dict['input']['beta']
    alpha = result_dict['input']['alpha']

    def rescale_xxp(x, xp):
        x_out = x/np.sqrt(emit_geo*beta)
        xp_out = (x*alpha + xp*beta)/np.sqrt(emit_geo*beta)
        return x_out, xp_out


    x, xp = phase_space_ellipse(emit_geo, beta, alpha, 1000)
    xs, xps = rescale_xxp(x, xp)

    sp_ellipse.plot(x, xp, color='black')
    sp_ellipse_norm.plot(xs, xps, color='black')

    for n_A, A in enumerate(result_dict['data']['A']):
        x = result_dict['data']['delta_xxp'][n_A,:,0]
        xp = result_dict['data']['delta_xxp'][n_A,:,1]
        sp_ellipse.plot(x, xp, ls='None', marker='.', label=A)

        xs, xps = rescale_xxp(x, xp)
        sp_ellipse_norm.plot(xs, xps, ls='None', marker='.', label=A)


    sp_ellipse.legend(title='A')





