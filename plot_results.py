from socket import gethostname
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import Normalize

LEFT = 0.1
RIGHT = 0.9
HSPACE = 0.4
WSPACE = 0.4
TOP = None
BOTTOM = None

if 'mpyubl38552' in gethostname():
    SMALL_SIZE = 7
    MEDIUM_SIZE = 8
    BIGGER_SIZE = 9

    plt.rc('font', size=SMALL_SIZE)
    plt.rc('axes', titlesize=SMALL_SIZE)
    plt.rc('axes', labelsize=SMALL_SIZE)
    plt.rc('xtick', labelsize=SMALL_SIZE)
    plt.rc('ytick', labelsize=SMALL_SIZE)
    plt.rc('legend', fontsize=SMALL_SIZE)
    plt.rc('figure', titlesize=MEDIUM_SIZE)
elif 'xfelbkr' in gethostname():
    SMALL_SIZE = 4
    MEDIUM_SIZE = 5
    BIGGER_SIZE = 6
    LINEWIDTH = 1
    MARKERSIZE = 2
    LEFT = 0.15
    RIGHT = 0.85
    BOTTOM = 0.75
    HSPACE = 0.75

    plt.rc('font', size=SMALL_SIZE)
    plt.rc('axes', titlesize=SMALL_SIZE)
    plt.rc('axes', labelsize=SMALL_SIZE)
    plt.rc('xtick', labelsize=SMALL_SIZE)
    plt.rc('ytick', labelsize=SMALL_SIZE)
    plt.rc('legend', fontsize=SMALL_SIZE)
    plt.rc('figure', titlesize=MEDIUM_SIZE)
    plt.rc('lines', linewidth=LINEWIDTH, markersize=MARKERSIZE)

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
        rec_point = result_dict['input']['correctors'][1]
        plot_handles = performance_figure(rec_point)
    fig, (sp_ellipse, sp_ellipse_norm, sp_A) = plot_handles

    emit_geo = result_dict['input']['emit_geo']
    beta = result_dict['input']['beta']
    alpha = result_dict['input']['alpha']
    A_range = result_dict['data']['A']
    phi_range = result_dict['data']['phi']

    def rescale_xxp(x, xp):
        x_out = x/np.sqrt(emit_geo*beta)
        xp_out = (x*alpha + xp*beta)/np.sqrt(emit_geo*beta)
        return x_out, xp_out


    x, xp = phase_space_ellipse(emit_geo, beta, alpha, 1000)
    xs, xps = rescale_xxp(x, xp)

    sp_ellipse.plot(x, xp, color='black')
    sp_ellipse_norm.plot(xs, xps, color='black')


    for n_A, A in enumerate(A_range):
        pulse_ene_mean = result_dict['data']['pulse_ene_mean'][n_A]
        notnan = ~np.isnan(pulse_ene_mean)
        if np.any(notnan):
            sp_A.plot(phi_range/np.pi*180, pulse_ene_mean*1e3, label=A, marker='.')

        x = result_dict['data']['delta_xxp'][n_A,:,0]
        xp = result_dict['data']['delta_xxp'][n_A,:,1]
        xs, xps = rescale_xxp(x, xp)

        for sp, xx, xxpp in [
                (sp_ellipse, x, xp),
                (sp_ellipse_norm, xs, xps),
                ]:
            sp.plot(xx, xxpp, ls='--', label=A)
            if np.any(notnan):
                sp.scatter(xx[notnan], xxpp[notnan], c=pulse_ene_mean[notnan])

    all_pulse_ene = result_dict['data']['pulse_ene_mean'].ravel()
    notnan = ~np.isnan(all_pulse_ene)
    if np.any(notnan):
        all_pulse_ene = all_pulse_ene[notnan]
        norm = Normalize(vmin=all_pulse_ene.min()*1e3, vmax=all_pulse_ene.max()*1e3)
        mappable = cm.ScalarMappable(norm=norm)
        fig.colorbar(mappable, ax=[sp_ellipse, sp_ellipse_norm])

    sp_A.axhline(result_dict['data']['init_pulse_energy']*1e3, color='black', ls='--', label=0)

    lim = sp_A.get_xlim()[0]
    sp_A.set_xlim(lim, -lim)

    for _sp in sp_ellipse, sp_A:
        _sp.legend(title='A')

def performance_figure(rec_point, figsize=[10, 12]):
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(left=LEFT, right=RIGHT, hspace=HSPACE, wspace=WSPACE, top=TOP, bottom=BOTTOM)
    sp_ellipse = plt.subplot(2, 2, 1)
    sp_ellipse.set_title('Phase space at %s' % rec_point)
    sp_ellipse.set_xlabel('$\Delta x$ ($\mu$m)')
    sp_ellipse.set_ylabel('$\Delta x\'$ ($\mu$rad)')

    sp_ellipse_norm = plt.subplot(2, 2, 2)
    sp_ellipse_norm.set_title('Norm. phase space at %s' % rec_point)
    sp_ellipse_norm.set_xlabel(r'$\Delta x/\sqrt{\epsilon\beta}$')
    sp_ellipse_norm.set_ylabel(r'$(\alpha\Delta x +\Delta x\'\beta)/\sqrt{\epsilon\beta}$')

    sp_A = plt.subplot(2, 2, 3)
    sp_A.set_title('Recorded pulse energies')
    sp_A.set_xlabel('Phase (deg)')
    sp_A.set_ylabel('Pulse energy (mJ)')

    return fig, (sp_ellipse, sp_ellipse_norm, sp_A)

