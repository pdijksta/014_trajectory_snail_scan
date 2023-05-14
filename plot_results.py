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
    RIGHT = 0.9
    BOTTOM = 0.15
    WSPACE = 0.5
    HSPACE = 0.55

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

    sp_ellipse.plot(x*1e6, xp*1e6, color='black')
    sp_ellipse_norm.plot(xs, xps, color='black')


    for n_A, A in enumerate(A_range):
        pulse_ene_mean = result_dict['data']['pulse_ene_mean'][n_A]
        notnan = ~np.isnan(pulse_ene_mean)
        if np.any(notnan):
            sp_A.plot(phi_range/np.pi*180, pulse_ene_mean*1e3, label=A, marker='.')

        x = result_dict['data']['delta_xxp'][n_A,:,0]
        xp = result_dict['data']['delta_xxp'][n_A,:,1]
        xs, xps = rescale_xxp(x, xp)

        for sp, xx, xxpp, factor in [
                (sp_ellipse, x, xp, 1e6),
                (sp_ellipse_norm, xs, xps, 1),
                ]:
            sp.plot(xx*factor, xxpp*factor, ls='--', label=A)
            if np.any(notnan):
                sp.scatter(xx[notnan]*factor, xxpp[notnan]*factor, c=pulse_ene_mean[notnan]*1e3, zorder=100)

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
    sp_ellipse_norm.set_ylabel(r'$(\alpha\Delta x +\beta\Delta x\')/\sqrt{\epsilon\beta}$')

    sp_A = plt.subplot(2, 2, 3)
    sp_A.set_title('Recorded pulse energies')
    sp_A.set_xlabel('Phase (deg)')
    sp_A.set_ylabel('Pulse energy (mJ)')

    return fig, (sp_ellipse, sp_ellipse_norm, sp_A)

def plot_orbit(result_dict, plot_handles=None):
    if plot_handles is None:
        correctors = result_dict['input']['correctors']
        plane = result_dict['input']['plane']
        maxA  = result_dict['data']['A'].max()
        plot_handles = orbit_figure(correctors, maxA, plane)
    fig, (sp_angle0, sp_angle1, sp_orbit) = plot_handles

    init_angles = result_dict['input']['init_corrector_vals']
    A_range = result_dict['data']['A']
    phi_range = result_dict['data']['phi']
    sp_angles = [sp_angle0, sp_angle1]
    for A_index, A in enumerate(A_range):
        for n_corr in range(2):
            angles = result_dict['data']['delta_corr_angles'][A_index,:,n_corr] + init_angles[n_corr]
            pulse_energies = result_dict['data']['pulse_ene_mean'][A_index]
            notnan = ~np.isnan(pulse_energies)
            if np.any(notnan):
                sp_angles[n_corr].plot(phi_range[notnan]/np.pi*180, angles[notnan]*1e6, ls='--', label='%.1f' % A)
                sp_angles[n_corr].scatter(phi_range[notnan]/np.pi*180, angles[notnan]*1e6, c=pulse_energies[notnan])

    all_pulse_ene = result_dict['data']['pulse_ene_mean'].ravel()
    notnan = ~np.isnan(all_pulse_ene)
    if np.any(notnan):
        all_pulse_ene = all_pulse_ene[notnan]
        norm = Normalize(vmin=all_pulse_ene.min()*1e3, vmax=all_pulse_ene.max()*1e3)
        mappable = cm.ScalarMappable(norm=norm)
        fig.colorbar(mappable, ax=sp_angles)

    for _sp in sp_angles:
        lim = _sp.get_xlim()[0]
        _sp.set_xlim(lim, -lim)
        _sp.legend(title='A')

    bpmz = [int(x.split('.')[1]) for x in result_dict['data']['bpm_names']]
    for n_phi, phi in enumerate(phi_range):

        orbit_mean = result_dict['data']['orbit_mean'][-1,n_phi]
        orbit_std = result_dict['data']['orbit_std'][-1,n_phi]

        if np.any(~np.isnan(orbit_mean)):
            sp_orbit.errorbar(bpmz, orbit_mean*1e3, yerr=orbit_std*1e3, label='%i' % (phi/np.pi*180))
    sp_orbit.legend(title='Phi')


def orbit_figure(correctors, maxA, plane, figsize=[10, 12]):
    fig = plt.figure(figsize=figsize)
    fig.subplots_adjust(left=LEFT, right=RIGHT, hspace=HSPACE, wspace=WSPACE, top=TOP, bottom=BOTTOM)
    spec = fig.add_gridspec(ncols=2, nrows=2)

    sp_angles = []
    for n_corr, corrector in enumerate(correctors):
        sp_angle = fig.add_subplot(spec[0,n_corr])
        sp_angle.set_title(corrector)
        sp_angle.set_xlabel('Phi (deg)')
        sp_angle.set_ylabel('Angle ($\mu$rad)')
        sp_angles.append(sp_angle)

    sp_orbit = fig.add_subplot(spec[1,:])
    sp_orbit.set_title('Orbit for A=%.1f' % maxA)
    sp_orbit.set_xlabel('s (m)')
    sp_orbit.set_ylabel('%s (mm)' % plane)

    return fig, (sp_angles[0], sp_angles[1], sp_orbit)

