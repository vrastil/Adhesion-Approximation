"""
'plot.py' modules serves for plotting results
"""
import math
import matplotlib
matplotlib.use('Agg')
from matplotlib import animation
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from matplotlib.colors import SymLogNorm
from matplotlib.ticker import FormatStrFormatter
import numpy as np
from itertools import izip

from . import power

matplotlib.rcParams['legend.numpoints'] = 1
matplotlib.rcParams.update({'font.size': 15})

def iter_data(zs, iterables, a_end=None, a_slice=1.5, skip_init=True, get_a=False):
    """ Generator: iterate through data in list 'iterables'
    yield list of values when a_i > a_slice*a_i-1 and a_i < a_slice*a_end
    stops when a_i > a_end, a_end is the last value in zs, if not specified
    return string representation of z; 'z = ' + str(zs[i]); or 'init'
    """
    if a_end is None:
        a_end = 1./(zs[-1]+1)
    a_ = 0
    my_it = [iter(x) for x in iterables]
    for z in zs:
        values = [next(x) for x in my_it]
        if z != 'init':
            a = 1./(1.+z)
            if ((a < a_slice * a_) or (a_slice * a > a_end)) and a != a_end:
                continue
            elif a > a_end:
                raise StopIteration()
            a_ = a
            lab = 'z = ' + str(z)
        elif skip_init:
            continue
        else:
            a = 0
            lab = 'init'
        if get_a:
            yield [lab] + values + [a]
        else:
            yield [lab] + values

def close_fig(filename, fig, save=True, show=False):
    """save and/or show figure, close figure"""
    if save:
        fig.savefig(filename)
    if show:
        plt.show()
    fig.clf()
    plt.close(fig)

def add_nyquist_info(ax, a_sim_info):
    """plot lines corresponding to particle, potential and analys nyquist wavelengtsh"""
    ls = iter([':', '-.', '--'])
    val_lab = {}
    for key, val in a_sim_info.k_nyquist.iteritems():
        if val in val_lab:
            val_lab[val] += ",\n" + " " * 8 + key
        else:
            val_lab[val] = r"$k_{Nq}$ (" + key
    for val, lab in val_lab.iteritems():
        ax.axvline(val, ls=next(ls), c='k', label=lab + r")")

def legend_manipulation(ax, figtext):
    #handles, labels = zip(*sorted(zip(*ax.get_legend_handles_labels()), key=lambda x: x[1], reverse=True))
    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles, labels, loc='upper left',
               bbox_to_anchor=(1.0, 1.0), fontsize=14)
    plt.draw()
    plt.figtext(0.5, 0.95, figtext,
                bbox={'facecolor': 'white', 'alpha': 0.2}, size=14, ha='center', va='top')
    plt.subplots_adjust(left=0.1, right=0.84, bottom=0.1, top=0.89)    

def get_a_init_from_zs(zs):
    """ from list of redshifts returns initial scale factor, i.e. value after 'init' """
    for z in zs:
        if z != 'init':
            return 1/(1.+z)

def plot_pwr_spec(data, zs, a_sim_info, Pk_list_extrap, err=False,
                  out_dir='auto', pk_type='dens', save=True, show=False):
    """" Plot power spectrum -- points and extrapolated values,
    show 'true' linear Pk at the initial and final redshift """
    if out_dir == 'auto':
        out_dir = a_sim_info.res_dir
    if pk_type == "dens":
        out_file = 'pwr_spec.png'
        suptitle = "Power spectrum"
    elif pk_type == "vel":
        out_file = 'vel_pwr_spec.png'
        suptitle = r"Power spectrum $(\nabla\cdot u)$"
    fig = plt.figure(figsize=(15, 11))
    ax = plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')

    for lab, Pkk, Pk_ext in iter_data(zs, [data, Pk_list_extrap]):
        k, P_k = Pkk[0], Pkk[1]
        ax.plot(k, P_k, 'o', ms=3, label=lab)
        # show 1 standard deviation
        if err:
            P_k_std = Pkk[2]
            ax.fill_between(k, P_k - P_k_std, P_k + P_k_std,
                             facecolor='darkgrey', alpha=0.5)
        k = np.geomspace(k[0]/5,k[-1]) # extra half a decade for lin-/nl-/extrpolated-pk
        ax.plot(k, [Pk_ext(k_) for k_ in k], 'k--')

    add_nyquist_info(ax, a_sim_info)

    # plot non/linear power spectra
    a_0 = 1./(1.+zs[-1])
    a_i = get_a_init_from_zs(zs)
    P_i = power.lin_pow_spec(a_i, k, a_sim_info.sim.cosmo)
    P_0 = power.lin_pow_spec(a_0, k, a_sim_info.sim.cosmo)
    if pk_type == "dens":
        P_0_nl = power.non_lin_pow_spec(a_0, k, a_sim_info.sim.cosmo)
        ax.plot(k, P_0_nl, '-')
    elif pk_type == "vel":
        P_i *= power.growth_change(a_i, a_sim_info.sim.cosmo)**2
        P_0 *= power.growth_change(a_0, a_sim_info.sim.cosmo)**2
    ax.plot(k, P_0, '-')
    ax.plot(k, P_i, '-')
    
    fig.suptitle(suptitle, y=0.99, size=20)
    ax.set_xlabel(r"$k [h/$Mpc$]$", fontsize=15)
    ax.set_ylabel(r"$P(k) [$Mpc$/h)^3$]", fontsize=15)

    # LEGEND manipulation
    legend_manipulation(ax, a_sim_info.info_tr())

    # close & save figure
    close_fig(out_dir + out_file, fig, save=save, show=show)

def plot_chi_pwr_spec(data_list_chi, zs_chi, a_sim_info, err=False, out_dir='auto', save=True, show=False):
    if out_dir == 'auto':
        out_dir = a_sim_info.res_dir
    suptitle = "Chameleon power spectrum"
    out_file = "pwr_spec_chi.png"

    fig = plt.figure(figsize=(15, 11))
    ax = plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')

    for lab, Pkk, a in iter_data(zs_chi, [data_list_chi], get_a=True):
        k, P_k = Pkk[0], Pkk[1]
        chi_bulk_a_n = power.chi_bulk_a_n(a, a_sim_info.chi_opt)
        P_k /= pow(chi_bulk_a_n, 2)
        lines = ax.plot(k, P_k, 'o', ms=3, label=lab)
        color = lines[0].get_color()
        P_a = power.chi_lin_pow_spec(a, k, a_sim_info.sim.cosmo, a_sim_info.chi_opt)
        ax.plot(k, P_a, '-', color=color)
        if err:
            P_k_std = Pkk[2] / pow(chi_bulk_a_n, 2)
            ax.fill_between(k, P_k - P_k_std, P_k + P_k_std,
                             facecolor='darkgrey', alpha=0.5)

    add_nyquist_info(ax, a_sim_info)
    
    fig.suptitle(suptitle, y=0.99, size=20)
    ax.set_xlabel(r"$k [h/$Mpc$]$", fontsize=15)
    ax.set_ylabel(r"$P(k) [$Mpc$/h)^3$]", fontsize=15)

    # LEGEND manipulation
    legend_manipulation(ax, a_sim_info.info_tr())

    # close & save figure
    close_fig(out_dir + out_file, fig, save=save, show=show)

def plot_corr_func_single(corr_data, lab, a_sim_info, corr_data_lin=None, corr_data_nl=None, out_dir='auto', save=True, show=False):
    if out_dir == 'auto':
        out_dir = a_sim_info.res_dir
    z_out = lab if lab == 'init' else 'z' + lab[4:]
    # *****************
    # first plot, xi(r)
    # *****************
    fig = plt.figure(figsize=(15, 11))
    
    # load all data, plot xi(r)
    r, xi = corr_data
    plt.plot(r, xi, 'o', ms=3, label=lab)
    if corr_data_lin is not None:
        r_lin, xi_lin = corr_data_lin
        plt.plot(r_lin, xi_lin, '-', label=r"$\Lambda$CDM (lin)")
    if corr_data_nl is not None:
        r_nl, xi_nl = corr_data_nl
        plt.plot(r_nl, xi_nl, '-', label=r"$\Lambda$CDM (nl)")
    # adjust figure, labels
    fig.suptitle("Correlation function", y=0.99, size=20)
    plt.xlabel(r"$r [$Mpc$/h]$", fontsize=15)
    plt.ylabel(r"$\xi(r)$", fontsize=15)
    plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0), fontsize=14)
    plt.yscale("log")
    plt.xscale("log")
    plt.draw()
    plt.figtext(0.5, 0.95, a_sim_info.info_tr(),
                bbox={'facecolor': 'white', 'alpha': 0.2}, size=14, ha='center', va='top')
    plt.subplots_adjust(left=0.1, right=0.84, bottom=0.1, top=0.89)
    # save & show (in jupyter)
    close_fig(out_dir + 'corr_func_%s.png' % z_out, fig, save=save, show=show)

    # **********************
    # second plot, r*r*xi(r)
    # **********************
    fig = plt.figure(figsize=(15, 11))
    # plot r*r*xi(r)
    plt.plot(r, r * r * xi, 'o', ms=3, label=lab)
    if corr_data_lin is not None:
        plt.plot(r_lin, r_lin * r_lin * xi_lin,
                 '-', label=r"$\Lambda$CDM (lin)")
    if corr_data_nl is not None:
        plt.plot(r_nl, r_nl * r_nl * xi_nl,
                 '-', label=r"$\Lambda$CDM (nl)")
    # adjust figure, labels
    fig.suptitle("Correlation function", y=0.99, size=20)
    plt.xlabel(r"$r [$Mpc$/h]$", fontsize=15)
    plt.ylabel(r"$r^2\xi(r)$", fontsize=15)
    plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0), fontsize=14)
    plt.yscale("linear")
    plt.draw()
    plt.figtext(0.5, 0.95, a_sim_info.info_tr(),
                bbox={'facecolor': 'white', 'alpha': 0.2}, size=14, ha='center', va='top')
    plt.subplots_adjust(left=0.1, right=0.84, bottom=0.1, top=0.89)
    # save & show (in jupyter)
    close_fig(out_dir + 'corr_r2_func_%s.png' % z_out, fig, save=save, show=show)

# correlation function stacked data, linear and emu corr. func in files
def plot_corr_func(corr_data_all, zs, a_sim_info, out_dir='auto', save=True, show=False):
    for lab, corr_par, corr_lin, corr_nl in iter_data(zs, [corr_data_all['par'],
                                                      corr_data_all['lin'], corr_data_all['nl']]):
        plot_corr_func_single(
            corr_par, lab, a_sim_info, corr_lin, corr_nl, out_dir, save, show)


def plot_pwr_spec_diff_from_data(data_list, zs, a_sim_info, out_dir='auto', pk_type='dens', ext_title='par', save=True, show=False):
    if out_dir == 'auto':
        out_dir = a_sim_info.res_dir
    if pk_type == "dens":
        out_file = 'pwr_spec_diff'
        suptitle = "Power spectrum difference"
    elif pk_type == "vel":
        out_file = 'vel_pwr_spec_diff'
        suptitle = r"Power spectrum difference $(\nabla\cdot u)$"
    elif pk_type == 'chi':
        out_file = 'pwr_spec_diff_chi'
        suptitle = "Chameleon power spectrum difference"
        # transform chameleon power spectrum to suppression
        for z, data in izip(zs, data_list):
            a, k, Pk = 1/(1.+z), data[0], data[1]
            data[1] = power.chi_trans_to_supp(a, k, Pk, a_sim_info.sim.cosmo, a_sim_info.chi_opt)
            if len(data) == 3:
                data[2] = power.chi_trans_to_supp(a, k, data[2], a_sim_info.sim.cosmo, a_sim_info.chi_opt)
        # transform supp (ref: lin) to supp (ref: init)
        power.chi_trans_to_init(data_list)
        ext_title = 'init'
            
    out_file += '_%s.png' % ext_title
    suptitle += ' (ref: %s)' % ext_title

    fig = plt.figure(figsize=(15, 11))
    ax = plt.gca()
    plt.xscale('log')

    ymin = ymax = 0
    # SMALL / MEDIUM / LARGE SCALE VALUES
    # half of nyquist wavelength, 7 subintervals
    k = data_list[-1][0]
    idx = (np.abs(k - 0.5 * a_sim_info.k_nyquist["particle"])).argmin() / 7
    cmap = cm.get_cmap('gnuplot')
    ax.axvspan(k[0 * idx], k[1 * idx], alpha=0.2, color=cmap(0.1))
    ax.axvspan(k[3 * idx], k[4 * idx], alpha=0.3, color=cmap(0.5))
    ax.axvspan(k[6 * idx], k[7 * idx], alpha=0.4, color=cmap(0.9))

    for lab, data, a in iter_data(zs, [data_list], get_a=True):
        k, P_k = data[0], data[1]
        P_k_std = data[2] if len(data) == 3 else None
        if P_k_std is None:
            plt.plot(k, P_k, 'o', ms=3, label=lab)
            ymax = max(ymax, np.max(P_k[0:7 * idx]))
            ymin = min(ymin, np.min(P_k[0:7 * idx]))
        else:
            plt.errorbar(k, P_k, fmt='o', yerr=P_k_std, ms=3, label=lab)
            ymax = max(ymax, np.max(P_k[0:7 * idx] + P_k_std[0:7 * idx]))
            ymin = min(ymin, np.min(P_k[0:7 * idx] - P_k_std[0:7 * idx]))
            
    add_nyquist_info(ax, a_sim_info)

    if pk_type != 'chi' and ymax > 1:
        ymax = 1
    ymax = math.ceil(ymax / 0.1) * 0.1
    ymin = math.floor(ymin / 0.1) * 0.1
    if ymax == ymin:
        ymax += 0.1
        ymin -= 0.1
    plt.ylim(ymin=ymin, ymax=ymax)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.grid(True)

    fig.suptitle(suptitle, y=0.99, size=20)
    plt.xlabel(r"$k [h/$Mpc$]$", fontsize=15)
    plt.ylabel(r"$\frac{P(k)-P_{lin}(k)}{P_{lin}(k)}$", fontsize=25)
    plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0), fontsize=14)
    plt.draw()
    plt.figtext(0.5, 0.95, a_sim_info.info_tr(),
                bbox={'facecolor': 'white', 'alpha': 0.2}, size=14, ha='center', va='top')
    plt.subplots_adjust(left=0.1, right=0.84, bottom=0.1, top=0.89)

    close_fig(out_dir + out_file, fig, save=save, show=show)


def plot_pwr_spec_diff_map_from_data(data_array, zs, a_sim_info, out_dir='auto', pk_type='dens', ext_title='', save=True, show=False):
    if out_dir == 'auto':
        out_dir = a_sim_info.res_dir
    if pk_type == "dens":
        out_file = 'pwr_spec_diff'
        suptitle = "Power spectrum difference"
    elif pk_type == "vel":
        out_file = 'vel_pwr_spec_diff'
        suptitle = r"Power spectrum difference $(\nabla\cdot u)$"
    elif pk_type == 'chi':
        out_file = 'pwr_spec_diff_chi'
        suptitle = "Chameleon power spectrum difference"
        # transform chameleon power spectrum to suppression
        for z, data in izip(zs, data_array):
            a, k, Pk = 1/(1.+z), data[0], data[1]
            data[1] = -1 + power.chi_trans_to_supp(a, k, Pk, a_sim_info.sim.cosmo, a_sim_info.chi_opt)
        # transform supp (ref: lin) to supp (ref: init)
        power.chi_trans_to_init(data_array)
        ext_title = 'init'

    out_file += '_%s_map.png' % ext_title
    suptitle += ' (ref: %s)' % ext_title

    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(1, 15, wspace=0.5)
    ax = plt.subplot(gs[0, : -1])
    cbar_ax = plt.subplot(gs[0, -1])

    ax.set_xscale('log')
    a = [1 / (1 + z) for z in zs]
    # hack around pcolormesh plotting edges
    if len(a) == 1:
        da = 2*a[0]
    else:
        da = (a[-1] - a[0]) / (len(a) - 1)
    a = np.array([a[0]-da/2] + [1 / (1 + z) + da/2 for z in zs])
    k = data_array[0][0]
    supp = data_array[:, 1, :] # extract Pk, shape = (zs, k)
    vmin = -1
    vmax = 1
    if pk_type != 'chi':
        linthresh = 0.2
        linscale = 1.0
    else:
        linthresh = 0.5
        linscale = 0.2

    if vmin < 0:
        ticks = [vmin, -linthresh, 0, linthresh, vmax]
    else:
        ticks = [vmin, linthresh, vmax]
    labels = [str(x) for x in ticks]
    labels[-1] = '> %i' % ticks[-1]

    im = ax.pcolormesh(k, a, supp, cmap='seismic', norm=SymLogNorm(linthresh=linthresh, linscale=linscale,
                                   vmin=vmin, vmax=vmax))
    cbar = fig.colorbar(im, cax=cbar_ax, ticks=ticks)
    cbar.ax.set_yticklabels(labels)

    if a_sim_info.k_nyquist is not None:
        ls = [':', '-.', '--']
        ls *= (len(a_sim_info.k_nyquist) - 1) / 3 + 1
        ls = iter(ls)
        val_set = set(a_sim_info.k_nyquist.itervalues())
        for val in val_set:
            ax.axvline(val, ls=next(ls), c='k')

    fig.suptitle(suptitle, y=0.99, size=20)
    ax.set_xlabel(r"$k [h/$Mpc$]$", fontsize=15)
    ax.set_ylabel(r"$a(t)$", fontsize=15)
    plt.draw()
    plt.figtext(0.5, 0.95, a_sim_info.info_tr(),
                bbox={'facecolor': 'white', 'alpha': 0.2}, size=14, ha='center', va='top')
    plt.subplots_adjust(left=0.1, right=0.84, bottom=0.1, top=0.89)

    close_fig(out_dir + out_file, fig, save=save, show=show)

def plot_supp(sim_infos, out_dir, suptitle='', save=True, show=False, scale='', show_k_lms=False, res=None):
    fig = plt.figure(figsize=(15, 11))
    cmap = plt.get_cmap('gist_ncar')
    colors = [cmap(i) for i in np.linspace(0, 1, len(sim_infos) + 1)]
    for i, a_sim_info in enumerate(sim_infos):
        a = a_sim_info.a
        if scale == 'large':
            supp = a_sim_info.supp[0][0]
            if res is not None:
                supp -= np.array(res.supp[0][0])
        elif scale == 'medium':
            supp = a_sim_info.supp[0][1]
            if res is not None:
                supp -= np.array(res.supp[0][1])
        elif scale == 'small':
            supp = a_sim_info.supp[0][2]
            if res is not None:
                supp -= np.array(res.supp[0][2])
        else:
            print "WARNING! Unknown scale ='%s'. Skipping." % scale
            return None
        plt.plot(a, supp, '-o', ms=3,
                 color=colors[i], label=a_sim_info.info_supp())
        del a, supp

    if show_k_lms:
        if scale == 'large':
            suptitle += '<%.2f, %.2f> h/Mpc' % (
                a_sim_info.supp[1][0][0], a_sim_info.supp[1][0][1])
        elif scale == 'medium':
            suptitle += '<%.2f, %.2f> h/Mpc' % (
                a_sim_info.supp[1][1][0], a_sim_info.supp[1][1][1])
        elif scale == 'small':
            suptitle += '<%.2f, %.2f> h/Mpc' % (
                a_sim_info.supp[1][2][0], a_sim_info.supp[1][2][1])

    #plt.ylim(ymin=-1, ymax=0)
    fig.suptitle("Power spectrum suppression" + suptitle, y=0.99, size=20)
    plt.xlabel(r"$a(t)$", fontsize=15)
    ylabel = r"$\langle{\frac{P(k)-P_{lin}(k)}{P_{lin}(k)}}\rangle$"
    if res is not None:
        ylabel += r', residual from $\nu=%.1f$' % res.nu
    plt.ylabel(ylabel, fontsize=25)
    plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0), fontsize=14)
    plt.subplots_adjust(left=0.1, right=0.7, bottom=0.1, top=0.89)
    close_fig(out_dir + 'supp.png', fig, save=save, show=show)


def plot_dens_histo(data_list, zs, a_sim_info, out_dir='auto', fix_N=1, fix_rho=0.0, save=True, show=False):
    if out_dir == 'auto':
        out_dir = a_sim_info.res_dir
    num_sub_x = 3
    num_sub_y = (len(zs) + num_sub_x - 1) / num_sub_x
    fig = plt.figure(figsize=(num_sub_x * 5, num_sub_y * 5.5))

    gs = gridspec.GridSpec(num_sub_y, num_sub_x, wspace=0.2,
                           hspace=0.3, left=0.1, right=0.84, bottom=0.1, top=0.89)

    for lab, data, gs_cur in iter_data(zs, [data_list, gs]):
        rho, count = data
        count *= fix_N
        rho += fix_rho
        xmin = -1
        xmax = rho[np.nonzero(count)[0][-1]] + 1
        ax = plt.subplot(gs_cur)
        ax.set_xlim(xmin=xmin, xmax=xmax)
        ax.hist(rho, bins=20, weights=count, facecolor='green',
                edgecolor='black', linewidth=0.8, normed=True)
        ax.set_yscale('log', nonposy='clip')
        ax.set_title(lab)

    fig.suptitle("Overdensity distribution", y=0.99, size=20)

    plt.figtext(0.5, 0.95, a_sim_info.info_tr(),
                bbox={'facecolor': 'white', 'alpha': 0.2}, size=14, ha='center', va='top')
    close_fig(out_dir + 'dens_histo.png', fig, save=save, show=show)


def plot_par_last_slice(files, files_t, zs, a_sim_info, out_dir='auto', save=True, show=False):
    if out_dir == 'auto':
        out_dir = a_sim_info.res_dir
    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(1, 1)
    ax = plt.subplot(gs[0])
    data = np.loadtxt(files[0])
    x, y = data[:, 0], data[:, 1]
    ax.set_xlim(0, np.max(x))
    ax.set_ylim(0, np.max(y))
    data = np.loadtxt(files[-1])
    x, y = data[:, 0], data[:, 1]

    num_track = len(np.loadtxt(files_t[0]))
    data = np.loadtxt(files_t[-1])
    x_t, y_t = data[:, 0], data[:, 1]
    num_steps = len(x_t) / num_track
    x_t = [x_t[i:i + num_steps] for i in xrange(0, len(x_t), num_steps)]
    y_t = [y_t[i:i + num_steps] for i in xrange(0, len(y_t), num_steps)]

    ax.plot(x, y, 'ob', ms=1)
    for i in xrange(num_track):
        ax.plot(x_t[i], y_t[i], '--or', ms=4, lw=1.5,
                markevery=(num_steps - 1, num_steps))

    plt.figtext(0.5, 0.94, a_sim_info.info_tr(),
                bbox={'facecolor': 'white', 'alpha': 0.2}, size=14, ha='center', va='top')
    ax.set_xlabel(r"$x [$Mpc$/h]$", fontsize=13)
    ax.set_ylabel(r"$z [$Mpc$/h]$", fontsize=13)
    fig.suptitle("Slice through simulation box (particles), z = %.2f" %
                 zs[-1], y=0.99, size=20)
    close_fig(out_dir + 'par_evol_last.png', fig, save=save, show=show)


def plot_par_evol(files, files_t, zs, a_sim_info, out_dir='auto', save=True):
    if out_dir == 'auto':
        out_dir = a_sim_info.res_dir
    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(1, 1)
    ax = plt.subplot(gs[0])
    data = np.loadtxt(files[0])
    x, y = data[:, 0], data[:, 1]
    data = np.loadtxt(files_t[0])
    ax.set_xlim(0, np.max(x))
    ax.set_ylim(0, np.max(y))

    num = len(zs)
    num_track = len(data)

    line, = ax.plot([], [], 'ob', ms=1, animated=True)
    lines_t = []
    for i in xrange(num_track):
        lines_t.append(ax.plot([], [], '--or', ms=4, lw=1.5, animated=True)[0])

    plt.figtext(0.5, 0.94, a_sim_info.info_tr(),
                bbox={'facecolor': 'white', 'alpha': 0.2}, size=14, ha='center', va='top')
    ax.set_xlabel(r"$x [$Mpc$/h]$", fontsize=13)
    ax.set_ylabel(r"$z [$Mpc$/h]$", fontsize=13)
    del x, y, data

    def animate(j):
        if j < num:
            i = j
        else:
            i = 2 * num - j - 1
        fig.suptitle("Slice through simulation box (particles), z = %.2f" %
                     zs[i], y=0.99, size=20)
        data = np.loadtxt(files[i])
        x, y = data[:, 0], data[:, 1]
        data = np.loadtxt(files_t[i])
        x_t, y_t = data[:, 0], data[:, 1]
        num_steps = len(x_t) / num_track
        x_t = [x_t[i:i + num_steps] for i in xrange(0, len(x_t), num_steps)]
        y_t = [y_t[i:i + num_steps] for i in xrange(0, len(y_t), num_steps)]

        line.set_data(x, y)
        for i, line_t in enumerate(lines_t):
            line_t.set_data(x_t[i], y_t[i])
            line_t.set_markevery((num_steps - 1, num_steps))
        del x, y, x_t, y_t, data
        return [line] + lines_t

    ani = animation.FuncAnimation(
        fig, animate, frames=2 * num, interval=250, blit=True)
    if save:
        ani.save(out_dir + 'par_evol.gif', writer='imagemagick')
    del ani
    fig.clf()
    plt.close(fig)

def plot_dens_one_slice(rho, z, a_sim_info, out_dir='auto', save=True, show=False):
    if out_dir == 'auto':
        out_dir = a_sim_info.res_dir
    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(1, 15, wspace=0.5)
    ax = plt.subplot(gs[0, : -1])
    cbar_ax = plt.subplot(gs[0, -1])

    plt.figtext(0.5, 0.94, a_sim_info.info_tr(),
                bbox={'facecolor': 'white', 'alpha': 0.2}, size=14, ha='center', va='top')
    ax.set_xlabel(r"$x [$Mpc$/h]$", fontsize=13)
    ax.set_ylabel(r"$z [$Mpc$/h]$", fontsize=13)
    L = int(np.sqrt(rho.shape[0]))
    rho.shape = L, L
    im = ax.imshow(rho, interpolation='bicubic', cmap='gnuplot',
                   norm=SymLogNorm(linthresh=1.0, linscale=1,
                                   vmin=-1, vmax=100), aspect='auto',
                   extent=[0, a_sim_info.box_opt["box_size"], 0, a_sim_info.box_opt["box_size"]])
    fig.suptitle("Slice through simulation box (overdensity), z = %.2f" %
                 z, y=0.99, size=20)
    cbar = fig.colorbar(im, cax=cbar_ax, ticks=[-1, 0, 1, 10, 100])
    cbar.ax.set_yticklabels(['-1', '0', '1', '10', '> 100'])
    close_fig(out_dir + 'dens_z%.2f.png' % z, fig, save=save, show=show)

def plot_dens_two_slices(files, zs, a_sim_info, out_dir='auto', save=True, show=False):
    half = len(files) / 2
    rho, z = np.loadtxt(files[half])[:, 2], zs[half]
    plot_dens_one_slice(rho, z, a_sim_info,
                        out_dir=out_dir, save=save, show=show)
    rho, z = np.loadtxt(files[-1])[:, 2], zs[-1]
    plot_dens_one_slice(rho, z, a_sim_info,
                        out_dir=out_dir, save=save, show=show)


def plot_dens_evol(files, zs, a_sim_info, out_dir='auto', save=True):
    if out_dir == 'auto':
        out_dir = a_sim_info.res_dir

    num = len(zs)
    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(1, 15, wspace=0.5)
    ax = plt.subplot(gs[0, : -1])
    cbar_ax = plt.subplot(gs[0, -1])
    plt.figtext(0.5, 0.94, a_sim_info.info_tr(),
                bbox={'facecolor': 'white', 'alpha': 0.2}, size=14, ha='center', va='top')
    ax.set_xlabel(r"$x [$Mpc$/h]$", fontsize=13)
    ax.set_ylabel(r"$z [$Mpc$/h]$", fontsize=13)

    def animate(j):
        if j < num:
            i = j
        else:
            i = 2 * num - j - 1
        rho = np.loadtxt(files[i])[:, 2]
        L = int(np.sqrt(rho.shape[0]))
        rho.shape = L, L
        im = ax.imshow(rho, interpolation='bicubic', cmap='gnuplot', animated=True,
                       norm=SymLogNorm(
                           linthresh=1.0, linscale=1, vmin=-1, vmax=100), aspect='auto',
                       extent=[0, a_sim_info.box_opt["box_size"], 0, a_sim_info.box_opt["box_size"]])
        fig.suptitle(
            "Slice through simulation box (overdensity), z = %.2f" % zs[i], y=0.99, size=20)
        cbar = fig.colorbar(im, cax=cbar_ax, ticks=[-1, 0, 1, 10, 100])
        cbar.ax.set_yticklabels(['-1', '0', '1', '10', '> 100'])
        del rho
        return [im]

    ani = animation.FuncAnimation(
        fig, animate, frames=2 * num, interval=250, blit=True)
    if save:
        ani.save(out_dir + 'dens_evol.gif', writer='imagemagick')
    del ani
    fig.clf()
    plt.close(fig)


def plot_supp_lms(supp, a, a_sim_info, out_dir='auto', pk_type='dens', suptitle='', save=True, show=False):
    if out_dir == 'auto':
        out_dir = a_sim_info.res_dir
    if pk_type == "dens":
        out_file = 'supp.png'
        suptitle = "Power spectrum suppression"
    elif pk_type == "vel":
        out_file = 'supp_vel.png'
        suptitle = r"Power spectrum suppression $(\nabla\cdot u)$"
    elif pk_type == 'chi':
        out_file = 'supp_chi.png'
        suptitle = "Chameleon power spectrum suppression"

    fig = plt.figure(figsize=(15, 11))
    ax = plt.gca()
    cmap = cm.get_cmap('gnuplot')

    for i, scale in enumerate(['Large', 'Medium', 'Small']):
        ax.errorbar(a, supp[i][0], fmt='-o', yerr=supp[i][1], ms=3,
                    color=cmap(0.1+i*0.4), lw=4-i*1.5,
                    label='%s-scale:\n' r'$\langle%.2f,%.2f\rangle$' % (scale, supp[i][2][0], supp[i][2][0]))

    ymin, ymax = ax.get_ylim()
    if pk_type != 'chi' and ymax > 1:
        ymax = 1

    ymax = math.ceil(ymax / 0.1) * 0.1
    ymin = math.floor(ymin / 0.1) * 0.1
    plt.ylim(ymin=ymin, ymax=ymax)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.grid(True)
        
    fig.suptitle(suptitle, y=0.99, size=20)
    plt.xlabel(r"$a(t)$", fontsize=15)
    plt.ylabel(
        r"$\langle{\frac{P(k)-P_{lin}(k)}{P_{lin}(k)}}\rangle$", fontsize=25)
    plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0), fontsize=14)
 #   plt.subplots_adjust(left=0.1, right=0.7, top=0.9, bottom=0.1)
    plt.subplots_adjust(left=0.1, right=0.84, bottom=0.1, top=0.89)
    plt.draw()
    plt.figtext(0.5, 0.95, a_sim_info.info_tr(),
                bbox={'facecolor': 'white', 'alpha': 0.2}, size=14, ha='center', va='top')
    close_fig(out_dir + out_file, fig, save=save, show=show)


def plot_all_single_supp(res, out_dir='/home/vrastil/Documents/GIT/Adhesion-Approximation/output/supp_comparison/',
                         Nm=0, Np=0, L=0, nu=0, rs=0, app=''):
    subfiles = res.get_subfiles(Nm=Nm, Np=Np, L=L, nu=nu, rs=rs, app=app)
    for a_sim_info in subfiles:
        res.load_k_supp(a_sim_info)
        plot_supp_lms(a_sim_info.supp[0], a_sim_info.a, a_sim_info,
                      k_lms=a_sim_info.supp[1], show=True, save=True)
