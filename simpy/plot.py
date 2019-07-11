"""
'plot.py' modules serves for plotting results
"""

from __future__ import print_function

import math
import matplotlib
matplotlib.use('Agg', warn=False)
from matplotlib import animation
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from matplotlib.colors import SymLogNorm
from matplotlib.ticker import FormatStrFormatter
from cycler import cycler

import numpy as np
try:
    from itertools import izip # python 2
except ImportError:
    izip = zip # python 3
from scipy.misc import derivative
#from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit

class MySpline(object):
    def __init__(self, x, y, func, p0=0):
        self.func = func
        self.popt, self.pcov = curve_fit(func, x, y, p0=p0)

    def __call__(self, x):
        return self.func(x, *self.popt)

from . import power
from . import utils as ut

# default values
matplotlib.rcParams['legend.numpoints'] = 1
matplotlib.rcParams['lines.linewidth'] = 3.0
matplotlib.rcParams['lines.markersize'] = 4.0
# matplotlib.rcParams['axes.prop_cycle'] = (
#     cycler('color', ['red', 'green', 'blue', 'orange', 'brown', 'violet', 'black']) +
#     cycler('ls', ['solid', 'dashed', 'dashdot', (0, (5, 10)), (0, (3, 5, 1, 5)), (0, (3, 1, 1, 1)), (0, (1, 1))])
# )
matplotlib.rcParams['axes.labelsize'] = 30
matplotlib.rcParams['xtick.labelsize'] = 20
matplotlib.rcParams['ytick.labelsize'] = 20
matplotlib.rcParams['legend.fontsize'] = 25
matplotlib.rcParams['font.size'] = 25

suptitle_size = 25
fig_size = (14, 9)
fig_size_map = (14, 14)
subplt_adj_sym = {'left' : 0.15, 'right' : 0.95, 'bottom' : 0.15, 'top' : 0.95}
report_dir = "/home/michal/Documents/GIT/FastSim/report/clanek/"

def iter_data(zs, iterables, a_end=None, a_slice=1.5, skip_init=True, get_a=False, only_last=False):
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
            if only_last and a != a_end:
                continue
            elif a > a_end:
                raise StopIteration()
            a_ = a
            lab = 'z = %.1f' % z
        elif skip_init or only_last:
            continue
        else:
            a = 0
            lab = 'init'
        if get_a:
            yield [lab] + values + [a]
        else:
            yield [lab] + values

def fig_suptitle(fig, suptitle="", y=0.99, size=suptitle_size):
    #fig.suptitle(suptitle, y=0.99, size=suptitle_size)
    pass

def close_fig(filename, fig, save=True, show=False, dpi=500, use_z_eff=False, format='eps', ext_legen=None):
    """save and/or show figure, close figure"""
    if use_z_eff:
        filename += '_z_eff'
    if save:
        if format == 'all':
            fig.savefig(filename + ".png", dpi=dpi, format="png")
            fig.savefig(filename + ".eps", dpi=dpi, format="eps")
        else:
            fig.savefig(filename + ".%s" % format, dpi=dpi, format=format)
        
        if ext_legen is not None:
            ext_legen.savefig(filename + "_legend.eps", bbox_inches='tight', dpi=dpi, format="eps")
    if show:
        plt.show()

    fig.clf()
    plt.close(fig)

    if ext_legen is not None:
        ext_legen.clf()
        plt.close(ext_legen)


def add_nyquist_info(ax, a_sim_info):
    """plot lines corresponding to particle, potential and analys nyquist wavelengtsh"""
    ls = iter([':', '-.', '--'])
    val_lab = {}
    for key, val in a_sim_info.k_nyquist.iteritems():
        if val in val_lab:
            val_lab[val] += ",\n" + " " * 8 + key
        else:
            val_lab[val] = r"$k_{Nq}$ (" + key
    for val, lab in val_lab.items():
        ax.axvline(val, ls=next(ls), c='k', label=lab + r")")

def trans_legend(handles, labels, ncol):
    n = len(labels)
    nrow = (n - 1) / ncol + 1
    order = []
    for i in range(n):
        row = i
        order.append()

def legend_manipulation(ax=None, figtext="", loc='upper left', bbox_to_anchor=(1.0, 1.0), ext_legen=False):
    ax = plt.gca() if ax is None else ax
    fig = plt.gcf()
    handles, labels = ax.get_legend_handles_labels()
    if ext_legen:
        if isinstance(ext_legen, dict):
            mlt_col = ext_legen['mlt_col']
            ncol = ext_legen['ncol']
            nrow = (len(labels) - 1) / ncol + 1
            ext_legen = plt.figure(figsize=(14, nrow*mlt_col))
            ax_leg = plt.gca()
            ax_leg.legend(handles, labels, loc='center', ncol=ncol, frameon=False, mode='expand')
            ax_leg.axes.get_xaxis().set_visible(False)
            ax_leg.axes.get_yaxis().set_visible(False)
            ext_legen.subplots_adjust(left=0.01, right=0.99, bottom=0.01, top=0.99)
            ext_legen.canvas.draw()
        else:
            ext_legen = False
    else:
        ax.legend(handles, labels, loc=loc,
                   bbox_to_anchor=bbox_to_anchor, fontsize=14)
    fig.canvas.draw()
    if figtext != "":
        plt.figtext(0.5, 0.95, figtext,
                bbox={'facecolor': 'white', 'alpha': 0.2}, size=14, ha='center', va='top', figure=fig)
    fig.subplots_adjust(left=0.1, right=0.84, bottom=0.1, top=0.89)
    return ext_legen


def adjust_extreme_values(ax, ymin, ymax):
    ymin_cur, ymax_cur = ax.get_ylim()
    if ymin_cur < ymin:
        ymin_cur = ymin
    if ymax_cur > ymax:
        ymax_cur = ymax
    ax.set_ylim(ymin_cur, ymax_cur)

def adjust_extreme_values_range(ax, x, y, x_max):
    y_valid = y[np.where(x < x_max)]
    y_min = np.min(y_valid)*0.8
    y_max = np.max(y_valid)*1.2
    ax.set_ylim(y_min, y_max)

def get_chi_label(si, single=False):
    label = r"$\Phi_{\rm scr}=%.1e,\ n=%.1f$" % (
                si.chi_opt["phi"], si.chi_opt["n"])
    if not single:
        label += r" (lin)" if si.chi_opt["linear"] else r" (nl)"

    return label

def get_chi_labels(sim_infos, single=False):
    return [get_chi_label(si, single=single) for si in sim_infos]

def plot_pwr_spec(data, zs, a_sim_info, Pk_list_extrap,
                  out_dir='auto', pk_type='dens', save=True, show=False, use_z_eff=False):
    """" Plot power spectrum -- points and extrapolated values,
    show 'true' linear Pk at the initial and final redshift """
    if out_dir == 'auto':
        out_dir = a_sim_info.res_dir
    if pk_type == "dens":
        out_file = 'pwr_spec'
        suptitle = "Power spectrum"
    elif pk_type == "vel":
        out_file = 'vel_pwr_spec'
        suptitle = r"Power spectrum $(\nabla\cdot u)$"
    fig = plt.figure(figsize=fig_size)
    ax = plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')

    for lab, Pkk, Pk_ext in iter_data(zs, [data, Pk_list_extrap]):
        k, P_k = Pkk[0], Pkk[1]
        ax.plot(k, P_k, 'o', label=lab)
        # show 1 standard deviation
        if len(Pkk) == 3:
            P_k_std = Pkk[2]
            ax.fill_between(k, P_k - P_k_std, P_k + P_k_std,
                             facecolor='darkgrey', alpha=0.5)
        k = np.geomspace(k[0]/5,k[-1]) # extra half a decade for lin-/nl-/extrpolated-pk
        ax.plot(k, [Pk_ext(k_) for k_ in k], 'k--')

    add_nyquist_info(ax, a_sim_info)

    # plot non/linear power spectra
    a_0 = 1./(1.+zs[-1])
    a_i = power.get_a_init_from_zs(zs)
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
    
    fig_suptitle(fig, suptitle)
    ax.set_xlabel(r"$k [h{\rm Mpc}^{-1}]$")
    ax.set_ylabel(r"$P(k) [(h^{-1}{\rm Mpc})^3]$")

    # LEGEND manipulation
    # legend_manipulation(ax, a_sim_info.info_tr())
    legend_manipulation(ax, "")

    # close & save figure
    close_fig(out_dir + out_file, fig, save=save, show=show, use_z_eff=use_z_eff)

def plot_pwr_spec_comparison(Pk_list_extrap, data, zs, labels, cosmo, out_dir='auto',
            save=True, show=False, use_z_eff=False, scale_to_lin=True, pk_type='dens', k_max=None, chi=False):
    """" Plot power spectrum -- points and extrapolated values,
    show 'true' linear Pk at the initial and final redshift """
    if out_dir == 'auto':
        out_dir = report_dir

    if pk_type == 'dens':
        out_file = 'pwr_spec'
    elif pk_type == 'vel':
        out_file = 'vel_pwr_spec'
    else:
        raise KeyError("Uknown pk_type: %s" % pk_type)
    suptitle = "Power spectrum"

    fig = plt.figure(figsize=fig_size)
    ax = plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')

    for _, Pkk, lab, Pk_extrap in iter_data(zs, [data, labels, Pk_list_extrap]):
        k, P_k, P_k_std = Pkk[0], Pkk[1], Pkk[2]

        if k_max is not None:
            idx = np.where(k < k_max)
            k = k[idx]
            P_k = P_k[idx]
            P_k_std = P_k_std[idx]

        if scale_to_lin:
            P_k_tmp = P_k / Pk_extrap.A_low
        else:
            P_k_tmp = P_k

        ax.plot(k, P_k_tmp, 'o', label=lab)   

        # show 1 standard deviation
        ax.fill_between(k, P_k_tmp - P_k_std, P_k_tmp + P_k_std,
                        facecolor='darkgrey', alpha=0.5)

    # plot non/linear power spectra
    k = np.geomspace(k[0],k[-1])
    a_0 = 1./(1.+zs[-1])
    P_0 = power.lin_pow_spec(a_0, k, cosmo)
    P_0_nl = power.non_lin_pow_spec(a_0, k, cosmo)

    if pk_type == 'vel' and not scale_to_lin:
        P_0 *= power.growth_change(a_0, cosmo)**2
        P_0_nl *= power.growth_change(a_0, cosmo)**2

    ax.plot(k, P_0, '-', label=r"$\Lambda$CDM (lin)")
    ax.plot(k, P_0_nl, '-',  label=r"$\Lambda$CDM (nl)")

    fig_suptitle(fig, suptitle, y=0.95)
    ax.set_xlabel(r"$k [h{\rm Mpc}^{-1}]$")
    ax.set_ylabel(r"$P(k) [(h^{-1}{\rm Mpc})^3]$")

    # LEGEND manipulation
    if chi:
        ext_legen = {'mlt_col' : 0.6, 'ncol' : 2}
    else:
        ext_legen = {'mlt_col' : 0.5, 'ncol' : 4}
    fig_leg = legend_manipulation(ax, "", ext_legen=ext_legen)
    fig.subplots_adjust(**subplt_adj_sym)

    # close & save figure
    close_fig(out_dir + out_file, fig, save=save, show=show, use_z_eff=use_z_eff, ext_legen=fig_leg)

def plot_pwr_spec_comparison_ratio_nl(Pk_list_extrap, data, zs, labels, cosmo, out_dir='auto',
        save=True, show=False, use_z_eff=False, scale_to_lin=True, pk_type='dens', k_max=None, chi=False):
    """" Plot power spectrum -- points and extrapolated values,
    show 'true' linear Pk at the initial and final redshift """
    if out_dir == 'auto':
        out_dir = report_dir
    if pk_type == 'dens':
        out_file = 'pwr_spec_ratio_nl'
    elif pk_type == 'vel':
        out_file = 'vel_pwr_spec_ratio_nl'
    else:
        raise KeyError("Uknown pk_type: %s" % pk_type)

    suptitle = "Power spectrum"

    fig = plt.figure(figsize=fig_size)
    ax = plt.gca()
    # ax.set_yscale('symlog', linthreshy=0.1, linscaley=1)
    # ax.set_yscale('log')
    ax.set_xscale('log')

    # get non-linear power spectra
    k = data[0][0]
    if k_max is not None:
        idx = np.where(k < k_max)
        k = k[idx]
    a_0 = 1./(1.+zs[-1])
    P_0_nl = power.non_lin_pow_spec(a_0, k, cosmo)

    for _, Pkk, lab, Pk_extrap in iter_data(zs, [data, labels, Pk_list_extrap]):
        k, P_k, P_k_std = Pkk[0], Pkk[1], Pkk[2]

        if k_max is not None:
            idx = np.where(k < k_max)
            k = k[idx]
            P_k = P_k[idx]
            P_k_std = P_k_std[idx]

        if scale_to_lin:
            P_k_tmp = P_k / Pk_extrap.A_low
        else:
            P_k_tmp = P_k

        ax.errorbar(k, P_k_tmp / P_0_nl, yerr=P_k_std/P_0_nl, fmt='o', label=lab)
        # ax.plot(k, P_k_tmp / P_0_nl - 1, 'o', label=lab)

    # plot linear power spectra
    k_ = np.geomspace(k[0],k[-1], num=200)
    a_0 = 1./(1.+zs[-1])
    P_0 = power.lin_pow_spec(a_0, k_, cosmo)
    P_0_nl_ = power.non_lin_pow_spec(a_0, k_, cosmo)

    ax.plot(k_, P_0 / P_0_nl_, '-', label=r"$\Lambda$CDM (lin)")
    
    fig_suptitle(fig, suptitle, y=0.95)
    ax.set_xlabel(r"$k [h{\rm Mpc}^{-1}]$")
    ax.set_ylabel(r"$P(k)/P_{\rm nl}(k)$")

    #ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.grid(True)

    # LEGEND manipulation
    legend_manipulation(ax, "", ext_legen=True)
    fig.subplots_adjust(**subplt_adj_sym)

    # close & save figure
    close_fig(out_dir + out_file, fig, save=save, show=show, use_z_eff=use_z_eff)


def plot_chi_pwr_spec(data_list_chi, zs_chi, a_sim_info, out_dir='auto', save=True, show=False, use_z_eff=False):
    if out_dir == 'auto':
        out_dir = a_sim_info.res_dir
    suptitle = "Chameleon power spectrum"
    out_file = "pwr_spec_chi"

    fig = plt.figure(figsize=fig_size)
    ax = plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')

    for lab, Pkk, a in iter_data(zs_chi, [data_list_chi], get_a=True):
        k, P_k = Pkk[0], Pkk[1]
        chi_bulk_a_n = power.chi_bulk_a_n(a, a_sim_info.chi_opt)
        P_k /= pow(chi_bulk_a_n, 2)
        lines = ax.plot(k, P_k, 'o', label=lab)
        color = lines[0].get_color()
        P_a = power.chi_lin_pow_spec(a, k, a_sim_info.sim.cosmo, a_sim_info.chi_opt)
        ax.plot(k, P_a, '-', color=color)
        if len(Pkk) == 3:
            P_k_std = Pkk[2] / pow(chi_bulk_a_n, 2)
            ax.fill_between(k, P_k - P_k_std, P_k + P_k_std,
                             facecolor='darkgrey', alpha=0.5)

    add_nyquist_info(ax, a_sim_info)
    fig_suptitle(fig, suptitle)
    ax.set_xlabel(r"$k [h{\rm Mpc}^{-1}]$")
    ax.set_ylabel(r"$P(k) [(h^{-1}{\rm Mpc})^3]$")

    # LEGEND manipulation
    # legend_manipulation(ax, a_sim_info.info_tr())
    legend_manipulation(ax, "")

    # close & save figure
    close_fig(out_dir + out_file, fig, save=save, show=show, use_z_eff=use_z_eff)


def plot_chi_fp_map(data_array, zs, chi_info, out_dir='auto', save=True, show=False, shading='flat',
                    max_nyquist=True, cut_low=False, vmin=1, vmax=3.0):
    #
    fig = plt.figure(figsize=fig_size_map)
    gs = gridspec.GridSpec(1, 60, wspace=0.5)
    ax = plt.subplot(gs[0, : -5])
    cbar_ax = plt.subplot(gs[0, -4 :])

    ax.set_xscale('log')
    a_z = [1 / (1 + z) for z in zs]
    
    # non-linear scale
    # k_nl = power.chi_psi_k_a(a_z, chi_info.sim.cosmo, chi_info.chi_opt)
    # chameleon mass
    m_chi = np.sqrt(power.chi_mass_sq(a_z, chi_info.sim.cosmo, chi_info.chi_opt))

    # hack around pcolormesh plotting edges
    if shading == 'flat':
        if len(a_z) == 1:
            da = 2*a_z[0]
        else:
            da = (a_z[-1] - a_z[0]) / (len(a_z) - 1)
        a = np.array([a_z[0]-da/2] + [1 / (1 + z) + da/2 for z in zs])
    
    k = data_array[0][0]
    supp = data_array[:, 1, :] # extract Pk, shape = (zs, k)

    linthresh = 1.05
    linscale = 10.0

    # cut
    idx_low = np.where(supp > 1.25)[0][-1] if cut_low else None
    idx_up = (np.abs(k - chi_info.k_nyquist["particle"])).argmin() if max_nyquist else None
    k = k[idx_low:idx_up]
    supp = supp[:, idx_low:idx_up]

    if not max_nyquist:
        add_nyquist_info(ax, chi_info)
            
    if vmin < 0:
        ticks = [vmin, -linthresh, 0, linthresh, vmax]
    else:
        ticks = [vmin, linthresh, vmax]
    labels = [str(x) for x in ticks]
    labels[-1] = '> %.1f' % ticks[-1]

    im = ax.pcolormesh(k, a, supp, cmap='Reds', shading=shading,
        norm=SymLogNorm(linthresh=linthresh, linscale=linscale, vmin=vmin, vmax=vmax))
    cbar = fig.colorbar(im, cax=cbar_ax, ticks=ticks)
    cbar.ax.set_yticklabels(labels)

    ax.set_xlabel(r"$k [h{\rm Mpc}^{-1}]$")
    ax.set_ylabel(r"$a(t)$")
    ax.tick_params(axis='both', which='major')
    plt.draw()

    # plot k_nl, keep ylim
    xmin, xmax = ax.get_xlim()
    ax.plot(m_chi, a_z, 'b-', lw=3)
    ax.set_xlim(xmin, xmax)

    plt.figtext(0.5, 0.95, "",
                bbox={'facecolor': 'white', 'alpha': 0.2}, size=14, ha='center', va='top')
    fig.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9)

    # close, save show
    if out_dir == 'auto':
        out_dir = report_dir
    bo = chi_info.box_opt
    out_file = 'chi_pwr_diff_map_fp_%im_%ip_%iM_%ib' % (bo["mesh_num"], bo["Ng"], bo["mesh_num_pwr"], bo["box_size"])

    if chi_info.chi_opt["linear"]:
        out_file += "_lin"
    else:
        out_file += "_nl"
    close_fig(out_dir + out_file, fig, save=save, show=show, format='png')

def plot_chi_fp_z(data_z, a_sim_info, labels, out_dir='auto', suptitle='auto', save=True, show=False, use_z_eff=False, max_nyquist=True):
    if out_dir == 'auto':
        out_dir = a_sim_info.res_dir
    bo = a_sim_info.box_opt
    out_file = 'chi_pwr_diff_fp_%im_%ip_%iM_%ib' % (bo["mesh_num"], bo["Ng"], bo["mesh_num_pwr"], bo["box_size"])
    # if suptitle == 'auto':
    #     suptitle = "Relative chameleon power spectrum"

    fig = plt.figure(figsize=fig_size)
    ax = plt.gca()
    plt.xscale('log')
    ymax = 1
    ymin = 0.95
    for data_chi, label in izip(data_z, labels): # each chi
        k = data_chi[0]
        Pk = data_chi[1]
        # std = data_chi[2]

        if max_nyquist:
            idx = (np.abs(k - a_sim_info.k_nyquist["particle"])).argmin()
            k = k[0:idx]
            Pk = Pk[0:idx]
            # std = std[0:idx]

        ymax = max(ymax, np.max(Pk))
        if "(lin)" in label:
            color = ax.get_lines()[-1].get_color()
            ls = ":"
        else:
            color = None
            ls = "-"

        # ax.errorbar(k, Pk, fmt='o', yerr=std, label=label)
        ax.plot(k, Pk, 'o', ls=ls, c=color, label=label)

    if not max_nyquist:
        add_nyquist_info(ax, a_sim_info)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.grid(True)

    ymax *= 1.1
    ax.set_ylim(ymin, ymax)

    fig_suptitle(fig, suptitle)
    plt.xlabel(r"$k [h{\rm Mpc}^{-1}]$")
    plt.ylabel(r"${P_\chi(k)}/{P_{\rm FP}(k)}$")
    #figtext = a_sim_info.info_tr().replace("FP: ", "")
    legend_manipulation(ax, figtext="", loc='upper left', bbox_to_anchor=(0.0,1.0))
    fig.subplots_adjust(**subplt_adj_sym)
    close_fig(out_dir + out_file, fig, save=save, show=show, use_z_eff=use_z_eff)


def plot_chi_fp_res_ax(ax, data_chi, si, ymax):
    k = data_chi[0]
    Pk = data_chi[1]
    # std = data_chi[2]
    ymax = max(ymax, np.max(Pk))
    label = r"$k_{nq} = %.1e$" % si.k_nyquist["potential"]
    if si.chi_opt["linear"]:
        label += " (lin)"
        color = ax.get_lines()[-1].get_color()
        ls = ":"
    else:
        label += " (nl)"
        color = None
        ls = "-"
    # ax.errorbar(k, Pk, fmt='o', yerr=std, label=label)
    ax.plot(k, Pk, 'o', ls=ls, label=label, c=color)
    return ymax

def plot_chi_fp_res(data_all, sim_infos, out_dir='auto', suptitle='auto', save=True, show=False):
    if out_dir == 'auto':
        out_dir = report_dir
    out_file = 'chi_resolution_eff'

    fig = plt.figure(figsize=fig_size)
    ax = plt.gca()
    ax.set_xscale('log')
    # ax.set_yscale('log')
    ax.set_yscale('symlog', linthreshy=2, linscaley=2)
    ymax = 1
    ymin = 0.95

    # go through all groups FP-CHI with the same resolution
    for data_chi_gr, si_gr in izip(data_all, sim_infos):
        # go through all chi data, i.e. linear and non-linear, plot first non-linear
        order = [1, 0] if si_gr[0].chi_opt["linear"] else [0, 1]

        for idx in order:
            ymax = plot_chi_fp_res_ax(ax, data_chi_gr[idx], si_gr[idx], ymax)

        # plot appropriate nyquist frequency
        color = ax.get_lines()[-1].get_color()
        ax.axvline(x=si_gr[0].k_nyquist["potential"], ls='--', c=color)

    #add_nyquist_info(ax, a_sim_info)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.grid(True)

    ymax *= 1.1
    ax.set_ylim(ymin, ymax)

    fig_suptitle(fig, suptitle)
    plt.xlabel(r"$k [h{\rm Mpc}^{-1}]$")
    plt.ylabel(r"${P_\chi(k)}/{P_{\rm FP}(k)}$")
    #figtext = a_sim_info.info_tr().replace("FP: ", "")
    legend_manipulation(ax, figtext="", loc='upper left', bbox_to_anchor=(0.0,1.0))
    fig.subplots_adjust(**subplt_adj_sym)
    close_fig(out_dir + out_file, fig, save=save, show=show)

def get_slope(k, P_k, dx=0.01,order=5):
    logk = np.log(k)
    logP_k = lambda logk : np.log(P_k(np.exp(logk)))
    return [derivative(logP_k, logk_, dx=dx, order=order) for logk_ in logk]


def plot_slope(data, zs, a_sim_info, Pk_list_extrap,
                  out_dir='auto', save=True, show=False, use_z_eff=False):
    """" Plot slope of power spectrum -- points and extrapolated values """
    if out_dir == 'auto':
        out_dir = a_sim_info.res_dir
    out_file = 'pwr_slope'
    suptitle = "Power spectrum slope"
    fig = plt.figure(figsize=fig_size)
    ax = plt.gca()
    ax.set_xscale('log')
    ax.set_ylim(-4,2)

    # get_slope = lambda k, P_k : [k_/P_k(k_)*derivative(P_k, k_, dx=k_/4) for k_ in k]

    for lab, Pkk, Pk_ext in iter_data(zs, [data, Pk_list_extrap], only_last=True):
        k, P_k = Pkk[0], Pkk[1]
        slope = np.diff(np.log(P_k))/np.diff(np.log(k))
        k_half = (k[1:] + k[:-1]) / 2.
        ax.plot(k_half, slope, 'o', label=lab)
        k = np.geomspace(k[0]/5,k[-1], num=400) # extra half a decade for lin-/nl-/extrpolated-pk
        slope = get_slope(k, Pk_ext, dx=0.2)
        ax.plot(k, slope, '--')

    add_nyquist_info(ax, a_sim_info)

    # plot non/linear power spectra
    a_0 = 1./(1.+zs[-1])
    P_0 = lambda x : power.lin_pow_spec(a_0, x, a_sim_info.sim.cosmo)
    P_0_nl = lambda x : power.non_lin_pow_spec(a_0, x, a_sim_info.sim.cosmo)
    slope = get_slope(k, P_0)
    ax.plot(k, slope, '-', label=r"$\Lambda$CDM (lin)")
    slope = get_slope(k, P_0_nl)
    ax.plot(k, slope, '-', label=r"$\Lambda$CDM (nl)")

    fig_suptitle(fig, suptitle)
    ax.set_xlabel(r"$k [h{\rm Mpc}^{-1}]$")
    ax.set_ylabel(r"d$\ln P(k)/$d$\ln k$]")

    # LEGEND manipulation
    # legend_manipulation(ax, a_sim_info.info_tr())
    legend_manipulation(ax, "")

    # close & save figure
    close_fig(out_dir + out_file, fig, save=save, show=show, use_z_eff=use_z_eff)
    

def plot_corr_func_universal(r, xi, r_lin, xi_lin, r_nl, xi_nl, lab, suptitle, ylabel,
                             figtext, out_dir, file_name, save, show, r2, extra_data=None, use_z_eff=False, chi=False):
    
    z_out = lab if lab == 'init' else 'z' + lab[4:]
    fig = plt.figure(figsize=fig_size)

    if extra_data is None: extra_data = []

    # check for r2 multiplier
    mlt = mlt_lin = mlt_nl = 1
    if r2:
        mlt = r*r
        if xi_lin is not None: mlt_lin = r_lin*r_lin
        if xi_nl is not None: mlt_nl = r_nl*r_nl
        ylabel = r"$r^2" + ylabel + r"(r)$"
        file_name = out_dir + '%s_r2_%s' % (file_name, z_out)
        plt.xscale("linear")
        plt.yscale("linear")
        for data in extra_data:
            data["mlt"] = data["r"]*data["r"]
    else:
        ylabel = r'$' + ylabel + r"(r)$"
        file_name = out_dir + '%s_%s' % (file_name, z_out)
        plt.xscale("log")
        plt.yscale("log")

    # plot all -- sim, lin, non-lin
    plt.plot(r, xi*mlt, 'o', label=lab)
    for data in extra_data:
        plt.plot(data["r"], data["xi"]*data["mlt"], 'o', label=data["lab"])
    if xi_lin is not None: plt.plot(r_lin, xi_lin*mlt_lin, '-', label=r"$\Lambda$CDM (lin)")
    if xi_nl is not None: plt.plot(r_nl, xi_nl*mlt_nl, '-', label=r"$\Lambda$CDM (nl)")

    # adjust figure, labels
    fig_suptitle(fig, suptitle)
    plt.xlabel(r"$r [h^{-1}{\rm Mpc}]$")
    plt.ylabel(ylabel)
    if chi:
        ext_legen = {'mlt_col' : 0.5, 'ncol' : 4}
    else:
        ext_legen = {'mlt_col' : 0.5, 'ncol' : 4}
    fig_leg = legend_manipulation(figtext="", ext_legen=ext_legen)
    fig.subplots_adjust(**subplt_adj_sym)

    # save & show (in jupyter)
    close_fig(file_name, fig, save=save, show=show, use_z_eff=use_z_eff, ext_legen=fig_leg)

def plot_corr_func_ratio(r, xi, r_lin, xi_lin, r_nl, xi_nl, lab, suptitle, ylabel,
                         figtext, out_dir, file_name, save, show, extra_data, peak_loc=None, use_z_eff=False, chi=False):
    # names
    z_out = lab if lab == 'init' else 'z' + lab[4:]
    ylabel = r'$' + ylabel + r"(r)/" + ylabel + r"_{\rm nl}(r)$"
    file_name = out_dir + '%s_ratio_%s' % (file_name, z_out)
    
    # check same lengths, validity of xi_n;
    if np.array_equal(r, r_nl):
        xi_an = xi_nl
        suptitle += r" $\Lambda$CDM (nl)"
    elif np.array_equal(r, r_lin):
        xi_an = xi_lin
        suptitle += r" $\Lambda$CDM (lin)"
    else:
        raise ValueError("Invalid values of radiues.")
    
    # figure
    fig = plt.figure(figsize=fig_size)
    ax = plt.gca()
    ax.yaxis.grid(True)
    ymin = 0.5
    ymax = 1.5
    ax.set_ylim(ymin,ymax)
    # ax.set_yscale('symlog', linthreshy=1e-3, linscaley=6)

    # plot ratio, linear
    if use_z_eff:
        sim = use_z_eff['sim']
        xi_an = power.corr_func(sim, z=use_z_eff['z'], non_lin=True)[1]
        xi_lin = power.corr_func(sim, z=use_z_eff['z'], non_lin=False)[1]

    ax.plot(r, xi_lin/xi_an, 'k-', label=r" $\Lambda$CDM (lin)")
    ax.plot(r, xi/xi_an, 'o', label=lab)

    # plot other data (if available)
    if extra_data is not None:
        for data in extra_data:
            xi_an = power.corr_func(sim, z=data['z_eff'], non_lin=True)[1]
            ax.plot(data['r'], data['xi']/xi_an, 'o', label=data['lab'])

    # plot BAO peak location (if available)
    if peak_loc is not None:
        ax.axvline(x=peak_loc, ls='--', color='k')

    # adjust figure, labels
    fig_suptitle(fig, suptitle)
    plt.xlabel(r"$r [h^{-1}{\rm Mpc}]$")
    plt.ylabel(ylabel)
    if chi:
        ext_legen = {'mlt_col' : 0.5, 'ncol' : 4}
    else:
        ext_legen = {'mlt_col' : 0.5, 'ncol' : 4}
    fig_leg = legend_manipulation(figtext="", ext_legen=ext_legen)
    fig.subplots_adjust(**subplt_adj_sym)

    # save & show (in jupyter)
    close_fig(file_name, fig, save=save, show=show, use_z_eff=use_z_eff, ext_legen=fig_leg)

def plot_corr_func_single(corr_data, lab, a_sim_info, corr_data_lin=None, corr_data_nl=None, out_dir='auto',
                          save=True, show=False, use_z_eff=False, is_sigma=False, only_r2=True, extra_data=None, peak_loc=None):
    if out_dir == 'auto':
        out_dir = a_sim_info.res_dir
    if is_sigma:
        suptitle = "Amplitude of density fluctuation"
        file_name = "sigma"
        ylabel = r"\sigma^2"
    else:
        suptitle = "Correlation function"
        file_name = "corr_func"
        ylabel = r"\xi"

    figtext = a_sim_info.info_tr()

    # modify labels if we are plotting multiple data
    if extra_data is not None:
        figtext = figtext.replace(a_sim_info.app + ": ", "")
        suptitle += ", " + lab
        lab = a_sim_info.app

    # get data
    r, xi = corr_data
    r_lin, xi_lin = corr_data_lin if corr_data_lin is not None else (None, None)
    r_nl, xi_nl = corr_data_nl if corr_data_nl is not None else (None, None)

    # first plot, xi(r)
    if not only_r2: plot_corr_func_universal(
        r, xi, r_lin, xi_lin, r_nl, xi_nl, lab, suptitle, ylabel, figtext,
        out_dir, file_name, save, show, False, extra_data, use_z_eff)

    # second plot, r*r*xi(r)
    plot_corr_func_universal(
        r, xi, r_lin, xi_lin, r_nl, xi_nl, lab, suptitle, ylabel, figtext,
        out_dir, file_name, save, show, True, extra_data, use_z_eff)

    # third plot, xi(r)/xi_lin/nl
    plot_corr_func_ratio(
        r, xi, r_lin, xi_lin, r_nl, xi_nl, lab, suptitle, ylabel, figtext,
        out_dir, file_name, save, show, extra_data, peak_loc, use_z_eff)

# correlation function stacked data, linear and emu corr. func in files
def plot_corr_func(corr_data_all, zs, a_sim_info, out_dir='auto', save=True, show=False, use_z_eff=False,
                   is_sigma=False, only_r2=True, extra_data=None, peak_loc=None):
    for lab, corr_par, corr_lin, corr_nl in iter_data(zs, [corr_data_all['par'],
                                                      corr_data_all['lin'], corr_data_all['nl']]):
        plot_corr_func_single(
            corr_par, lab, a_sim_info,
            corr_data_lin=corr_lin, corr_data_nl=corr_nl, out_dir=out_dir,
            save=save, show=show, is_sigma=is_sigma, only_r2=only_r2,
            extra_data=extra_data, peak_loc=peak_loc, use_z_eff=use_z_eff)

def plot_peak_uni(a_sim_info, ax, bao_type, idx, use_z_eff=False, ls=None, get_last_col=False, fp_comp=False, single=False):
    # load all available data (GSL integration could have failed)
    peak_data = [x for x in a_sim_info.data["corr_func"]["par_peak"] if x["z"] != "init"]
    zs = [x["z"] for x in peak_data]
    a = [1./(1+z) for z in zs]

    # location / amplitude / width of the BAO peak, label
    data = np.array([x["popt"][idx] for x in peak_data])
    data_err = np.array([x["perr"][idx] for x in peak_data])

    label = a_sim_info.app
    if label == 'CHI':
        label = get_chi_label(a_sim_info, single=single)

    # get last used color
    color = ax.get_lines()[-1].get_color() if get_last_col else None
    
    # comparison to the non-linear prediction at z_eff
    if use_z_eff:
        corr = [power.corr_func(use_z_eff['sim'], z=z, non_lin=True) for z in use_z_eff['z']]
        data_nl = np.array([power.get_bao_peak(x)["popt"][idx] for x in corr])

    # comparison to FP
    elif fp_comp:
        peak_data_nl = [x for x in fp_comp.data["corr_func"]["par_peak"] if x["z"] != "init"]
        data_nl = np.array([x["popt"][idx] for x in peak_data_nl])

    # comparison to the non-linear prediction
    else:
        peak_data_nl = [x for x in a_sim_info.data["corr_func"]["nl_peak"] if x["z"] != "init"]
        data_nl = np.array([x["popt"][idx] for x in peak_data_nl])

    # plot simulation peak
    # ax.set_yscale('symlog', linthreshy=0.01, linscaley=0.5)
    if not single:
        label += ' (%s)' % bao_type
    ax.errorbar(a, data / data_nl - 1, yerr=data_err / data_nl, ls=ls, label=label, color=color)
    

def plot_peak_loc(a_sim_info, ax, use_z_eff=False, get_last_col=False, fp_comp=False, single=False):
    """ plot peak location to the given axis """
    plot_peak_uni(a_sim_info, ax, "loc", 1, use_z_eff=use_z_eff, ls='-', get_last_col=get_last_col, fp_comp=fp_comp, single=single)

def plot_peak_amp(a_sim_info, ax, use_z_eff=False, get_last_col=False, fp_comp=False, single=False):
    """ plot peak amplitude to the given axis """
    plot_peak_uni(a_sim_info, ax, "amp", 0, use_z_eff=use_z_eff, ls=':', get_last_col=get_last_col, fp_comp=fp_comp, single=single)

def plot_peak_width(a_sim_info, ax, use_z_eff=False, get_last_col=False, fp_comp=False, single=False):
    """ plot peak amplitude to the given axis """
    plot_peak_uni(a_sim_info, ax, "width", 2, use_z_eff=use_z_eff, ls='--', get_last_col=get_last_col, fp_comp=fp_comp, single=single)

def plot_corr_peak(sim_infos, out_dir='auto', save=True, show=False, use_z_eff=False, plot_loc=True, plot_amp=True, plot_width=True, fp_comp=False, single=False, chi=False):
    # output
    if out_dir == 'auto':
        if len(sim_infos) == 1:
            out_dir = sim_infos[0].res_dir
        else:
            out_dir = report_dir
    out_file = "corr_peak"
    if fp_comp:
        out_file += '_fp_chi'
    if plot_loc:
        out_file += "_loc"
    if plot_amp:
        out_file += "_amp"
    if plot_width:
        out_file += "_width"

    # figure
    fig = plt.figure(figsize=fig_size)
    ax = plt.gca()
    ax.yaxis.grid(True)

    for i, a_sim_info in enumerate(sim_infos):
        # use effective redshift
        z_eff = use_z_eff[i] if use_z_eff else False

        # peak location
        if plot_loc:
            plot_peak_loc(a_sim_info, ax, use_z_eff=z_eff, fp_comp=fp_comp, single=single)

        # peak amplitude
        if plot_amp:
            get_last_col = plot_loc
            plot_peak_amp(a_sim_info, ax, use_z_eff=z_eff, get_last_col=get_last_col, fp_comp=fp_comp, single=single)

        # peak width
        if plot_width:
            get_last_col = plot_loc or plot_amp
            plot_peak_width(a_sim_info, ax, use_z_eff=z_eff, get_last_col=get_last_col, fp_comp=fp_comp, single=single)

    # linear prediction
    if single:
        sim = a_sim_info.sim
        a = np.linspace(1./201, 1, num=100)
        zs = 1/a - 1
        if plot_loc: idx = 1
        elif plot_amp: idx = 0
        elif plot_width: idx = 2
        corr = [power.corr_func(sim, z=z, non_lin=False) for z in zs]
        corr_nl = [power.corr_func(sim, z=z, non_lin=True) for z in zs]
        data_lin = np.array([power.get_bao_peak(x)["popt"][idx] for x in corr])
        data_nl = np.array([power.get_bao_peak(x)["popt"][idx] for x in corr_nl])
        ax.plot(a, data_lin/data_nl - 1, 'k-', label=r"$\Lambda$CDM (lin)")

    # labels
    plt.xlabel(r"$a$")
    fig_suptitle(fig, "Relative BAO peak location and amplitude")

    # LEGEND manipulation
    if chi:
        ext_legen = {'mlt_col' : 0.6, 'ncol' : 2}
    else:
        ext_legen = {'mlt_col' : 0.5, 'ncol' : 4}
    fig_leg = legend_manipulation(figtext="", ext_legen=ext_legen)
    fig.subplots_adjust(**subplt_adj_sym)

    # close & save figure
    close_fig(out_dir + out_file, fig, save=save, show=show, use_z_eff=use_z_eff, ext_legen=fig_leg)

def plot_eff_time_ax(a_sim_info, ax, a_eff_type="Pk"):
    # ZA and TZA do not have non-linear power spectra
    if not a_sim_info.data["eff_time"][a_eff_type]:
        return

    # extract variables
    a = a_sim_info.data["eff_time"][a_eff_type]['a']
    # D_eff = a_sim_info.data["eff_time"][a_eff_type]['D_eff']
    D_eff_ratio = a_sim_info.data["eff_time"][a_eff_type]['D_eff_ratio']
    a_err = a_sim_info.data["eff_time"][a_eff_type]['a_err']
    label = a_sim_info.app # +  '$: L = %i$ Mpc/h' % a_sim_info.box_opt["box_size"]

    # plot
    if a_eff_type == "sigma_R" or a_eff_type == "Pk":
        ax.plot(a, D_eff_ratio, 'o-', label=label)
    elif a_eff_type == "Pk_nl":
        ax.errorbar(a, D_eff_ratio, 'o-', yerr=a_err, label=label)
    color = ax.get_lines()[-1].get_color()

    # spline
    try:
        func = lambda x, d, b, c : d + b*x*np.exp(-c*x)
        spl = MySpline(a, D_eff_ratio, func, p0=(1,0,0))
        a_spl = np.linspace(a[0], a[-1], 100)
        ax.plot(a_spl, spl(a_spl), '--', color=color)
    except:
        pass


def plot_eff_time(stack_infos, out_dir='auto', a_eff_type="Pk", save=True, show=False, use_z_eff=False, verbose=False):
    # plot everything
    if a_eff_type == 'all':
        plot_eff_time(stack_infos, out_dir=out_dir, a_eff_type="sigma_R", save=save, show=show, use_z_eff=use_z_eff)
        plot_eff_time(stack_infos, out_dir=out_dir, a_eff_type="Pk", save=save, show=show, use_z_eff=use_z_eff)
        # plot_eff_time(stack_infos, out_dir=out_dir, a_eff_type="Pk_nl", save=save, show=show, use_z_eff=use_z_eff)
        return

    # output
    if out_dir == 'auto':
        if len(stack_infos) == 1:
            out_dir = stack_infos[0].res_dir
        else:
            out_dir = report_dir

    # figure
    fig = plt.figure(figsize=fig_size)
    ax = plt.gca()
    
    if verbose:
        ut.print_info("Effective time: '%s'" % a_eff_type)
    for stack_info in stack_infos:
        plot_eff_time_ax(stack_info, ax, a_eff_type)
    
    ax.set_ylabel(r'$D_{eff}/D_{GR}$')
    ax.set_xlabel(r'$a$')
    ax.legend()
    fig.subplots_adjust(**subplt_adj_sym)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.grid(True)

    # adjust extreme values
    adjust_extreme_values(ax, 0.9, 1.1)
    
    # save & show (in jupyter)
    close_fig(out_dir + 'D_eff_' + a_eff_type, fig, save=save, show=show, use_z_eff=use_z_eff)

def plot_pwr_spec_nl_amp_ax(a_sim_info, ax, ymin=0, ymax=0.4):
    # extract variables
    a = power.get_a_fom_zs(a_sim_info.data["pk_nl_amp"]['z'])
    A = a_sim_info.data["pk_nl_amp"]['A']
    A_err = a_sim_info.data["pk_nl_amp"]['A_err']
    label = a_sim_info.app +  '$: L = %i$ Mpc/h' % a_sim_info.box_opt["box_size"]

    # plot
    ax.errorbar(a, A, yerr=A_err, label=label)

    # ymin = min(ymin, np.min(A))
    # ymax = min(ymax, np.max(A))
    ax.set_ylim(ymin=ymin*0.9, ymax=ymax*1.1)

def plot_pwr_spec_nl_amp(stack_infos, out_dir='auto', save=True, show=False):
    # output
    if out_dir == 'auto':
        if len(stack_infos) == 1:
            out_dir = stack_infos[0].res_dir
        else:
            out_dir = report_dir

    # figure
    fig = plt.figure(figsize=fig_size)
    ax = plt.gca()
    
    for stack_info in stack_infos:
        plot_pwr_spec_nl_amp_ax(stack_info, ax)
    
    ax.set_ylabel(r'$A_{\rm nl}$')
    ax.set_xlabel(r'$a$')
    ax.legend()
    fig.subplots_adjust(**subplt_adj_sym)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.grid(True)
    
    # save & show (in jupyter)
    close_fig(out_dir + 'pwr_spec_nl_amp', fig, save=save, show=show)

def plot_eff_growth_rate_ax(a_sim_info, ax, a_eff_type="Pk"):
    # ZA and TZA do not have non-linear power spectra
    if not a_sim_info.data["eff_time"][a_eff_type]:
        return

    # extract variables
    a = np.array(a_sim_info.data["eff_time"][a_eff_type]['a'])
    a_spl = np.linspace(a[0], a[-1], 100)
    D_eff = a_sim_info.data["eff_time"][a_eff_type]['D_eff']
    D_eff_ratio = a_sim_info.data["eff_time"][a_eff_type]['D_eff_ratio']

    f = np.diff(np.log(D_eff)) /np.diff(np.log(a))
    
    # smooth before derivative
    func = lambda x, d, b, c : d + b*x*np.exp(-c*x)
    spl_ = MySpline(a, D_eff_ratio, func, p0=(1,0,0))
    spl = lambda x : spl_(x) * power.growth_factor(x, a_sim_info.sim.cosmo)

    # diff in midpoints
    a_spl = (a_spl[1:] + a_spl[:-1])/2.
    a = (a[1:] + a[:-1])/2.

    slope = get_slope(a_spl, spl, dx=0.05)

    # linear prediction for growth rate
    f_lin = power.growth_rate(a, a_sim_info.sim.cosmo)
    f_lin_spl = power.growth_rate(a_spl, a_sim_info.sim.cosmo)

    # plot
    label = a_sim_info.app # +  '$: L = %i$ Mpc/h' % a_sim_info.box_opt["box_size"]
    ax.plot(a, f/f_lin, 'o-', label=label)
    color = ax.get_lines()[-1].get_color()
    ax.plot(a_spl, slope/f_lin_spl, '--', color=color)

def plot_eff_growth_rate(stack_infos, out_dir='auto', a_eff_type="Pk", save=True, show=False, use_z_eff=False):
    # plot everything
    if a_eff_type == 'all':
        plot_eff_growth_rate(stack_infos, out_dir=out_dir, a_eff_type="sigma_R", save=save, show=show, use_z_eff=use_z_eff)
        plot_eff_growth_rate(stack_infos, out_dir=out_dir, a_eff_type="Pk", save=save, show=show, use_z_eff=use_z_eff)
        # plot_eff_growth_rate(stack_infos, out_dir=out_dir, a_eff_type="Pk_nl", save=save, show=show, use_z_eff=use_z_eff)
        return

    # output
    if out_dir == 'auto':
        if len(stack_infos) == 1:
            out_dir = stack_infos[0].res_dir
        else:
            out_dir = report_dir

    # figure
    fig = plt.figure(figsize=fig_size)
    ax = plt.gca()
    
    ut.print_info("Growth rate: '%s'" % a_eff_type)
    for stack_info in stack_infos:
        plot_eff_growth_rate_ax(stack_info, ax, a_eff_type)
    
    ax.set_ylabel(r'$f/f_{GR}$')
    ax.set_xlabel(r'$a$')
    ax.legend()
    fig.subplots_adjust(**subplt_adj_sym)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.grid(True)

    # adjust extreme values
    adjust_extreme_values(ax, 0.8, 1.1)
    
    # save & show (in jupyter)
    close_fig(out_dir + 'f_eff_' + a_eff_type, fig, save=save, show=show, use_z_eff=use_z_eff)

def plot_pwr_spec_diff_from_data_ax(ax, data_list, zs, a_sim_info, show_scales=True, pk_type='dens', show_nyquist=True, max_nyquist=False):
    if pk_type == 'chi':
        # transform chameleon power spectrum to suppression
        for z, data in izip(zs, data_list):
            a, k, Pk = 1/(1.+z), data[0], data[1]
            data[1] = power.chi_trans_to_supp(a, k, Pk, a_sim_info.sim.cosmo, a_sim_info.chi_opt)
            if len(data) == 3:
                data[2] = power.chi_trans_to_supp(a, k, data[2], a_sim_info.sim.cosmo, a_sim_info.chi_opt)
        # transform supp (ref: lin) to supp (ref: init)
        power.chi_trans_to_init(data_list)

    ax.set_xscale('log')

    ymin = ymax = 0

    # SMALL / MEDIUM / LARGE SCALE VALUES
    # half of nyquist wavelength, 7 subintervals
    k = data_list[-1][0]
    idx = (np.abs(k - 0.5 * a_sim_info.k_nyquist["particle"])).argmin() / 7
    cmap = cm.get_cmap('gnuplot')
    if show_scales:
        ax.axvspan(k[0 * idx], k[1 * idx], alpha=0.2, color=cmap(0.1))
        ax.axvspan(k[3 * idx], k[4 * idx], alpha=0.3, color=cmap(0.5))
        ax.axvspan(k[6 * idx], k[7 * idx], alpha=0.4, color=cmap(0.9))

    for lab, data, a in iter_data(zs, [data_list], get_a=True):
        k, P_k = data[0], data[1]
        P_k_std = data[2] if len(data) == 3 else None

        if P_k_std is None:
            plt.plot(k, P_k, 'o', label=lab)
            ymax = max(ymax, np.max(P_k[0:7 * idx]))
            ymin = min(ymin, np.min(P_k[0:7 * idx]))
        else:
            plt.errorbar(k, P_k, fmt='o', yerr=P_k_std, label=lab)
            ymax = max(ymax, np.max(P_k[0:7 * idx] + P_k_std[0:7 * idx]))
            ymin = min(ymin, np.min(P_k[0:7 * idx] - P_k_std[0:7 * idx]))

    if show_nyquist:
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

def plot_pwr_spec_diff_from_data_mlt(data_lists, zs, sim_infos, out_dir='auto', show_scales=False, pk_type='dens',
                                 ext_title='init', save=True, show=False, use_z_eff=False):
    # output
    if out_dir == 'auto':
        if len(sim_infos) == 1:
            out_dir = sim_infos[0].res_dir
        else:
            out_dir = report_dir

    if pk_type == "dens":
        out_file = 'pwr_spec_diff'
        suptitle = "Power spectrum difference"
    elif pk_type == "vel":
        out_file = 'vel_pwr_spec_diff'
        suptitle = r"Power spectrum difference $(\nabla\cdot u)$"
    elif pk_type == 'chi':
        out_file = 'pwr_spec_diff_chi'
        suptitle = "Chameleon power spectrum difference"
        ext_title = 'init'
            
    out_file += '_%s' % ext_title
    suptitle += ' (ref: %s)' % ext_title

    # size
    x_fig_size, y_fig_size = fig_size
    x_fig_size *= 2
    y_fig_size *= 2
    fig = plt.figure(figsize=(x_fig_size, y_fig_size))

    gs = gridspec.GridSpec(2, 2)
    for gs_cur, data_list, a_sim_info in izip(gs, data_lists, sim_infos):
        ax = plt.subplot(gs_cur)
        plot_pwr_spec_diff_from_data_ax(ax, data_list, zs, a_sim_info, show_scales=show_scales, pk_type=pk_type)
    
    # adjust in case of overflow
    ymin, ymax = ax.get_ylim()
    if ymin < -1: ymin = -1
    ax.set_ylim(ymin, ymax)

    fig_suptitle(fig, suptitle)

    ax.set_xlabel(r"$k [h{\rm Mpc}^{-1}]$")
    ax.set_ylabel(r"$P(k)/P_{\rm lin}(k)-1$")

    # legend_manipulation(ax, a_sim_info.info_tr())
    legend_manipulation(ax, "")
    close_fig(out_dir + out_file, fig, save=save, show=show, use_z_eff=use_z_eff)

def plot_pwr_spec_diff_from_data(data_list, zs, a_sim_info, out_dir='auto', show_scales=True, pk_type='dens',
                                 ext_title='par', save=True, show=False, use_z_eff=False, add_app=False, show_nyquist=True, max_nyquist=False, chi=False):
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
        ext_title = 'init'
            
    out_file += '_%s' % ext_title
    suptitle += ' (ref: %s)' % ext_title

    if add_app:
        out_file += '_%s' % a_sim_info.app

    fig = plt.figure(figsize=fig_size)
    ax = plt.gca()
    plot_pwr_spec_diff_from_data_ax(ax, data_list, zs, a_sim_info, show_scales=show_scales, pk_type=pk_type,
                                    show_nyquist=show_nyquist, max_nyquist=max_nyquist)
    
    # adjust in case of overflow
    ymin, ymax = ax.get_ylim()
    if ymin < -1: ymin = -1
    ax.set_ylim(ymin, ymax)

    fig_suptitle(fig, suptitle)
    ax.set_xlabel(r"$k [h{\rm Mpc}^{-1}]$")
    ax.set_ylabel(r"$P(k)/P_{\rm lin}(k)-1$")
    # legend_manipulation(ax, a_sim_info.info_tr())
    if chi:
        ext_legen = {'mlt_col' : 0.5, 'ncol' : 4}
    else:
        ext_legen = {'mlt_col' : 0.7, 'ncol' : 4}
    fig_leg = legend_manipulation(ax, "", ext_legen=ext_legen)
    close_fig(out_dir + out_file, fig, save=save, show=show, use_z_eff=use_z_eff, ext_legen=fig_leg)


def plot_pwr_spec_diff_map_from_data(data_array, zs, a_sim_info, out_dir='auto', add_app=False, pk_type='dens',
                                     ext_title='', save=True, show=False, shading='flat',
                                     show_nyquist=True, use_z_eff=False, vmin=-1, vmax=1):
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

    if ext_title == '':
        out_file += '_map'
    else:
        out_file += '_%s_map' % ext_title

    suptitle += ' (ref: %s)' % ext_title

    if add_app:
        out_file += '_%s' % a_sim_info.app

    fig = plt.figure(figsize=fig_size_map)
    gs = gridspec.GridSpec(1, 60, wspace=0.5)
    ax = plt.subplot(gs[0, : -5])
    cbar_ax = plt.subplot(gs[0, -4 :])

    ax.set_xscale('log')
    a = [1 / (1 + z) for z in zs]

    # hack around pcolormesh plotting edges
    if shading == 'flat':
        if len(a) == 1:
            da = 2*a[0]
        else:
            da = (a[-1] - a[0]) / (len(a) - 1)
        a = np.array([a[0]-da/2] + [1 / (1 + z) + da/2 for z in zs])
    
    k = data_array[0][0]
    supp = data_array[:, 1, :] # extract Pk, shape = (zs, k)
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

    im = ax.pcolormesh(k, a, supp, cmap='seismic', shading=shading,
        norm=SymLogNorm(linthresh=linthresh, linscale=linscale, vmin=vmin, vmax=vmax))
    cbar = fig.colorbar(im, cax=cbar_ax, ticks=ticks)
    cbar.ax.set_yticklabels(labels)

    if show_nyquist and a_sim_info.k_nyquist is not None:
        ls = [':', '-.', '--']
        ls *= (len(a_sim_info.k_nyquist) - 1) / 3 + 1
        ls = iter(ls)
        val_set = set(a_sim_info.k_nyquist.itervalues())
        for val in val_set:
            ax.axvline(val, ls=next(ls), c='k')

    fig_suptitle(fig, suptitle)
    ax.set_xlabel(r"$k [h{\rm Mpc}^{-1}]$")
    ax.set_ylabel(r"$a(t)$")
    plt.draw()
    plt.figtext(0.5, 0.95, "",
                bbox={'facecolor': 'white', 'alpha': 0.2}, size=14, ha='center', va='top')
    fig.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9)

    close_fig(out_dir + out_file, fig, save=save, show=show, use_z_eff=use_z_eff, format='png')

def plot_supp(sim_infos, out_dir, suptitle='', save=True, show=False, use_z_eff=False, scale='', show_k_lms=False, res=None):
    fig = plt.figure(figsize=fig_size)
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
            print("WARNING! Unknown scale ='%s'. Skipping." % scale)
            return None
        plt.plot(a, supp, '-o',
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
    fig_suptitle(fig, "Power spectrum suppression" + suptitle)
    plt.xlabel(r"$a(t)$")
    ylabel = r"$\langle{P(k)/P_{\rm lin}(k)-1}\rangle$"
    if res is not None:
        ylabel += r', residual from $\nu=%.1f$' % res.nu
    plt.ylabel(ylabel)
    legend_manipulation()
    close_fig(out_dir + 'supp', fig, save=save, show=show, use_z_eff=use_z_eff)


def plot_dens_histo(data_list, zs, a_sim_info, out_dir='auto', fix_N=1, fix_rho=0.0, save=True, show=False, use_z_eff=False):
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

    fig_suptitle(fig, "Overdensity distribution")

    plt.figtext(0.5, 0.95, a_sim_info.info_tr(),
                bbox={'facecolor': 'white', 'alpha': 0.2}, size=14, ha='center', va='top')
    close_fig(out_dir + 'dens_histo', fig, save=save, show=show, use_z_eff=use_z_eff)


def plot_par_last_slice(files, files_t, zs, a_sim_info, out_dir='auto', save=True, show=False, use_z_eff=False):
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
    num_steps = len(x_t) // num_track
    x_t = [x_t[i:i + num_steps] for i in xrange(0, len(x_t), num_steps)]
    y_t = [y_t[i:i + num_steps] for i in xrange(0, len(y_t), num_steps)]

    ax.plot(x, y, 'ob', ms=1)
    for i in xrange(num_track):
        ax.plot(x_t[i], y_t[i], '--or', ms=4, lw=1.5,
                markevery=(num_steps - 1, num_steps))

    plt.figtext(0.5, 0.94, a_sim_info.info_tr(),
                bbox={'facecolor': 'white', 'alpha': 0.2}, size=14, ha='center', va='top')
    ax.set_xlabel(r"$x [h^{-1}{\rm Mpc}]$")
    ax.set_ylabel(r"$z [h^{-1}{\rm Mpc}]$")
    fig_suptitle(fig, "Slice through simulation box (particles), z = %.2f" % zs[-1])
    close_fig(out_dir + 'par_evol_last', fig, save=save, show=show, use_z_eff=use_z_eff)


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
    ax.set_xlabel(r"$x [h^{-1}{\rm Mpc}]$")
    ax.set_ylabel(r"$z [h^{-1}{\rm Mpc}]$")
    del x, y, data

    def animate(j):
        if j < num:
            i = j
        else:
            i = 2 * num - j - 1
        fig_suptitle(fig, "Slice through simulation box (particles), z = %.2f" % zs[i])
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

def plot_dens_one_slice(rho, z, a_sim_info, out_dir='auto', save=True, show=False, use_z_eff=False):
    if out_dir == 'auto':
        out_dir = a_sim_info.res_dir
    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(1, 15, wspace=0.5)
    ax = plt.subplot(gs[0, : -1])
    cbar_ax = plt.subplot(gs[0, -1])

    plt.figtext(0.5, 0.94, a_sim_info.info_tr(),
                bbox={'facecolor': 'white', 'alpha': 0.2}, size=14, ha='center', va='top')
    ax.set_xlabel(r"$x [h^{-1}{\rm Mpc}]$")
    ax.set_ylabel(r"$z [h^{-1}{\rm Mpc}]$")
    L = int(np.sqrt(rho.shape[0]))
    rho.shape = L, L
    im = ax.imshow(rho, interpolation='bicubic', cmap='gnuplot',
                   norm=SymLogNorm(linthresh=1.0, linscale=1,
                                   vmin=-1, vmax=100), aspect='auto',
                   extent=[0, a_sim_info.box_opt["box_size"], 0, a_sim_info.box_opt["box_size"]])
    fig_suptitle(fig, "Slice through simulation box (overdensity), z = %.2f" % z)
    cbar = fig.colorbar(im, cax=cbar_ax, ticks=[-1, 0, 1, 10, 100])
    cbar.ax.set_yticklabels(['-1', '0', '1', '10', '> 100'])
    close_fig(out_dir + 'dens_z%.2f' % z, fig, save=save, show=show, use_z_eff=use_z_eff, format='png')

def plot_dens_two_slices(files, zs, a_sim_info, out_dir='auto', save=True, show=False, use_z_eff=False):
    half = len(files) // 2
    rho, z = np.loadtxt(files[half])[:, 2], zs[half]
    plot_dens_one_slice(rho, z, a_sim_info,
                        out_dir=out_dir, save=save, show=show, use_z_eff=use_z_eff)
    rho, z = np.loadtxt(files[-1])[:, 2], zs[-1]
    plot_dens_one_slice(rho, z, a_sim_info,
                        out_dir=out_dir, save=save, show=show, use_z_eff=use_z_eff)


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
    ax.set_xlabel(r"$x [h^{-1}{\rm Mpc}]$")
    ax.set_ylabel(r"$z [h^{-1}{\rm Mpc}]$")

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
        fig_suptitle(fig, "Slice through simulation box (overdensity), z = %.2f" % zs[i])
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


def plot_chi_evol(zs, a_sim_info, chi_opt=None, out_dir='auto', save=True, show=False, use_z_eff=False, k_scr=False):
    """" Plot evolution of chameleon background values -- Compton wavelength and screening potential """
    if out_dir == 'auto':
        out_dir = a_sim_info.res_dir
    out_file = 'chi_evol'
    suptitle = "Evolution of Chameleon"

    N = 4 if k_scr else 3

    # adjust size for 4 subplots
    x_fig_size, y_fig_size = fig_size
    y_fig_size *= N/3.

    fig = plt.figure(figsize=(x_fig_size, y_fig_size))
    cosmo = a_sim_info.sim.cosmo
    if chi_opt is None:
        chi_opt = [a_sim_info.chi_opt]
        
    ax1 = plt.subplot(N, 1, 1)
    ax2 = plt.subplot(N, 1, 2, sharex=ax1)
    ax3 = plt.subplot(N, 1, 3, sharex=ax1)
    
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax3.set_yscale('log')

    ax1.set_ylabel(r"$\lambda_C\, [h^{-1}{\rm Mpc}]$")
    ax2.set_ylabel(r"$\chi/M_{\rm Pl}$")
    ax3.set_ylabel(r"$\phi_{\rm scr}$")

    ax1.tick_params(axis='both', which='major')
    ax2.tick_params(axis='both', which='major')
    ax3.tick_params(axis='both', which='major')

    if k_scr:
        ax4 = plt.subplot(N, 1, 4, sharex=ax1)
        ax4.set_yscale('log')
        ax4.set_ylabel(r"$k_{\rm scr}\, [h{\rm Mpc}^{-1}]$")
        ax4.set_xlabel(r"z")
        ax4.tick_params(axis='both', which='major')
    else:
        ax3.set_xlabel(r"z")
    
    zs = [z for z in zs if z != 'init']
    a = [1./(z+1) for z in zs]
    
    linestyle_cycler = iter(['-','--',':','-.'])
    for chi in chi_opt:
        # chameleon parameters
        wavelengths = np.array([power.chi_compton_wavelength(a_, cosmo, chi) for a_ in a])
        psi_a = [power.chi_psi_a(a_, chi) for a_ in a]
        chi_a = [power.chi_bulk_a(a_, chi, CHI_A_UNITS=False) for a_ in a]

        # plots
        ls = linestyle_cycler.next()
        ax1.plot(zs, 1/wavelengths, ls=ls, label=r"$\Phi_{\rm scr} = 10^{%i}$, $n=%.1f$" % (np.log10(chi["phi"]), chi["n"]))
        ax2.plot(zs, chi_a, ls=ls)
        ax3.plot(zs, psi_a, ls=ls)
        if k_scr:
            k_scr_data = power.chi_psi_k_a(a, cosmo, chi)
            ax4.plot(zs, k_scr_data)
    
    fig_suptitle(fig, suptitle)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    
    # legend
    legend_manipulation(ax=ax1, loc='upper right')

    # subplots
    fig.subplots_adjust(hspace=0, **subplt_adj_sym)

    # close & save figure
    close_fig(out_dir + out_file, fig, save=save, show=show, use_z_eff=use_z_eff)


def plot_supp_lms(supp, a, a_sim_info, out_dir='auto', pk_type='dens', suptitle='', save=True, show=False,
                  add_app=False, scale_in_leg=True, use_z_eff=False):
    if out_dir == 'auto':
        out_dir = a_sim_info.res_dir
    if pk_type == "dens":
        out_file = 'supp'
        suptitle = "Power spectrum suppression"
    elif pk_type == "vel":
        out_file = 'supp_vel'
        suptitle = r"Power spectrum suppression $(\nabla\cdot u)$"
    elif pk_type == 'chi':
        out_file = 'supp_chi'
        suptitle = "Chameleon power spectrum suppression"

    fig = plt.figure(figsize=fig_size)
    ax = plt.gca()
    cmap = cm.get_cmap('gnuplot')

    for i, scale in enumerate(['Large', 'Medium', 'Small']):
        if scale_in_leg:
            label = '%s-scale:\n' r'$\langle%.2f,%.2f\rangle$' % (scale, supp[i][2][0], supp[i][2][0])
        else:
            label = None
            print('%s-scale: %.4f,%.4f' % (scale, supp[i][2][0], supp[i][2][0]))
        ax.errorbar(a, supp[i][0], fmt='-o', yerr=supp[i][1],
                    color=cmap(0.1+i*0.4), lw=5-i*1.5,
                    label=label)

    ymin, ymax = ax.get_ylim()
    if pk_type != 'chi' and ymax > 1:
        ymax = 1

    ymax = math.ceil(ymax / 0.1) * 0.1
    ymin = math.floor(ymin / 0.1) * 0.1
    plt.ylim(ymin=ymin, ymax=ymax)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.grid(True)

    fig_suptitle(fig, suptitle)
    plt.xlabel(r"$a(t)$")
    plt.ylabel(
        r"$\langle{P(k)/P_{\rm lin}(k)-1}\rangle$")
    
    # legend
    if scale_in_leg:
        # legend_manipulation(figtext=a_sim_info.info_tr())
        legend_manipulation(figtext="")

    if add_app:
        out_file += '_%s' % a_sim_info.app
    close_fig(out_dir + out_file, fig, save=save, show=show, use_z_eff=use_z_eff)


def plot_all_single_supp(res, out_dir='/home/michal/Documents/GIT/Adhesion-Approximation/output/supp_comparison/',
                         Nm=0, Np=0, L=0, nu=0, rs=0, app=''):
    subfiles = res.get_subfiles(Nm=Nm, Np=Np, L=L, nu=nu, rs=rs, app=app)
    for a_sim_info in subfiles:
        res.load_k_supp(a_sim_info)
        plot_supp_lms(a_sim_info.supp, a_sim_info.a, a_sim_info, show=True)


from matplotlib.patches import Ellipse

def get_err_ell(ax, opt, cov):
    if opt.shape != (2,): raise IndexError("'opt' argument has wrong shape")
    if cov.shape != (2,2): raise IndexError("'cov' argument has wrong shape")
    x, y = opt[0], opt[1]
    lambda_, v = np.linalg.eig(cov)
    lambda_ = np.sqrt(lambda_)
    height = lambda_[1]*2
    width = lambda_[0]*2
    angle = np.rad2deg(np.arccos(v[0, 0]))
    ell = Ellipse(xy=(x, y), width=width, height=height, angle=angle,
                  edgecolor='k', facecolor='none')
    ax.add_artist(ell)
    ax.plot(x, y, 'ko')