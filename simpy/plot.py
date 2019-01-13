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
import numpy as np
try:
    from itertools import izip # python 2
except ImportError:
    izip = zip # python 3
from scipy.misc import derivative

from . import power

matplotlib.rcParams['legend.numpoints'] = 1
label_size = 20
suptitle_size = 25
# fig_size = (15, 11)
fig_size = (14, 9)
subplt_adj_sym = {'left' : 0.15, 'right' : 0.95, 'bottom' : 0.15, 'top' : 0.95}
matplotlib.rcParams.update({'font.size': 15})
report_dir = "/home/michal/Documents/GIT/FastSim/report/plots/"

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

def close_fig(filename, fig, save=True, show=False, dpi=100, use_z_eff=False):
    """save and/or show figure, close figure"""
    if use_z_eff:
        filename += 'z_eff.png'
    else:
        filename += '.png'
    if save:
        fig.savefig(filename, dpi=dpi)
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

def legend_manipulation(ax=None, figtext="", loc='upper left', bbox_to_anchor=(1.0, 1.0)):
    ax = plt.gca() if ax is None else ax
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles, labels, loc=loc,
               bbox_to_anchor=bbox_to_anchor, fontsize=14)
    plt.draw()
    if figtext != "":
        plt.figtext(0.5, 0.95, figtext,
                bbox={'facecolor': 'white', 'alpha': 0.2}, size=14, ha='center', va='top')
    plt.subplots_adjust(left=0.1, right=0.84, bottom=0.1, top=0.89)

def get_a_init_from_zs(zs):
    """ from list of redshifts returns initial scale factor, i.e. value after 'init' """
    for z in zs:
        if z != 'init':
            return 1/(1.+z)

def plot_pwr_spec(data, zs, a_sim_info, Pk_list_extrap, err=False,
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
    
    fig_suptitle(fig, suptitle)
    ax.set_xlabel(r"$k [h/$Mpc$]$", fontsize=label_size)
    ax.set_ylabel(r"$P(k) [$Mpc$/h)^3$]", fontsize=label_size)

    # LEGEND manipulation
    # legend_manipulation(ax, a_sim_info.info_tr())
    legend_manipulation(ax, "")

    # close & save figure
    close_fig(out_dir + out_file, fig, save=save, show=show, use_z_eff=use_z_eff)

def plot_pwr_spec_comparison(data, zs, labels, cosmo,
                  out_dir='auto', save=True, show=False, use_z_eff=False):
    """" Plot power spectrum -- points and extrapolated values,
    show 'true' linear Pk at the initial and final redshift """
    if out_dir == 'auto':
        out_dir = report_dir
    out_file = 'pwr_spec'
    suptitle = "Power spectrum"

    fig = plt.figure(figsize=fig_size)
    ax = plt.gca()
    ax.set_yscale('log')
    ax.set_xscale('log')

    for _, Pkk, lab in iter_data(zs, [data, labels]):
        k, P_k = Pkk[0], Pkk[1]
        ax.plot(k, P_k, 'o', ms=3, label=lab)
        # show 1 standard deviation
        P_k_std = Pkk[2]
        ax.fill_between(k, P_k - P_k_std, P_k + P_k_std,
                        facecolor='darkgrey', alpha=0.5)
        k = np.geomspace(k[0],k[-1])

    # plot non/linear power spectra
    a_0 = 1./(1.+zs[-1])
    
    P_0 = power.lin_pow_spec(a_0, k, cosmo)
    P_0_nl = power.non_lin_pow_spec(a_0, k, cosmo)
    ax.plot(k, P_0, '-', label=r"$\Lambda$CDM (lin)")
    ax.plot(k, P_0_nl, '-',  label=r"$\Lambda$CDM (nl)")
    
    fig_suptitle(fig, suptitle, y=0.95)
    ax.set_xlabel(r"$k [h/$Mpc$]$", fontsize=label_size)
    ax.set_ylabel(r"$P(k) [$Mpc$/h)^3$]", fontsize=label_size)

    # LEGEND manipulation
    legend_manipulation(ax, "", loc='best')
    plt.subplots_adjust(**subplt_adj_sym)

    # close & save figure
    close_fig(out_dir + out_file, fig, save=save, show=show, use_z_eff=use_z_eff)


def plot_chi_pwr_spec(data_list_chi, zs_chi, a_sim_info, err=False, out_dir='auto', save=True, show=False, use_z_eff=False):
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
        lines = ax.plot(k, P_k, 'o', ms=3, label=lab)
        color = lines[0].get_color()
        P_a = power.chi_lin_pow_spec(a, k, a_sim_info.sim.cosmo, a_sim_info.chi_opt)
        ax.plot(k, P_a, '-', color=color)
        if err:
            P_k_std = Pkk[2] / pow(chi_bulk_a_n, 2)
            ax.fill_between(k, P_k - P_k_std, P_k + P_k_std,
                             facecolor='darkgrey', alpha=0.5)

    add_nyquist_info(ax, a_sim_info)
    fig_suptitle(fig, suptitle)
    ax.set_xlabel(r"$k [h/$Mpc$]$", fontsize=label_size)
    ax.set_ylabel(r"$P(k) [$Mpc$/h)^3$]", fontsize=label_size)

    # LEGEND manipulation
    # legend_manipulation(ax, a_sim_info.info_tr())
    legend_manipulation(ax, "")

    # close & save figure
    close_fig(out_dir + out_file, fig, save=save, show=show, use_z_eff=use_z_eff)


def plot_chi_fp_map(data, zs, a_sim_info):
    pass

def plot_chi_fp_z(data_z, a_sim_info, phi_s, out_dir='auto', suptitle='auto', save=True, show=False, use_z_eff=False):
    if out_dir == 'auto':
        out_dir = a_sim_info.res_dir
    out_file = 'chi_pwr_diff_fp'
    # if suptitle == 'auto':
    #     suptitle = "Relative chameleon power spectrum"

    fig = plt.figure(figsize=fig_size)
    ax = plt.gca()
    plt.xscale('log')
    ymax = 1
    ymin = 0.95
    for data_chi, phi in izip(data_z, phi_s): # each chi
        k = data_chi[0]
        Pk = data_chi[1]
        std = data_chi[2]
        ymax = max(ymax, np.max(Pk))
        ax.errorbar(k, Pk, fmt='o', yerr=std, ms=3, label=r"$\Phi_{scr}=%.1e$" % phi)

    add_nyquist_info(ax, a_sim_info)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.grid(True)

    ymax *= 1.1
    ax.set_ylim(ymin, ymax)

    fig_suptitle(fig, suptitle)
    plt.xlabel(r"$k [h/$Mpc$]$", fontsize=label_size)
    plt.ylabel(r"${P_\chi(k)}/{P_{FP}(k)}$", fontsize=label_size)
    #figtext = a_sim_info.info_tr().replace("FP: ", "")
    legend_manipulation(ax, figtext="", loc='upper left', bbox_to_anchor=(0.0,1.0))
    plt.subplots_adjust(**subplt_adj_sym)
    close_fig(out_dir + out_file, fig, save=save, show=show, use_z_eff=use_z_eff)

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
        ax.plot(k_half, slope, 'o', ms=3, label=lab)
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
    ax.set_xlabel(r"$k [h/$Mpc$]$", fontsize=label_size)
    ax.set_ylabel(r"d$\ln P(k)/$d$\ln k$]", fontsize=label_size)

    # LEGEND manipulation
    # legend_manipulation(ax, a_sim_info.info_tr())
    legend_manipulation(ax, "")

    # close & save figure
    close_fig(out_dir + out_file, fig, save=save, show=show, use_z_eff=use_z_eff)
    

def plot_corr_func_universal(r, xi, r_lin, xi_lin, r_nl, xi_nl, lab, suptitle, ylabel,
                             figtext, out_dir, file_name, save, show, r2, extra_data=None, use_z_eff=False):
    
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
    plt.plot(r, xi*mlt, 'o', ms=3, label=lab)
    for data in extra_data:
        plt.plot(data["r"], data["xi"]*data["mlt"], 'o', ms=3, label=data["lab"])
    if xi_lin is not None: plt.plot(r_lin, xi_lin*mlt_lin, '-', label=r"$\Lambda$CDM (lin)")
    if xi_nl is not None: plt.plot(r_nl, xi_nl*mlt_nl, '-', label=r"$\Lambda$CDM (nl)")

    # adjust figure, labels
    fig_suptitle(fig, suptitle)
    plt.xlabel(r"$r [$Mpc$/h]$", fontsize=label_size)
    plt.ylabel(ylabel, fontsize=label_size)
    legend_manipulation(figtext="", loc='best')
    plt.subplots_adjust(**subplt_adj_sym)

    # save & show (in jupyter)
    close_fig(file_name, fig, save=save, show=show, use_z_eff=use_z_eff)

def plot_corr_func_ratio(r, xi, r_lin, xi_lin, r_nl, xi_nl, lab, suptitle, ylabel,
                         figtext, out_dir, file_name, save, show, extra_data, peak_loc=None, use_z_eff=False):
    # names
    z_out = lab if lab == 'init' else 'z' + lab[4:]
    ylabel = r'$' + ylabel + r"(r)$"
    file_name = out_dir + '%s_%s' % (file_name, z_out)
    
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
    ax.set_ylim(-0.5,0.5)

    # plot ratio
    plt.plot(r, xi/xi_an - 1, 'o', ms=3, label=lab)

    # plot other data (if available)
    if extra_data is not None:
        for data in extra_data:
            plt.plot(data['r'], data['xi']/xi_an - 1, 'o', ms=3, label=data['lab'])

    # plot BAO peak location (if available)
    if peak_loc is not None:
        ax.axvline(x=peak_loc, ls='--', color='k')

    # adjust figure, labels
    fig_suptitle(fig, suptitle)
    plt.xlabel(r"$r [$Mpc$/h]$", fontsize=label_size)
    plt.ylabel(ylabel, fontsize=label_size)
    legend_manipulation(figtext="", loc='best')
    plt.subplots_adjust(**subplt_adj_sym)

    # save & show (in jupyter)
    close_fig(file_name, fig, save=save, show=show, use_z_eff=use_z_eff)

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
def plot_corr_func(corr_data_all, zs, a_sim_info, out_dir='auto', save=True, show=False, use_z_eff=False, is_sigma=False, only_r2=True, extra_data=None, peak_loc=None):
    for lab, corr_par, corr_lin, corr_nl in iter_data(zs, [corr_data_all['par'],
                                                      corr_data_all['lin'], corr_data_all['nl']]):
        plot_corr_func_single(
            corr_par, lab, a_sim_info,
            corr_data_lin=corr_lin, corr_data_nl=corr_nl, out_dir=out_dir,
            save=save, show=show, is_sigma=is_sigma, only_r2=only_r2,
            extra_data=extra_data, peak_loc=peak_loc)

def plot_peak_loc(a, a_sim_info, ax, cut, only_nl=True):
    """ plot peak location to the given axis """
    # location of the BAO peak, label
    loc = np.array(a_sim_info.data["corr_func"]["par_peak"][0][cut])
    label = a_sim_info.app

    # comparison to the non-linear prediction
    loc_nl = np.array(a_sim_info.data["corr_func"]["nl_peak"][0][cut])
    ax.plot(a, loc / loc_nl, label=label + ' (loc)')

    # comparison to the inear prediction
    if not only_nl:
        loc_lin = np.array(a_sim_info.data["corr_func"]["lin_peak"][0][cut])
        ax.plot(a, loc / loc_lin, label=label + ' (loc, lin)')

def plot_peak_amp(a, a_sim_info, ax, cut, only_nl=True):
    """ plot peak amplitude to the given axis """
    # amplitude of the BAO peak
    amp = np.array(a_sim_info.data["corr_func"]["par_peak"][1][cut])
    label = a_sim_info.app

    # comparison to the non-linear prediction
    amp_nl = np.array(a_sim_info.data["corr_func"]["nl_peak"][1][cut])
    ax.plot(a, amp / amp_nl, label=label + ' (amp)')

    # comparison to the inear prediction
    if not only_nl:
        amp_lin = np.array(a_sim_info.data["corr_func"]["lin_peak"][1][cut])
        ax.plot(a, amp / amp_lin, label=label + ' (amp, lin)')

def plot_corr_peak(zs, sim_infos, out_dir='auto', save=True, show=False, use_z_eff=False):
    # output
    if out_dir == 'auto':
        out_dir = report_dir
    out_file = "corr_peak"

    # get rid of init zs
    cut = slice(1, None, None) if zs[0] == 'init' else slice(None, None, None)
    a = 1/(1.+np.array(zs[cut]))

    # figure
    fig = plt.figure(figsize=fig_size)
    ax = plt.gca()

    for a_sim_info in sim_infos:
        # peak location
        plot_peak_loc(a, a_sim_info, ax, cut)

        # peak amplitude
        plot_peak_amp(a, a_sim_info, ax, cut)    

    # labels
    plt.xlabel(r"$a$", fontsize=label_size)
    fig_suptitle(fig, "Relative BAO peak location and amplitude")

    # LEGEND manipulation
    legend_manipulation(ax, "", loc='lower left', bbox_to_anchor=(0.0, 0.0))
    plt.subplots_adjust(**subplt_adj_sym)

    # close & save figure
    close_fig(out_dir + out_file, fig, save=save, show=show, use_z_eff=use_z_eff)

def plot_eff_time_ax(a_sim_info, ax, a_eff_type="sigma_R"):
    # extract variables
    a = a_sim_info.data["eff_time"][a_eff_type]['a']
    D_eff_ratio = a_sim_info.data["eff_time"][a_eff_type]['D_eff_ratio']
    a_err = a_sim_info.data["eff_time"][a_eff_type]['a_err']

    # plot
    if a_eff_type == "sigma_R":
        label = a_sim_info.app +  '$: L = %i$ Mpc/h' % a_sim_info.box_opt["box_size"]
        ax.plot(a, D_eff_ratio, label=label)
    elif a_eff_type == "Pk":
        ax.errorbar(a, D_eff_ratio, yerr=a_err, label=a_sim_info.info_tr())
        ax.set_ylim(ymin=0.8)


def plot_eff_time(stack_infos, out_dir='auto', a_eff_type="sigma_R", save=True, show=False, use_z_eff=False):
    if out_dir == 'auto':
        out_dir = report_dir

    # figure
    fig = plt.figure(figsize=fig_size)
    ax = plt.gca()
    
    for stack_info in stack_infos:
        plot_eff_time_ax(stack_info, ax, a_eff_type)
    
    ax.set_ylabel(r'$D_{eff}/D_{GR}$', fontsize=label_size)
    ax.set_xlabel(r'$a$', fontsize=label_size)
    ax.legend()
    plt.subplots_adjust(**subplt_adj_sym)
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax.yaxis.grid(True)
    
    # save & show (in jupyter)
    close_fig(out_dir + 'D_eff_' + a_eff_type + '_', fig, save=save, show=show, use_z_eff=use_z_eff)


def plot_pwr_spec_diff_from_data(data_list, zs, a_sim_info, out_dir='auto', pk_type='dens', ext_title='par', save=True, show=False, use_z_eff=False):
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
            
    out_file += '_%s' % ext_title
    suptitle += ' (ref: %s)' % ext_title

    fig = plt.figure(figsize=fig_size)
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

    fig_suptitle(fig, suptitle)
    plt.xlabel(r"$k [h/$Mpc$]$", fontsize=label_size)
    plt.ylabel(r"$\frac{P(k)-P_{lin}(k)}{P_{lin}(k)}$", fontsize=25)
    # legend_manipulation(ax, a_sim_info.info_tr())
    legend_manipulation(ax, "")
    close_fig(out_dir + out_file, fig, save=save, show=show, use_z_eff=use_z_eff)


def plot_pwr_spec_diff_map_from_data(data_array, zs, a_sim_info, out_dir='auto', pk_type='dens', ext_title='', save=True, show=False, use_z_eff=False,
                                    vmin=-1, vmax=1):
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

    out_file += '_%s_map' % ext_title
    suptitle += ' (ref: %s)' % ext_title

    fig = plt.figure(figsize=(8, 8))
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

    fig_suptitle(fig, suptitle)
    ax.set_xlabel(r"$k [h/$Mpc$]$", fontsize=label_size)
    ax.set_ylabel(r"$a(t)$", fontsize=label_size)
    plt.draw()
    plt.figtext(0.5, 0.95, "",
                bbox={'facecolor': 'white', 'alpha': 0.2}, size=14, ha='center', va='top')
    plt.subplots_adjust(left=0.1, right=0.84, bottom=0.1, top=0.89)
    plt.subplots_adjust(left=0.1, right=0.9, bottom=0.1, top=0.9)

    close_fig(out_dir + out_file, fig, save=save, show=show, use_z_eff=use_z_eff)

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
    fig_suptitle(fig, "Power spectrum suppression" + suptitle)
    plt.xlabel(r"$a(t)$", fontsize=label_size)
    ylabel = r"$\langle{\frac{P(k)-P_{lin}(k)}{P_{lin}(k)}}\rangle$"
    if res is not None:
        ylabel += r', residual from $\nu=%.1f$' % res.nu
    plt.ylabel(ylabel, fontsize=25)
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
    ax.set_xlabel(r"$x [$Mpc$/h]$", fontsize=label_size)
    ax.set_ylabel(r"$z [$Mpc$/h]$", fontsize=label_size)
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
    ax.set_xlabel(r"$x [$Mpc$/h]$", fontsize=label_size)
    ax.set_ylabel(r"$z [$Mpc$/h]$", fontsize=label_size)
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
    ax.set_xlabel(r"$x [$Mpc$/h]$", fontsize=label_size)
    ax.set_ylabel(r"$z [$Mpc$/h]$", fontsize=label_size)
    L = int(np.sqrt(rho.shape[0]))
    rho.shape = L, L
    im = ax.imshow(rho, interpolation='bicubic', cmap='gnuplot',
                   norm=SymLogNorm(linthresh=1.0, linscale=1,
                                   vmin=-1, vmax=100), aspect='auto',
                   extent=[0, a_sim_info.box_opt["box_size"], 0, a_sim_info.box_opt["box_size"]])
    fig_suptitle(fig, "Slice through simulation box (overdensity), z = %.2f" % z)
    cbar = fig.colorbar(im, cax=cbar_ax, ticks=[-1, 0, 1, 10, 100])
    cbar.ax.set_yticklabels(['-1', '0', '1', '10', '> 100'])
    close_fig(out_dir + 'dens_z%.2f' % z, fig, save=save, show=show, use_z_eff=use_z_eff)

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
    ax.set_xlabel(r"$x [$Mpc$/h]$", fontsize=label_size)
    ax.set_ylabel(r"$z [$Mpc$/h]$", fontsize=label_size)

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


def plot_chi_evol(zs, a_sim_info, chi_opt=None, out_dir='auto', save=True, show=False, use_z_eff=False):
    """" Plot evolution of chameleon background values -- Compton wavelength and screening potential """
    if out_dir == 'auto':
        out_dir = a_sim_info.res_dir
    out_file = 'chi_evol'
    suptitle = "Evolution of Chameleon"
    fig = plt.figure(figsize=fig_size)
    cosmo = a_sim_info.sim.cosmo
    if chi_opt is None:
        chi_opt = [a_sim_info.chi_opt]
        
    ax1 = plt.subplot(311)
    ax2 = plt.subplot(313, sharex=ax1)
    ax3 = plt.subplot(312, sharex=ax1)
    
    ax1.set_yscale('log')
    ax2.set_yscale('log')
    ax3.set_yscale('log')
    
    zs = [z for z in zs if z != 'init']
    a = [1./(z+1) for z in zs]
    
    for chi in chi_opt:
        wavelengths = [power.chi_compton_wavelength(a_, cosmo, chi) for a_ in a]
        psi_a = [power.chi_psi_a(a_, chi) for a_ in a]
        chi_a = [power.chi_bulk_a(a_, chi, CHI_A_UNITS=False) for a_ in a]
        ax1.plot(zs, wavelengths, '-', label=r"$\Phi_{scr} = 10^{%i}$, $n=%.1f$" % (np.log10(chi["phi"]), chi["n"]))
        ax2.plot(zs, psi_a, '-')
        ax3.plot(zs, chi_a, '-')
    
    fig_suptitle(fig, suptitle)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax3.get_xticklabels(), visible=False)
    
    ax1.set_ylabel(r"$\lambda_C [$Mpc$/h]$", fontsize=label_size)
    ax2.set_ylabel(r"$\phi_{scr}$", fontsize=label_size)
    ax3.set_ylabel(r"$\chi/M_{pl}$", fontsize=label_size)
    ax2.set_xlabel(r"z", fontsize=label_size)
    
    # legend
    legend_manipulation(ax=ax1, loc='upper right')

    # subplots
    plt.subplots_adjust(hspace=0, **subplt_adj_sym)

    # close & save figure
    close_fig(out_dir + out_file, fig, save=save, show=show, use_z_eff=use_z_eff)


def plot_supp_lms(supp, a, a_sim_info, out_dir='auto', pk_type='dens', suptitle='', save=True, show=False, use_z_eff=False):
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

    fig_suptitle(fig, suptitle)
    plt.xlabel(r"$a(t)$", fontsize=label_size)
    plt.ylabel(
        r"$\langle{\frac{P(k)-P_{lin}(k)}{P_{lin}(k)}}\rangle$", fontsize=25)
    
    # legend
    # legend_manipulation(figtext=a_sim_info.info_tr())
    legend_manipulation(figtext="")

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
    ax.plot(x, y, 'ko', ms=3)