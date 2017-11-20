import math
import matplotlib
matplotlib.use('Agg')
from matplotlib import animation
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
# from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from scipy import interpolate
from scipy.integrate import quad
from scipy.misc import derivative

from . import get_files_in_traverse_dir


def trans_fce(k, Omega_m=1, h=0.67, q_log=2.34, q_1=3.89, q_2=16.2, q_3=5.47, q_4=6.71):
    if k == 0:
        return 1
    q = k / (Omega_m * h)
    tmp_log = np.log(1 + q_log * q) / (q_log * q)
    tmp_pol = 1 + q_1 * q + (q_2 * q)**2 + (q_3 * q)**3 + (q_4 * q)**4
    return tmp_log * tmp_pol**(-1. / 4)


def hubble_param(a, cosmo):
    Om = cosmo["Omega_b"] + cosmo["Omega_c"]
    OL = 1 - Om
    return np.sqrt(Om * abs(a)**(-3) + OL)


def growth_factor(a, cosmo=None):
    if cosmo is None:
        return a
    if a == 0:
        return 0

    def f(a_, cosmo_): return (a_ * hubble_param(a_, cosmo_))**(-3)
    D, err = quad(f, 0, a, args=(cosmo,))
    if "D_norm" not in cosmo:
        cosmo["D_norm"], err = quad(f, 0, 1, args=(cosmo,))
    return np.sign(a) * D * hubble_param(a, cosmo) / cosmo["D_norm"]


def growth_change(a, cosmo=None):
    if cosmo is None:
        return 1
    return derivative(growth_factor, a, dx=1e-6, args=(cosmo,))


def pwr_spec_prim(k, A=187826, n=1):
    return A * k**n


def pwr_spec(k, cosmo, ccl_cosmo=None, z=0, q_log=2.34, q_1=3.89, q_2=16.2, q_3=5.47, q_4=6.71):
    supp = np.exp(-k * k / cosmo["smoothing_k"]) if cosmo["smoothing_k"] else 1
    a = 1. / (z + 1.)
    supp *= growth_factor(a, cosmo)**2
    if ccl_cosmo is None:
        A, n, Omega_m, h = cosmo["A"], cosmo["index"], cosmo["Omega_m"], cosmo["h"]
        P_prim = pwr_spec_prim(k, A=A, n=n)
        T_k = trans_fce(k, Omega_m=Omega_m, h=h, q_log=q_log,
                        q_1=q_1, q_2=q_2, q_3=q_3, q_4=q_4)
        return supp * P_prim * T_k**2
    else:
        import pyccl as ccl
        return supp * ccl.linear_matter_power(ccl_cosmo, k * cosmo["h"], 1.) * cosmo["h"]**3


def plot_pwr_spec(pwr_spec_files, zs, a_sim_info, pwr_spec_files_extrap=None, pwr_spec_files_emu=None,
                  out_dir='auto', pk_type='dens', save=True, show=False):
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
    plt.yscale('log')
    plt.xscale('log')

    a_end = 1 / (1 + zs[-1])
    a_ = 0
    lab = 'init'
    init_idx_cor = 0  # if there is pwr_spec file 'init', extrap files have one file less
    for i, pwr in enumerate(pwr_spec_files):
        if zs[i] != 'init':
            a = 1 / (1 + zs[i])
            if ((a < 1.5 * a_) or (1.5 * a > a_end)) and a != a_end:
                continue
            a_ = a
            lab = 'z = ' + str(zs[i])
        else:
            init_idx_cor += 1

        data = np.loadtxt(pwr)
        k, P_k = data[:, 0], data[:, 1]
        plt.plot(k, P_k, 'o', ms=3, label=lab)
        del k, P_k, data

        if pwr_spec_files_extrap is not None and zs[i] != 'init':
            data = np.loadtxt(pwr_spec_files_extrap[i - init_idx_cor])
            k, P_k = data[:, 0], data[:, 1]
            plt.plot(k, P_k, 'k--')
            del k, P_k, data

    if pwr_spec_files_emu is not None:
        data = np.loadtxt(pwr_spec_files_emu[-1])
        k, P_k = data[:, 0], data[:, 1]
        plt.plot(k, P_k, '-')
        del k, P_k, data

    if a_sim_info.k_nyquist is not None:
        ls = [':', '-.', '--']
        ls *= (len(a_sim_info.k_nyquist) - 1) / 3 + 1
        ls = iter(ls)
        val_lab = {}
        for key, val in a_sim_info.k_nyquist.iteritems():
            if val in val_lab:
                val_lab[val] += ",\n" + " " * 8 + key
            else:
                val_lab[val] = r"$k_{Nq}$ (" + key
        for val, lab in val_lab.iteritems():
            ax.axvline(val, ls=next(ls), c='k', label=lab + r")")

    data = np.loadtxt(pwr_spec_files[-1])
    k = data[:, 0]
    k = np.logspace(np.log10(k[0]), np.log10(k[-1]), num=50)
    del data

    if pk_type == "dens":
        P_0 = [pwr_spec(k_, a_sim_info.cosmo,
                        ccl_cosmo=a_sim_info.ccl_cosmo, z=zs[-1]) for k_ in k]
        P_i = [pwr_spec(k_, a_sim_info.cosmo, ccl_cosmo=a_sim_info.ccl_cosmo,
                        z=zs[init_idx_cor]) for k_ in k]
        plt.plot(k, P_0, '-')
        plt.plot(k, P_i, '-')
    elif pk_type == "vel":
        P_0 = [pwr_spec(k_, a_sim_info.cosmo,
                        ccl_cosmo=a_sim_info.ccl_cosmo, z=0) for k_ in k]
        a_0 = 1. / (zs[init_idx_cor] + 1.)
        a_i = 1. / (zs[-1] + 1.)
        P_i = np.array(P_0) * growth_change(a_i, a_sim_info.cosmo)**2
        P_0 = np.array(P_0) * growth_change(a_0, a_sim_info.cosmo)**2
        plt.plot(k, P_i, '-')
        plt.plot(k, P_0, '-')

    fig.suptitle(suptitle, y=0.99, size=20)
    plt.xlabel(r"$k [h/$Mpc$]$", fontsize=15)
    plt.ylabel(r"$P(k) [$Mpc$/h)^3$]", fontsize=15)

    # LEGEND manipulation
    handles, labels = ax.get_legend_handles_labels()
    hl = sorted(zip(handles, labels), key=lambda x: x[1], reverse=True)
    handles, labels = zip(*hl)
    # pade_lab_index = labels.index("Pade app.")
    # labels.insert(-1, labels.pop(pade_lab_index))
    # handles.insert(-1, handles.pop(pade_lab_index))
    plt.legend(handles, labels, loc='upper left',
               bbox_to_anchor=(1.0, 1.0), fontsize=14)

    plt.draw()
    plt.figtext(0.5, 0.95, a_sim_info.info_tr(),
                bbox={'facecolor': 'white', 'alpha': 0.2}, size=14, ha='center', va='top')
    plt.subplots_adjust(left=0.1, right=0.84, bottom=0.1, top=0.89)

    if save:
        plt.savefig(out_dir + out_file)
    if show:
        plt.show()
    plt.close(fig)


def plot_pwr_spec_stacked(data_list, zs, a_sim_info, Pk_list_extrap=None, pwr_spec_files_emu=None,
                          out_dir='auto', pk_type='dens', save=True, show=False):
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
    plt.yscale('log')
    plt.xscale('log')

    a_end = 1 / (1 + zs[-1])
    a_ = 0
    lab = 'init'
    for i, data in enumerate(data_list):
        if zs[i] != 'init':
            a = 1 / (1 + zs[i])
            if ((a < 1.5 * a_) or (1.5 * a > a_end)) and a != a_end:
                continue
            a_ = a
            lab = 'z = ' + str(zs[i])

        k, P_k, P_k_std = data[0], data[1], data[2]
        # plt.errorbar(k, P_k, yerr=P_k_std, fmt='o', ms=3, label=lab)
        plt.plot(k, P_k, 'o', ms=3, label=lab)
        plt.fill_between(k, P_k - P_k_std, P_k + P_k_std,
                         facecolor='darkgrey', alpha=0.5)

        if Pk_list_extrap is not None:
            plt.plot(k, Pk_list_extrap[i].eval_list(k), 'k--')

        del k, P_k, P_k_std, data

    if pwr_spec_files_emu is not None:
        data = np.loadtxt(pwr_spec_files_emu[-1])
        k, P_k = data[:, 0], data[:, 1]
        plt.plot(k, P_k, '-')
        del k, P_k, data

    if a_sim_info.k_nyquist is not None:
        ls = [':', '-.', '--']
        ls *= (len(a_sim_info.k_nyquist) - 1) / 3 + 1
        ls = iter(ls)
        val_lab = {}
        for key, val in a_sim_info.k_nyquist.iteritems():
            if val in val_lab:
                val_lab[val] += ",\n" + " " * 8 + key
            else:
                val_lab[val] = r"$k_{Nq}$ (" + key
        for val, lab in val_lab.iteritems():
            ax.axvline(val, ls=next(ls), c='k', label=lab + r")")

    data = data_list[-1]
    k = data[0]
    k = np.logspace(np.log10(k[0]), np.log10(k[-1]), num=50)
    del data

    if pk_type == "dens":
        z0 = zs[0] if zs[0] != "init" else zs[1]
        P_0 = [pwr_spec(k_, a_sim_info.cosmo,
                        ccl_cosmo=a_sim_info.ccl_cosmo, z=zs[-1]) for k_ in k]
        P_i = [pwr_spec(k_, a_sim_info.cosmo,
                        ccl_cosmo=a_sim_info.ccl_cosmo, z=z0) for k_ in k]
        plt.plot(k, P_0, '-')
        plt.plot(k, P_i, '-')
    elif pk_type == "vel":
        P_0 = [pwr_spec(k_, a_sim_info.cosmo,
                        ccl_cosmo=a_sim_info.ccl_cosmo, z=0) for k_ in k]
        a_0 = 1. / (zs[0] + 1.)
        a_i = 1. / (zs[-1] + 1.)
        P_i = np.array(P_0) * growth_change(a_i, a_sim_info.cosmo)**2
        P_0 = np.array(P_0) * growth_change(a_0, a_sim_info.cosmo)**2
        plt.plot(k, P_i, '-')
        plt.plot(k, P_0, '-')

    fig.suptitle(suptitle, y=0.99, size=20)
    plt.xlabel(r"$k [h/$Mpc$]$", fontsize=15)
    plt.ylabel(r"$P(k) [$Mpc$/h)^3$]", fontsize=15)

    # LEGEND manipulation
    handles, labels = ax.get_legend_handles_labels()
    hl = sorted(zip(handles, labels), key=lambda x: x[1], reverse=True)
    handles, labels = zip(*hl)
    # pade_lab_index = labels.index("Pade app.")
    # labels.insert(-1, labels.pop(pade_lab_index))
    # handles.insert(-1, handles.pop(pade_lab_index))
    plt.legend(handles, labels, loc='upper left',
               bbox_to_anchor=(1.0, 1.0), fontsize=14)

    plt.draw()
    plt.figtext(0.5, 0.95, a_sim_info.info_tr(),
                bbox={'facecolor': 'white', 'alpha': 0.2}, size=14, ha='center', va='top')
    plt.subplots_adjust(left=0.1, right=0.84, bottom=0.1, top=0.89)

    if save:
        plt.savefig(out_dir + out_file)
    if show:
        plt.show()
    plt.close(fig)


def plot_corr_func_from_data_single(corr_data, z, a_sim_info, corr_data_lin=None, corr_data_emu=None, out_dir='auto', save=True, show=False):
    if out_dir == 'auto':
        out_dir = a_sim_info.res_dir

    # *****************
    # first plot, xi(r)
    # *****************
    fig = plt.figure(figsize=(15, 11))
    if z == 'init':
        lab = z
        z_out = z
    else:
        lab = 'z = ' + str(z)
        z_out = 'z%.2f' % z
    # load all data, plot xi(r)
    r, xi = corr_data[0], corr_data[1]
    plt.plot(r, xi, 'o', ms=3, label=lab)
    if corr_data_lin is not None:
        r_lin, xi_lin = corr_data_lin[0], corr_data_lin[1]
        plt.plot(r_lin, xi_lin, '-', label=r"$\Lambda$CDM (lin)")
    if corr_data_emu is not None:
        r_emu, xi_emu = corr_data_emu[0], corr_data_emu[1]
        plt.plot(r_emu, xi_emu, '-', label=r"$\Lambda$CDM (emu)")
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
    if save:
        plt.savefig(out_dir + 'corr_func_%s.png' % z_out)
    if show:
        plt.show()
    plt.close(fig)

    # **********************
    # second plot, r*r*xi(r)
    # **********************
    fig = plt.figure(figsize=(15, 11))
    # plot r*r*xi(r)
    plt.plot(r, r * r * xi, 'o', ms=3, label=lab)
    if corr_data_lin is not None:
        plt.plot(r_lin, r_lin * r_lin * xi_lin,
                 '-', label=r"$\Lambda$CDM (lin)")
    if corr_data_emu is not None:
        plt.plot(r_emu, r_emu * r_emu * xi_emu,
                 '-', label=r"$\Lambda$CDM (emu)")
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
    if save:
        plt.savefig(out_dir + 'corr_r2_func_%s.png' % z_out)
    if show:
        plt.show()
    plt.close(fig)


# correlation function stacked data, linear and emu corr. func in files
def plot_corr_func_from_data(corr_data_list, zs, a_sim_info, corr_func_files_lin=None, corr_func_files_emu=None, zs_emu=None, out_dir='auto', save=True, show=False):

    a_end = 1 / (1 + zs[-1])
    a_ = 0
    j = 0
    init_idx_cor = 0  # if there is pwr_spec file 'init', extrap files have one file less
    for i, z in enumerate(zs):
        # input correlation function handled separately
        if z == "init":
            init_idx_cor += 1
            continue
        # load emulator correlation BEFORE disregarding the plot (near times),
        # still need to update j in this case
        if corr_func_files_emu is not None and z == zs_emu[j]:
            corr_data_emu = np.transpose(np.loadtxt(corr_func_files_emu[j]))
            j += 1
        else:
            corr_data_emu = None

        a = 1 / (1 + z)
        # no need to plot EVERY correlation function
        if a < 2 * a_ and z > 1 and a != a_end:
            continue
        a_ = a
        corr_data = corr_data_list[i]
        if corr_func_files_lin is not None:
            corr_data_lin = np.transpose(np.loadtxt(
                corr_func_files_lin[i - init_idx_cor]))
        else:
            corr_data_lin = None

        plot_corr_func_from_data_single(
            corr_data, z, a_sim_info, corr_data_lin, corr_data_emu, out_dir, save, show)
        del corr_data, corr_data_lin, corr_data_emu
    # input correlation function
    if zs[0] == "init" and zs_emu[-1] == 0:
        corr_data_emu = np.transpose(np.loadtxt(corr_func_files_emu[-1]))
        corr_data_lin = np.transpose(np.loadtxt(corr_func_files_lin[-1]))
        corr_data = corr_data_list[0]
        plot_corr_func_from_data_single(
            corr_data, zs[0], a_sim_info, corr_data_lin, corr_data_emu, out_dir, save, show)


# correlation function directly from siulation, no stacking, all data in files
def plot_corr_func(corr_func_files, zs, a_sim_info, corr_func_files_lin=None, corr_func_files_emu=None, zs_emu=None, out_dir='auto', save=True, show=False):
    corr_data_list = [np.transpose(np.loadtxt(
        corr_func_files[i])) for i in xrange(len(zs))]
    plot_corr_func_from_data(corr_data_list, zs, a_sim_info, corr_func_files_lin,
                             corr_func_files_emu, zs_emu, out_dir, save, show)


def plot_pwr_spec_diff_from_data(data_list, zs, a_sim_info, out_dir='auto', pk_type='dens', ext_title='', save=True, show=False):
    if out_dir == 'auto':
        out_dir = a_sim_info.res_dir
    if pk_type == "dens":
        out_file = 'pwr_spec_diff'
        suptitle = "Power spectrum difference"
    elif pk_type == "vel":
        out_file = 'vel_pwr_spec_diff'
        suptitle = r"Power spectrum difference $(\nabla\cdot u)$"

    if ext_title == 'input':
        ext_title = ' (ref: input)'
        out_file += '_input'
    if ext_title == 'hybrid':
        ext_title = ' (ref: hybrid)'
        out_file += '_hybrid'
    if ext_title == 'par':
        ext_title = ' (ref: particle)'
        out_file += '_par'
    out_file += '.png'
    suptitle += ext_title

    fig = plt.figure(figsize=(15, 11))
    ax = plt.gca()
    plt.xscale('log')

    # SMALL / MEDIUM / LARGE SCALE VALUES
    k = data_list[-1][0]
    # half of nyquist wavelength, 7 subintervals
    idx = (np.abs(k - 0.5 * a_sim_info.k_nyquist["particle"])).argmin() / 7
    cmap = cm.get_cmap('gnuplot')
    ax.axvspan(k[0 * idx], k[1 * idx], alpha=0.2, color=cmap(0.1))
    ax.axvspan(k[3 * idx], k[4 * idx], alpha=0.3, color=cmap(0.5))
    ax.axvspan(k[6 * idx], k[7 * idx], alpha=0.4, color=cmap(0.9))
    del k

    a_end = 1 / (1 + zs[-1])
    a_ = 0
    ymin = ymax = 0
    for i, data in enumerate(data_list):
        a = 1 / (1 + zs[i])
        if ((a < 1.5 * a_) or (1.5 * a > a_end)) and a != a_end:
            continue
        a_ = a
        lab = 'z = ' + str(zs[i])
        k, P_k = data[0], data[1]
        try:
            P_k_std = data[2]
        except IndexError:
            plt.plot(k, P_k, 'o', ms=3, label=lab)
        else:
            plt.errorbar(k, P_k, fmt='o', yerr=P_k_std, ms=3, label=lab)
            del P_k_std

        ymax = max(ymax, np.max(P_k[0:7 * idx]))
        ymin = min(ymin, np.min(P_k[0:7 * idx]))
        del k, P_k

    if a_sim_info.k_nyquist is not None:
        ls = [':', '-.', '--']
        ls *= (len(a_sim_info.k_nyquist) - 1) / 3 + 1
        ls = iter(ls)
        val_lab = {}
        for key, val in a_sim_info.k_nyquist.iteritems():
            if val in val_lab:
                val_lab[val] += ",\n" + " " * 8 + key
            else:
                val_lab[val] = r"$k_{Nq}$ (" + key
        for val, lab in val_lab.iteritems():
            ax.axvline(val, ls=next(ls), c='k', label=lab + r")")

    if ymax > 1:
        ymax = 1
    ymax = math.ceil(ymax / 0.1) * 0.1
    ymin = math.floor(ymin / 0.1) * 0.1
    ax.yaxis.grid(True)
    plt.ylim(ymin=ymin, ymax=ymax)

    fig.suptitle(suptitle, y=0.99, size=20)
    plt.xlabel(r"$k [h/$Mpc$]$", fontsize=15)
    plt.ylabel(r"$\frac{P(k)-P_{lin}(k)}{P_{lin}(k)}$", fontsize=25)
    plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0), fontsize=14)
    plt.draw()
    plt.figtext(0.5, 0.95, a_sim_info.info_tr(),
                bbox={'facecolor': 'white', 'alpha': 0.2}, size=14, ha='center', va='top')
    plt.subplots_adjust(left=0.1, right=0.84, bottom=0.1, top=0.89)
    if save:
        plt.savefig(out_dir + out_file)
    if show:
        plt.show()
    plt.close(fig)


def plot_pwr_spec_diff_map_from_data(data_list, zs, a_sim_info, out_dir='auto', pk_type='dens', ext_title='', save=True, show=False):
    if out_dir == 'auto':
        out_dir = a_sim_info.res_dir
    if pk_type == "dens":
        out_file = 'pwr_spec_diff'
        suptitle = "Power spectrum difference"
    elif pk_type == "vel":
        out_file = 'vel_pwr_spec_diff'
        suptitle = r"Power spectrum difference $(\nabla\cdot u)$"

    if ext_title == 'input':
        ext_title = ' (ref: input)'
        out_file += '_input'
    if ext_title == 'hybrid':
        ext_title = ' (ref: hybrid)'
        out_file += '_hybrid'
    if ext_title == 'par':
        ext_title = ' (ref: particle)'
        out_file += '_par'
    out_file += '_map.png'
    suptitle += ext_title

    fig = plt.figure(figsize=(10, 10))
    gs = gridspec.GridSpec(1, 15, wspace=0.5)
    ax = plt.subplot(gs[0, : -1])
    cbar_ax = plt.subplot(gs[0, -1])
    ax.set_xscale('log')

    a = [1 / (1 + z) for z in zs]
    supp = np.array(data_list)[:,1,:] # extract Pk, shape = (zs, k)
    k = np.array(data_list[0][0])

    im = ax.pcolor(k, a, supp, cmap='seismic', vmin=-0.2, vmax=0.2)
    cbar = fig.colorbar(im, cax=cbar_ax)

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
    if save:
        plt.savefig(out_dir + out_file)
    if show:
        plt.show()
    plt.close(fig)


def plot_pwr_spec_diff(pwr_spec_diff_files, zs, a_sim_info, out_dir='auto', pk_type='dens', ext_title='', save=True, show=False):
    data_list = [np.transpose(np.loadtxt(a_file))
                 for a_file in pwr_spec_diff_files]
    plot_pwr_spec_diff_from_data(
        data_list, zs, a_sim_info, out_dir, pk_type, ext_title, save, show)


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
    if save:
        plt.savefig(out_dir + 'supp.png')
    if show:
        plt.show()
    plt.close(fig)


def plot_dens_histo(dens_bin_files, zs, a_sim_info, out_dir='auto', fix_N=1, fix_rho=0.0, save=True, show=False):
    if out_dir == 'auto':
        out_dir = a_sim_info.res_dir
    num_sub_x = 3
    num_sub_y = (len(zs) + num_sub_x - 1) / num_sub_x
    fig = plt.figure(figsize=(num_sub_x * 5, num_sub_y * 5.5))

    gs = gridspec.GridSpec(num_sub_y, num_sub_x, wspace=0.2,
                           hspace=0.3, left=0.1, right=0.84, bottom=0.1, top=0.89)

    for i, rho_fl in enumerate(dens_bin_files):
        lab = 'z = ' + str(zs[i])
        data = np.loadtxt(rho_fl)
        rho, count = data[:, 0], data[:, 1]
        count *= fix_N
        rho += fix_rho
        xmin = -1
        xmax = rho[np.nonzero(count)[0][-1]] + 1
        ax = plt.subplot(gs[i])
        # plt.subplot(num_sub_y, num_sub_x, i + 1)
        ax.set_xlim(xmin=xmin, xmax=xmax)
        ax.hist(rho, bins=20, weights=count, facecolor='green',
                edgecolor='black', linewidth=0.8, normed=True)
        ax.set_yscale('log', nonposy='clip')
        ax.set_title(lab)
        del rho, count, data

    fig.suptitle("Overdensity distribution", y=0.99, size=20)

    plt.figtext(0.5, 0.95, a_sim_info.info_tr(),
                bbox={'facecolor': 'white', 'alpha': 0.2}, size=14, ha='center', va='top')
    # plt.subplots_adjust(wspace=0.2, hspace=0.3)
    # plt.subplots_adjust(wspace=0.2, hspace=0.3, left=0.1, right=0.84, top=0.9)
    if save:
        plt.savefig(out_dir + 'dens_histo.png')
    if show:
        plt.show()
    plt.close(fig)


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
    if save:
        plt.savefig(out_dir + 'par_evol_last.png')
    if show:
        plt.show()
    plt.close(fig)


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
    plt.close(fig)
    del ani


def plot_dens_one_slice(rho, z, a_sim_info, out_dir='auto', save=True, show=False):
    if out_dir == 'auto':
        out_dir = a_sim_info.res_dir
    from matplotlib.colors import SymLogNorm
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
    if save:
        plt.savefig(out_dir + 'dens_z%.2f.png' % z)
    if show:
        plt.show()
    plt.close(fig)


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

    from matplotlib.colors import SymLogNorm
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
    plt.close(fig)
    del ani


def plot_supp_lms(supp_lms, a, a_sim_info, k_lms, supp_std_lms, out_dir='auto', pk_type='dens', suptitle='', save=True, show=False):
    if out_dir == 'auto':
        out_dir = a_sim_info.res_dir
    if pk_type == "dens":
        out_file = 'supp.png'
        suptitle = "Power spectrum suppression"
    elif pk_type == "vel":
        out_file = 'supp_vel.png'
        suptitle = r"Power spectrum suppression $(\nabla\cdot u)$"
    fig = plt.figure(figsize=(15, 11))
    supp_l, supp_m, supp_s = supp_lms
    supp_large_std, supp_medium_std, supp_small_std = supp_std_lms
    k_l, k_m, k_s = k_lms

    cmap = cm.get_cmap('gnuplot')
    plt.errorbar(a, supp_l, fmt='-o', yerr=supp_large_std, ms=3, color=cmap(0.1), lw=4,
                 label='Large-scale:\n' r'$\langle%.2f,%.2f\rangle$' % (k_l[0], k_l[1]))
    plt.errorbar(a, supp_m, fmt='-o', yerr=supp_medium_std, ms=3, color=cmap(0.5), lw=2.5,
                 label='Medium-scale:\n'r'$\langle%.2f,%.2f\rangle$' % (k_m[0], k_m[1]))
    plt.errorbar(a, supp_s, fmt='-o', yerr=supp_small_std, ms=3, color=cmap(0.9), lw=1,
                 label='Small-scale:\n'r'$\langle%.2f,%.2f\rangle$' % (k_s[0], k_s[1]))
    # plt.plot(a, supp_l, '-o', ms=3, label='Large-scale:\n' r'$\langle%.2f,%.2f\rangle$ h/Mpc' % (k_l[0], k_l[1]))
    # plt.plot(a, supp_m, '-o', ms=3, label='Medium-scale:\n'r'$\langle%.2f,%.2f\rangle$ h/Mpc' % (k_m[0], k_m[1]))
    # plt.plot(a, supp_s, '-o', ms=3, label='Small-scale:\n'r'$\langle%.2f,%.2f\rangle$ h/Mpc' % (k_s[0], k_s[1]))

    ax = plt.gca()
    ax.yaxis.grid(True)
    ymin, ymax = ax.get_ylim()
    if ymax > 1:
        ymax = 1
    ymax = math.ceil(ymax / 0.1) * 0.1
    ymin = math.floor(ymin / 0.1) * 0.1
    plt.ylim(ymin=ymin, ymax=ymax)

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
    if save:
        plt.savefig(out_dir + out_file)
    if show:
        plt.show()
    plt.close(fig)


def plot_all_single_supp(res, out_dir='/home/vrastil/Documents/GIT/Adhesion-Approximation/output/supp_comparison/',
                         Nm=0, Np=0, L=0, nu=0, rs=0, app=''):
    subfiles = res.get_subfiles(Nm=Nm, Np=Np, L=L, nu=nu, rs=rs, app=app)
    for a_sim_info in subfiles:
        res.load_k_supp(a_sim_info)
        plot_supp_lms(a_sim_info.supp[0], a_sim_info.a, a_sim_info,
                      k_lms=a_sim_info.supp[1], show=True, save=True)
