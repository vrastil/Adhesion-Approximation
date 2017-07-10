import matplotlib
matplotlib.use('Agg')
from matplotlib import animation
import matplotlib.pyplot as plt
import numpy as np

from . import get_files_in_traverse_dir

def trans_fce(k, Omega_0=1, h=0.67, q_log=2.34, q_1=3.89, q_2=16.1, q_3=5.4, q_4=6.71):
    if k == 0: return 1
    q = k/(Omega_0*h)
    tmp_log = np.log(1 + q_log*q)/(q_log*q)
    tmp_pol = 1 + q_1*q + (q_2*q)**2 + (q_3*q)**3 + (q_4*q)**4
    return tmp_log*tmp_pol**(-1./4)

def pwr_spec_prim(k, A=187826, n=1):
    return A*k**n

def pwr_spec(k, z=0, A=187826, n=1, Omega_0=1, h=0.67, q_log=2.34, q_1=3.89, q_2=16.1, q_3=5.4, q_4=6.71):
    P_prim = pwr_spec_prim(k, A=A, n=n)
    T_k = trans_fce(k, Omega_0=Omega_0, h=h, q_log=q_log, q_1=q_1, q_2=q_2, q_3=q_3, q_4=q_4)
    return P_prim*T_k**2 / (z+1)**2

def plot_pwr_spec(pwr_spec_files, zs, a_sim_info, out_dir, save=True, show=False):
    fig = plt.figure(figsize=(14, 8))
    plt.yscale('log')
    plt.xscale('log')
    a_ = 0

    lab = 'init'
    for i, pwr in enumerate(pwr_spec_files):
        if zs[i] != 'init':
            a = 1/(1+zs[i])
            if (a < 2.4*a_) and a != 1: continue
            a_ = a
            lab = 'z = ' + str(zs[i])
        data = np.loadtxt(pwr)
        k, P_k = data[:, 0], data[:, 1]
        plt.plot(k, P_k, 'o', ms=3, label=lab)
    k = np.logspace(np.log10(k[0]), np.log10(k[-1]), num=20)
    P_0 = [pwr_spec(k_) for k_ in k]
    P_i = [pwr_spec(k_, z=200) for k_ in k]
    plt.plot(k, P_0, '-')
    plt.plot(k, P_i, '-')

    fig.suptitle("Power spectrum", y=0.95, size=20)
    plt.xlabel(r"$k [h/$Mpc$]$", fontsize=15)
    plt.ylabel(r"$P(k) [$Mpc$/h)^3$]", fontsize=15)
    leg = plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0), fontsize=14)
    plt.draw()
    box = leg.get_frame().get_bbox()
    x = box.bounds[0]/fig.bbox.bounds[2]
    y = box.bounds[1]/fig.bbox.bounds[3]
    plt.figtext(x - 0.045, y - 0.01, a_sim_info.info(),
                bbox={'facecolor':'white', 'alpha':0.2}, size=14, ha='left', va='top')
    plt.subplots_adjust(left=0.1, right=0.85, top=0.9, bottom=0.1)
    if show: plt.show()
    if save: plt.savefig(out_dir + 'pwr_spec.png')
    plt.close(fig)

def plot_pwr_spec_diff(pwr_spec_diff_files, zs, a_sim_info, out_dir, save=True, show=False):
    fig = plt.figure(figsize=(14, 8))
    plt.xscale('log')
    a_ = 0

    lab = 'init'
    for i, pwr in enumerate(pwr_spec_diff_files):
        if zs[i] != 'init':
            a = 1/(1+zs[i])
            if (a < 2.4*a_) and a != 1: continue
            a_ = a
            lab = 'z = ' + str(zs[i])
        data = np.loadtxt(pwr)
        k, P_k = data[:, 0], data[:, 1]
        plt.plot(k, P_k, 'o', ms=3, label=lab)

    ymax = 0.2
    ymin = np.min(P_k)
    for y in np.arange(ymax-0.2, ymin, -0.2):
        plt.axhline(y=y, color='black', lw=0.2)
    plt.ylim(ymin=ymin, ymax=ymax)
    fig.suptitle("Power spectrum difference", y=0.95, size=20)
    plt.xlabel(r"$k [h/$Mpc$]$", fontsize=15)
    plt.ylabel(r"$\frac{P(k)-P_{lin}(k)}{P_{lin}(k)}$", fontsize=25)
    leg = plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0), fontsize=14)
    plt.draw()
    box = leg.get_frame().get_bbox()
    x = box.bounds[0]/fig.bbox.bounds[2]
    y = box.bounds[1]/fig.bbox.bounds[3]
    plt.figtext(x - 0.045, y - 0.01, a_sim_info.info(),
                bbox={'facecolor':'white', 'alpha':0.2}, size=14, ha='left', va='top')
    plt.subplots_adjust(left=0.1, right=0.85, top=0.9, bottom=0.1)
    if show: plt.show()
    if save: plt.savefig(out_dir + 'pwr_spec_diff.png')
    plt.close(fig)

def plot_supp(sim_infos, out_dir, suptitle='', save=True, show=False):
    fig = plt.figure(figsize=(14, 8))
    for a_sim_info in sim_infos:
        try:
            supp_fl = get_files_in_traverse_dir(a_sim_info.dir + 'supp/', '*.dat')[0][0]
        except ValueError:
            print "WARNING! Missing data in '%s'. Skipping step." % (a_sim_info.dir + 'supp/')
        else:
            data = np.loadtxt(supp_fl)
            a, supp = data[:, 0], data[:, 1]
            plt.plot(a, supp, '-o', ms=3, label=a_sim_info.info_supp())

    #plt.ylim(ymin=-1, ymax=0)
    fig.suptitle("Power spectrum suppression" + suptitle, y=0.95, size=20)
    plt.xlabel(r"$a(t)$", fontsize=15)
    plt.ylabel(r"$\langle{\frac{P(k)-P_{lin}(k)}{P_{lin}(k)}}\rangle$", fontsize=25)
    plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0), fontsize=14)
    plt.subplots_adjust(left=0.1, right=0.7, top=0.9, bottom=0.1)
    if show: plt.show()
    if save: plt.savefig(out_dir + 'supp.png')
    plt.close(fig)

def plot_dens_histo(dens_bin_files, zs, a_sim_info, out_dir, fix_N=1, fix_rho=0.0, save=True, show=False):
    num_sub_x = 3
    num_sub_y = (len(zs) + num_sub_x - 1) / num_sub_x
    fig = plt.figure(figsize=(num_sub_x*5, num_sub_y*6))

    for i, rho_fl in enumerate(dens_bin_files):
        lab = 'z = ' + str(zs[i])
        data = np.loadtxt(rho_fl)
        rho, count = data[:, 0], data[:, 1]
        count *= fix_N
        rho += fix_rho
        xmin = -1
        xmax = rho[np.nonzero(count)[0][-1]] +1
        plt.subplot(num_sub_y, num_sub_x, i+1)
        plt.xlim(xmin=xmin, xmax=xmax)
        plt.hist(rho, bins=20, weights=count, facecolor='green', edgecolor='black', linewidth=0.8)
        plt.yscale('log', nonposy='clip')
        plt.title(lab)

    fig.suptitle("Overdensity distribution", y=0.97, size=20)

    plt.figtext(0.5, 0.94, a_sim_info.info_tr(),
                bbox={'facecolor':'white', 'alpha':0.2}, size=14, ha='center', va='top')
    plt.subplots_adjust(wspace=0.2, hspace=0.3)
    if show: plt.show()
    if save: plt.savefig(out_dir + 'dens_histo.png')
    plt.close(fig)

def plot_par_last_slice(files, files_t, zs, a_sim_info, out_dir, save=True, show=False):
    fig = plt.figure(figsize=(10, 10))
    data = np.loadtxt(files[0])
    x, y = data[:, 0], data[:, 1]
    ax = plt.axes(xlim=(0, np.max(x)), ylim=(0, np.max(y)))
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
        ax.plot(x_t[i], y_t[i], '--or', ms=4, lw=1.5, markevery=(num_steps-1, num_steps))

    plt.figtext(0.5, 0.94, a_sim_info.info_tr(),
                bbox={'facecolor':'white', 'alpha':0.2}, size=14, ha='center', va='top')
    plt.xlabel(r"$x [$Mpc$/h]$", fontsize=13)
    plt.ylabel(r"$z [$Mpc$/h]$", fontsize=13)
    plt.title("Slice through simulation box (particles), z = %.2f" % zs[-1], y=1.09, size=20)
    if show: plt.show()
    if save: plt.savefig(out_dir + 'par_evol_last.png')
    plt.close(fig)

def plot_par_evol(files, files_t, zs, a_sim_info, out_dir, save=True):
    data = np.loadtxt(files[0])
    x, y = data[:, 0], data[:, 1]
    data = np.loadtxt(files_t[0])

    num = len(zs)
    num_track = len(data)

    fig = plt.figure(figsize=(10, 10))
    ax = plt.axes(xlim=(0, np.max(x)), ylim=(0, np.max(y)))
    line, = ax.plot([], [], 'ob', ms=1, animated=True)
    lines_t = []
    for i in xrange(num_track):
        lines_t.append(ax.plot([], [], '--or', ms=4, lw=1.5, animated=True)[0])

    plt.figtext(0.5, 0.94, a_sim_info.info_tr(),
                bbox={'facecolor':'white', 'alpha':0.2}, size=14, ha='center', va='top')
    plt.xlabel(r"$x [$Mpc$/h]$", fontsize=13)
    plt.ylabel(r"$z [$Mpc$/h]$", fontsize=13)
    del x, y, data

    def animate(j):
        if j < num: i = j
        else: i = 2*num - j - 1
        plt.title("Slice through simulation box (particles), z = %.2f" % zs[i], y=1.09, size=20)
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
            line_t.set_markevery((num_steps-1, num_steps))
        del x, y, x_t, y_t, data
        return [line] + lines_t

    ani = animation.FuncAnimation(fig, animate, frames=2*num, interval=250, blit=True)
    if save: ani.save(out_dir + 'par_evol.gif', writer='imagemagick')
    plt.close(fig)
    del ani

def plot_dens_one_slice(rho, z, out_dir, a_sim_info, save=True, show=False):
    from matplotlib.colors import SymLogNorm
    fig, ax = plt.subplots(figsize=(10, 10))
    plt.figtext(0.5, 0.9, a_sim_info.info_tr(),
                bbox={'facecolor':'white', 'alpha':0.2}, size=14, ha='center', va='top')
    plt.xlabel(r"$x [$Mpc$/h]$", fontsize=13)
    plt.ylabel(r"$z [$Mpc$/h]$", fontsize=13)
    cbar_ax = fig.add_axes([0.85, 0.155, 0.05, 0.695])
    fig.subplots_adjust(right=0.82)
    L = int(np.sqrt(rho.shape[0]))
    rho.shape = L, L
    im = ax.imshow(rho, interpolation='bicubic', cmap='gnuplot',
                   norm=SymLogNorm(linthresh=1.0, linscale=1, vmin=-1, vmax=100),
                   extent=[0, a_sim_info.box, 0, a_sim_info.box])
    fig.suptitle("Slice through simulation box (overdensity), z = %.2f" % z, y=0.95, size=20)
    fig.colorbar(im, cax=cbar_ax)
    if show: plt.show()
    if save: plt.savefig(out_dir + 'dens_z%.2f.png' % z)
    plt.close(fig)

def plot_dens_two_slices(files, zs, a_sim_info, out_dir, save=True, show=False):
    half = len(files) / 2
    rho, z = np.loadtxt(files[half])[:, 2], zs[half]
    plot_dens_one_slice(rho, z, out_dir, a_sim_info, save=save, show=show)
    rho, z = np.loadtxt(files[-1])[:, 2], zs[-1]
    plot_dens_one_slice(rho, z, out_dir, a_sim_info, save=save, show=show)

def plot_dens_evol(files, zs, a_sim_info, out_dir, save=True):
    from matplotlib.colors import SymLogNorm
    num = len(zs)

    fig, ax = plt.subplots(figsize=(10, 10))
    plt.figtext(0.5, 0.9, a_sim_info.info_tr(),
                bbox={'facecolor':'white', 'alpha':0.2}, size=14, ha='center', va='top')
    plt.xlabel(r"$x [$Mpc$/h]$", fontsize=13)
    plt.ylabel(r"$z [$Mpc$/h]$", fontsize=13)
    cbar_ax = fig.add_axes([0.85, 0.155, 0.05, 0.695])
    fig.subplots_adjust(right=0.82)

    def animate(j):
        if j < num: i = j
        else: i = 2*num - j - 1
        rho = np.loadtxt(files[i])[:, 2]
        L = int(np.sqrt(rho.shape[0]))
        rho.shape = L, L
        im = ax.imshow(rho, interpolation='bicubic', cmap='gnuplot', animated=True,
                       norm=SymLogNorm(linthresh=1.0, linscale=1, vmin=-1, vmax=100),
                       extent=[0, a_sim_info.box, 0, a_sim_info.box])
        fig.suptitle("Slice through simulation box (overdensity), z = %.2f" % zs[i], y=0.95, size=20)
        fig.colorbar(im, cax=cbar_ax)
        del rho
        return [im]

    ani = animation.FuncAnimation(fig, animate, frames=2*num, interval=250, blit=True)
    if save: ani.save(out_dir + 'dens_evol.gif', writer='imagemagick')
    plt.close(fig)
    del ani

def plot_supp_lms(supp_lms, a, a_sim_info, out_dir, k_lms=None, suptitle='', save=True, show=False):
    fig = plt.figure(figsize=(14, 8))
    supp_l, supp_m, supp_s = supp_lms
    k_l, k_m, k_s = k_lms
    plt.plot(a, supp_l, '-o', ms=3, label=r'Large-scale: $\langle%.2f,%.2f\rangle$ h/Mpc' % (k_l[0], k_l[1]))
    plt.plot(a, supp_m, '-o', ms=3, label=r'Medium-scale: $\langle%.2f,%.2f\rangle$ h/Mpc' % (k_m[0], k_m[1]))
    plt.plot(a, supp_s, '-o', ms=3, label=r'Small-scale: $\langle%.2f,%.2f\rangle$ h/Mpc' % (k_s[0], k_s[1]))

    fig.suptitle("Power spectrum suppression", y=0.95, size=20)
    plt.xlabel(r"$a(t)$", fontsize=15)
    plt.ylabel(r"$\langle{\frac{P(k)-P_{lin}(k)}{P_{lin}(k)}}\rangle$", fontsize=25)
    leg = plt.legend(loc='upper left', bbox_to_anchor=(1.0, 1.0), fontsize=14)
    plt.subplots_adjust(left=0.1, right=0.7, top=0.9, bottom=0.1)
    plt.draw()
    box = leg.get_frame().get_bbox()
    x = box.bounds[0]/fig.bbox.bounds[2]
    y = box.bounds[1]/fig.bbox.bounds[3]
    plt.figtext(x + 0.005, y - 0.01, a_sim_info.info(),
                bbox={'facecolor':'white', 'alpha':0.2}, size=14, ha='left', va='top')
    if show: plt.show()
    if save: plt.savefig(out_dir + 'supp.png')
    plt.close(fig)

