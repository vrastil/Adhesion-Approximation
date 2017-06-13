import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

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

def plot_pwr_spec(pwr_spec_files, zs, a_sim_info):
    fig = plt.figure(figsize=(12, 8))
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
    plt.figtext(x + 0.006, y - 0.02, a_sim_info.info(), bbox={'facecolor':'white', 'alpha':0.2}, size=14, ha='left', va='top')
    plt.show()

def plot_pwr_spec_diff(pwr_spec_diff_files, zs, a_sim_info):
    fig = plt.figure(figsize=(12, 8))
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
    plt.figtext(x + 0.006, y - 0.02, a_sim_info.info(), bbox={'facecolor':'white', 'alpha':0.2}, size=14, ha='left', va='top')
    plt.show()

