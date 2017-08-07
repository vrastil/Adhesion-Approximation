from datetime import datetime
import sys
import gc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import json

from . import *
from . import plot

class SimInfo(object):
    def __init__(self, *args):
        self.results = False
        if len(args) == 2:
            self.load_file(*args)
        else:
            self.num_g = 0
            self.num_p = 0
            self.num_m = 0
            self.box = 0
            self.nu = 0
            self.rs = 0
            self.app = ''
            self.dir = ''
            self.res_dir = ''

    def info(self):
        info = ''
        info += '%s:\n' % self.app
        info += '$N_m = %i$\n' % self.num_m
        info += '$N_p = %i^3$\n' % self.num_p
        info += '$L = %i$ Mpc/h\n' % self.box
        if self.app == 'AA': info += r'$\nu = %.1f$' % self.nu
        if self.app == 'FP_pp': info += r'$r_s = %.1f$' % self.rs
        return info

    def info_tr(self):
        info = self.info()
        return info.replace('\n', '  ')

    def info_supp(self):
        info = '%s: $res = %.2f$' % (self.app, float(self.num_m) / self.box)
        if self.app == 'AA': info += r', $\nu = %.1f$' % self.nu
        if self.app == 'FP_pp': info += r', $r_s = %.1f$' % self.rs
        return info

    def load_file(self, a_file, run_date):
        with open(a_file) as data_file:
            data = json.loads(data_file.read())

        self.num_g = data["Ng"]
        self.num_p = data["par_num"]**(1./3)
        self.num_m = data["mesh_num"]
        self.box = data["box_size"]
        self.nu = data["viscosity"]
        self.rs = data["cut_radius"]
        self.app = data["app"]
        self.results = data["results"]
        if self.results is None:
            self.results = {}

        self.file = a_file
        self.dir = a_file.replace(a_file.split("/")[-1], '')
        self.res_dir = self.dir + 'results/'
        create_dir(self.res_dir)

    def done(self, key):
        with open(self.file) as data_file:
            data = json.loads(data_file.read())
        if data["results"] is None:
            data["results"] = {}

        data["results"][key] = True

        with open(self.file, 'w') as outfile:
            json.dump(data, outfile, indent=2)

def load_k_supp(files, k_l_lim=[0,10], k_m_lim=[20,27], k_s_lim=[35,40]):
    """
    k_min = 2*PI / L
    k_max = 2*PI / L * N_m
    log_bin : 100 values between k_min and k_max (logspace)

    large scale :: k = k_min * (log_bin^[0:10])
    medium scale :: k = k_min * (log_bin^[15:25])
    small scale :: k = k_min * (log_bin^[30:40]) """

    supp_large = []
    supp_medium = []
    supp_small = []

    for a_file in files:
        data = np.loadtxt(a_file)
        P_diff = data[:, 1]
        supp_large.append(np.mean(P_diff[k_l_lim]))
        supp_medium.append(np.mean(P_diff[k_m_lim]))
        supp_small.append(np.mean(P_diff[k_s_lim]))
        del data, P_diff

    data = np.loadtxt(files[-1])
    k = data[:, 0]
    k_l = [k[k_l_lim[0]], k[k_l_lim[1]]]
    k_m = [k[k_m_lim[0]], k[k_m_lim[1]]]
    k_s = [k[k_s_lim[0]], k[k_s_lim[1]]]
    return (supp_large, supp_medium, supp_small), (k_l, k_m, k_s)

def analyze_run(a_sim_info, rerun=False, skip_ani=False):
    if rerun or not a_sim_info.results:
        plt.rcParams['legend.numpoints'] = 1

        # Power spectrum
        zs, files = try_get_zs_files(a_sim_info, 'pwr_spec/')
        if zs is not None:
            print 'Plotting power spectrum...'
            plot.plot_pwr_spec(files, zs, a_sim_info)
            a_sim_info.done("pwr_spec")
        del zs, files

        # Power spectrum difference
        zs, files = try_get_zs_files(a_sim_info, 'pwr_diff/')
        if zs is not None:
            print 'Plotting power spectrum difference...'
            plot.plot_pwr_spec_diff(files, zs, a_sim_info)
            a_sim_info.done("pwr_spec_diff")
        # Power spectrum suppression
            print 'Plotting power spectrum suppression...'
            a = [1./(z+1) for z in zs]
            supp_lms, k_lms = load_k_supp(files)
            plot.plot_supp_lms(supp_lms, a, a_sim_info, k_lms=k_lms)
            a_sim_info.done("pwr_spec_supp")
        del zs, files

        # Density distribution
        zs, files = try_get_zs_files(a_sim_info, 'rho_bin/')
        if zs is not None:
            print 'Plotting density distribution...'
            zs, files = slice_zs_files(zs, files)
            plot.plot_dens_histo(files, zs, a_sim_info)
            a_sim_info.done("dens_hist")
        del zs, files

        # Particles evolution
        zs, files = try_get_zs_files(a_sim_info, 'par_cut/', a_file='par*.dat')
        zs_t, files_t = try_get_zs_files(a_sim_info, 'par_cut/', a_file='track*.dat')
        if zs is not None:
            if zs != zs_t: print "ERROR! 'par_cut' files differ from 'track_par_pos' files. Skipping step."
            else:
                # last slice
                print 'Plotting slice through simulation box (particles)...'
                plot.plot_par_last_slice(files, files_t, zs, a_sim_info)
                a_sim_info.done("par_slice")
                # animation
                if not skip_ani:
                    print 'Plotting particles evolution...'
                    plot.plot_par_evol(files, files_t, zs, a_sim_info)
                    a_sim_info.done("par_ani")
        del zs, files, zs_t, files_t

        # Density evolution
        zs, files = try_get_zs_files(a_sim_info, 'rho_map/', a_file='*.dat')
        if zs is not None:
            # two slices
            print 'Plotting two slices through simulation box (overdensity)...'
            plot.plot_dens_two_slices(files, zs, a_sim_info)
            a_sim_info.done("dens_slice")
            # animation
            if not skip_ani:
                print 'Plotting density evolution...'
                plot.plot_dens_evol(files, zs, a_sim_info)
                a_sim_info.done("dens_ani")
        del zs, files

        # Clean
        plt.close("all")
        gc.collect()
    else:
        print 'Run already analyzed!'

def analyze_all(out_dir='/home/vrastil/Documents/GIT/Adhesion-Approximation/output/',
                rerun=False, skip_ani=False, only=None):
    files = get_files_in_traverse_dir(out_dir, 'sim_param.json')
    sim_infos = []
    for args in files:
        sim_infos.append(SimInfo(*args))

    if only is not None:
        sim_infos = sim_infos[only]

    for a_sim_info in sim_infos:
        print 'Analyzing run %s' % a_sim_info.info_tr()
        analyze_run(a_sim_info, rerun=rerun, skip_ani=skip_ani)
    print 'All runs analyzed!'

if __name__ == 'main':
    out_dir = sys.argv[1]
    analyze_all(out_dir)
