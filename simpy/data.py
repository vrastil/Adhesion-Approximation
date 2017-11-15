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

RESULTS_KEYS = ["pwr_spec", "pwr_spec_diff", "pwr_spec_supp", "dens_hist",
                "vel_pwr_spec", "vel_pwr_spec_diff", "vel_pwr_spec_supp",
                "par_slice", "par_ani", "dens_slice", "dens_ani", "corr_func"]

class SimInfo(object):
    def __init__(self, a_file):
        if a_file.endswith('.json'):
            self.load_file(a_file)
        else:
            print "WARNING! Invalid simulation parameters file '%s'." % a_file

        if self.cosmo["pwr_type"] >= 4:
            import pyccl as ccl
            self.ccl_cosmo = ccl.Cosmology(
                Omega_c=self.cosmo["Omega_c"], Omega_b=self.cosmo["Omega_b"], h=self.cosmo["h"],
                sigma8=self.cosmo["sigma8"], n_s=self.cosmo["index"],
                transfer_function=self.cosmo["transfer_function_method"],
                matter_power_spectrum=self.cosmo["matter_power_spectrum_method"],
                mass_function=self.cosmo["mass_function_method"])

    def info(self):
        info = ''
        info += '%s:\n' % self.app
        info += '$N_m = %i$\n' % self.box_opt["mesh_num"]
        info += '$N_M = %i$\n' % self.box_opt["mesh_num_pwr"]
        info += '$N_p = %i^3$\n' % self.box_opt["par_num"]
        info += '$L = %i$ Mpc/h\n' % self.box_opt["box_size"]
        if self.app == 'AA': info += r'$\nu = %.1f$ (Mpc/h)$^2$' % self.app_opt["viscosity"]
        if self.app == 'FP_pp': info += r'$r_s = %.1f$' % self.app_opt["cut_radius"]
        return info

    def info_tr(self):
        info = self.info()
        return info.replace('\n', '  ')

    def info_supp(self):
        info = '%s: $res = %.2f$' % (self.app, float(self.box_opt["mesh_num"]) / self.box_opt["box_size"])
        if self.app == 'AA': info += r', $\nu = %.1f$' % self.app_opt["viscosity"]
        if self.app == 'FP_pp': info += r', $r_s = %.1f$' % self.app_opt["cut_radius"]
        return info

    def load_file(self, a_file):
        with open(a_file) as data_file:
            data = json.loads(data_file.read())

        for key in data:
            setattr(self, key, data[key])

        if self.results is None:
            self.results = {}
        for key in RESULTS_KEYS:
            if key not in self.results:
                self.results[key] = False

        self.file = a_file
        self.dir = a_file.replace(a_file.split("/")[-1], '')
        self.res_dir = self.dir + 'results/'
        create_dir(self.res_dir)

    def rerun(self, rerun, key, skip, zs):
        if zs is None or key in skip:
            return False
        elif rerun == "all":
            return True
        else:
            return not self.results[key] or key in rerun

    def done(self, key):
        with open(self.file) as data_file:
            data = json.loads(data_file.read())
        if data["results"] is None:
            data["results"] = {}

        data["results"][key] = True

        with open(self.file, 'w') as outfile:
            json.dump(data, outfile, indent=2)

def load_k_supp(files, k_nyquist_par):
    """
    Divide available k-values into 7 subinterval from
    k_min = 2*PI / L to k_max = 67% k_nyquist_par
    large scale :: k = 1st subinterval
    medium scale :: k = 4rd subinterval
    small scale :: k = 7th subinterval """

    supp_large = []
    supp_medium = []
    supp_small = []

    supp_large_std = []
    supp_medium_std = []
    supp_small_std = []

    for a_file in files:
        data = np.loadtxt(a_file)
        k = data[:, 0]
        P_diff = data[:, 1]

        idx = (np.abs(k-0.5*k_nyquist_par)).argmin() / 7

        supp_large.append(np.mean(P_diff[0*idx:1*idx]))
        supp_large_std.append(np.std(P_diff[0*idx:1*idx]))

        supp_medium.append(np.mean(P_diff[3*idx:4*idx]))
        supp_medium_std.append(np.std(P_diff[3*idx:4*idx]))

        supp_small.append(np.mean(P_diff[6*idx:7*idx]))        
        supp_small_std.append(np.std(P_diff[6*idx:7*idx]))

        del data, k, P_diff

    data = np.loadtxt(files[-1])
    k = data[:, 0]
    k_l = [k[0*idx], k[1*idx]]
    k_m = [k[3*idx], k[4*idx]]
    k_s = [k[6*idx], k[7*idx]]
    return (supp_large, supp_medium, supp_small), (supp_large_std, supp_medium_std, supp_small_std), (k_l, k_m, k_s)

def analyze_run(a_sim_info, rerun=None, skip=None):
    plt.rcParams['legend.numpoints'] = 1

    if skip is None:
        skip = []
    elif skip == "ani":
        skip = ["par_ani", "dens_ani"]

    if rerun is None:
        rerun = []

    # Power spectrum
    key = "pwr_spec"
    zs, files_emu = try_get_zs_files(a_sim_info, 'pwr_spec/', a_file='*emu*.dat')
    zs, files_extrap = try_get_zs_files(a_sim_info, 'pwr_spec/', a_file='*extrap*.dat')
    zs, files = try_get_zs_files(a_sim_info, 'pwr_spec/', a_file='*par*.dat')
    zs_init, file_init = try_get_zs_files(a_sim_info, 'pwr_spec/', a_file='*init*.dat')
    if a_sim_info.rerun(rerun, key, skip, zs):
        if zs_init is not None:
            # insert at the begginnig -- important for iteration
            zs.insert(0, zs_init[0])
            files.insert(0, file_init[0])
        print 'Plotting power spectrum...'
        plot.plot_pwr_spec(files, zs, a_sim_info, pwr_spec_files_extrap=files_extrap, pwr_spec_files_emu=files_emu)
        a_sim_info.done(key)
    del zs, zs_init, files, files_extrap, files_emu, file_init
    
    # Power spectrum difference -- input
    key = "pwr_spec_diff"
    print 'Plotting power spectrum difference...'
    zs, files = try_get_zs_files(a_sim_info, 'pwr_diff/', a_file='*input*')
    if a_sim_info.rerun(rerun, key, skip, zs):
        plot.plot_pwr_spec_diff(files, zs, a_sim_info, ext_title='input')
        a_sim_info.done(key)
    # Power spectrum difference -- hybrid
    zs, files = try_get_zs_files(a_sim_info, 'pwr_diff/', a_file='*hybrid*')
    if a_sim_info.rerun(rerun, key, skip, zs):
        plot.plot_pwr_spec_diff(files, zs, a_sim_info, ext_title='hybrid')
        a_sim_info.done(key)
    # Power spectrum difference -- input
    zs, files = try_get_zs_files(a_sim_info, 'pwr_diff/', a_file='*par*')
    if a_sim_info.rerun(rerun, key, skip, zs):
        plot.plot_pwr_spec_diff(files, zs, a_sim_info, ext_title='par')
        a_sim_info.done(key)
    # Power spectrum suppression
    key = "pwr_spec_supp"
    if a_sim_info.rerun(rerun, key, skip, zs):
        print 'Plotting power spectrum suppression...'
        a = [1./(z+1) for z in zs]
        supp_lms, supp_std_lms, k_lms = load_k_supp(files, a_sim_info.k_nyquist["particle"])
        plot.plot_supp_lms(supp_lms, a, a_sim_info, k_lms, supp_std_lms)
        a_sim_info.done(key)
    del zs, files

    # Velocity power spectrum
    key = "vel_pwr_spec"
    zs, files = try_get_zs_files(a_sim_info, 'vel_pwr_spec/')
    if a_sim_info.rerun(rerun, key, skip, zs):
        print 'Plotting velocity power spectrum...'
        plot.plot_pwr_spec(files, zs, a_sim_info, pk_type='vel')
        a_sim_info.done(key)
    del zs, files
    # Velocity power spectrum difference
    key = "vel_pwr_spec_diff"
    zs, files = try_get_zs_files(a_sim_info, 'vel_pwr_diff/')
    if a_sim_info.rerun(rerun, key, skip, zs):
        print 'Plotting velocity power spectrum difference...'
        plot.plot_pwr_spec_diff(files, zs, a_sim_info, pk_type='vel')
        a_sim_info.done(key)
    # Velocity power spectrum suppression
    key = "vel_pwr_spec_supp"
    if a_sim_info.rerun(rerun, key, skip, zs):
        print 'Plotting power spectrum suppression...'
        a = [1./(z+1) for z in zs]
        supp_lms, supp_std_lms, k_lms = load_k_supp(files, a_sim_info.k_nyquist["particle"])
        plot.plot_supp_lms(supp_lms, a, a_sim_info, k_lms, supp_std_lms, pk_type='vel')
        a_sim_info.done(key)
    del zs, files

    # Correlation function
    key = "corr_func"
    zs_emu, files_emu = try_get_zs_files(a_sim_info, 'corr_func/', a_file='*gsl*emu*.dat')
    zs, files_lin = try_get_zs_files(a_sim_info, 'corr_func/', a_file='*gsl*lin*.dat')
    zs, files = try_get_zs_files(a_sim_info, 'corr_func/', a_file='*gsl*par*.dat')
    if a_sim_info.rerun(rerun, key, skip, zs):
        print 'Plotting correlation function...'
        plot.plot_corr_func(files, zs, a_sim_info, corr_func_files_lin=files_lin, corr_func_files_emu=files_emu, zs_emu=zs_emu)
        a_sim_info.done(key)
    del zs, zs_emu, files, files_lin

    # Density distribution
    key = "dens_hist"
    zs, files = try_get_zs_files(a_sim_info, 'rho_bin/')
    if a_sim_info.rerun(rerun, key, skip, zs):
        print 'Plotting density distribution...'
        zs, files = slice_zs_files(zs, files)
        plot.plot_dens_histo(files, zs, a_sim_info)
        a_sim_info.done(key)
    del zs, files

    # Particles evolution
    zs, files = try_get_zs_files(a_sim_info, 'par_cut/', a_file='par*.dat')
    zs_t, files_t = try_get_zs_files(a_sim_info, 'par_cut/', a_file='track*.dat')
    if zs != zs_t: print "ERROR! 'par_cut' files differ from 'track_par_pos' files. Skipping step."
    else:
        # last slice
        key = "par_slice"
        if a_sim_info.rerun(rerun, key, skip, zs):
            print 'Plotting slice through simulation box (particles)...'
            plot.plot_par_last_slice(files, files_t, zs, a_sim_info)
            a_sim_info.done(key)
        # animation
        key = "par_ani"
        if a_sim_info.rerun(rerun, key, skip, zs):
            print 'Plotting particles evolution...'
            plot.plot_par_evol(files, files_t, zs, a_sim_info)
            a_sim_info.done(key)
    del zs, files, zs_t, files_t

    # Density evolution
    zs, files = try_get_zs_files(a_sim_info, 'rho_map/')
    # two slices
    key = "dens_slice"
    if a_sim_info.rerun(rerun, key, skip, zs):
        print 'Plotting two slices through simulation box (overdensity)...'
        plot.plot_dens_two_slices(files, zs, a_sim_info)
        a_sim_info.done(key)
    # animation
    key = "dens_ani"
    if a_sim_info.rerun(rerun, key, skip, zs):
        print 'Plotting density evolution...'
        plot.plot_dens_evol(files, zs, a_sim_info)
        a_sim_info.done(key)
    del zs, files

    # Clean
    plt.close("all")
    gc.collect()

def analyze_all(out_dir='/home/vrastil/Documents/GIT/Adhesion-Approximation/output/',
                rerun=None, skip=None, only=None):
    files = get_files_in_traverse_dir(out_dir, 'sim_param.json')
    sim_infos = []
    for a_file, subdir in files:
        sim_infos.append(SimInfo(a_file))

    if only is not None:
        sim_infos = sim_infos[only]

    for a_sim_info in sim_infos:
        print 'Analyzing run %s' % a_sim_info.info_tr()
        analyze_run(a_sim_info, rerun=rerun, skip=skip)
    print 'All runs analyzed!'

if __name__ == 'main':
    out_dir = sys.argv[1]
    analyze_all(out_dir)
