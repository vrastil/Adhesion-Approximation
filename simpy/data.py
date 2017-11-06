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
    def __init__(self, *args):
        self.results = False
        if len(args) == 2 and args[0].endswith('.json'):
            self.load_file(*args)
        elif len(args) == 2 and  args[0].endswith('.log'):
            self.load_file_log(*args)
        else:
            print "WARNING! Invalid simulation parameters file '%s'." % args[0]

        self.cosmo = None
        if hasattr(self, 'pwr'):
            if self.pwr["pwr_type"] >= 4:
                import pyccl as ccl
                self.cosmo = ccl.Cosmology(
                    Omega_c=self.pwr["Omega_c"], Omega_b=self.pwr["Omega_b"], h=self.pwr["h"],
                    sigma8=self.pwr["sigma8"], n_s=self.pwr["ns"],
                    transfer_function=self.ccl["transfer_function_method"],
                    matter_power_spectrum=self.ccl["matter_power_spectrum_method"],
                    mass_function=self.ccl["mass_function_method"])

    def info(self):
        info = ''
        info += '%s:\n' % self.app
        info += '$N_m = %i$\n' % self.num_m
        info += '$N_M = %i$\n' % self.num_M
        info += '$N_p = %i^3$\n' % self.num_p
        info += '$L = %i$ Mpc/h\n' % self.box
        if self.app == 'AA': info += r'$\nu = %.1f$ (Mpc/h)$^2$' % self.nu
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
        self.num_p = data["par_num"]
        self.num_m = data["mesh_num"]
        self.num_M = data["mesh_num_pwr"]
        self.box = data["box_size"]
        self.nu = data["viscosity"]
        self.rs = data["cut_radius"]
        self.app = data["app"]
        self.pwr = {}
        self.pwr["A"] = data["A"]
        self.pwr["ns"] = data["index"]
        self.pwr["k2_G"] = data["smoothing_k"]
        self.pwr["sigma8"] = data["sigma8"]
        self.pwr["Omega_c"] = data["Omega_c"]
        self.pwr["Omega_b"] = data["Omega_b"]
        self.pwr["h"] = data["h"]
        self.pwr["pwr_type"] = data["pwr_type"]
        if "k_pade" in data: self.k_pade = data["k_pade"]
        else: self.k_pade = None
        if "k_nyquist" in data: self.k_nyquist = data["k_nyquist"]
        else: self.k_nyquist = None
        if "seed" in data:
            self.seed = data["seed"]
        self.ccl = {}
        self.ccl["transfer_function_method"] = data["transfer_function_method"]
        self.ccl["matter_power_spectrum_method"] = data["matter_power_spectrum_method"]
        self.ccl["mass_function_method"] = data["mass_function_method"]

        self.results = data["results"]
        if self.results is None:
            self.results = {}
        for key in RESULTS_KEYS:
            if key not in self.results:
                self.results[key] = False

        self.file = a_file
        self.dir = a_file.replace(a_file.split("/")[-1], '')
        self.res_dir = self.dir + 'results/'
        create_dir(self.res_dir)

    def load_file_log(self, a_file, run_date):
        with open(a_file, 'r') as f:
            content = f.readlines()

        for line in content:
            if line.startswith('Ng'):
                self.num_g = int(line.replace('Ng:\t\t', '').rstrip())
            elif line.startswith('Num_par'):
                self.num_p = int(line.replace('Num_par:\t', '').rstrip()[:-2])
            elif line.startswith('Num_mesh'):
                self.num_m = int(line.replace('Num_mesh:\t', '').rstrip()[:-2])
            elif line.startswith('Box size'):
                self.box = int(line.replace('Box size:\t', '').rstrip()[:-6])
            elif 'nu = ' in line:
                self.nu = float(line[line.index('nu =') + 5:-6])
            elif 'rs = ' in line:
                self.rs = float(line[line.index('rs =') + 5:line.index(',')])
            elif line.startswith('Results: Done'):
                self.results = True

        self.num_M = self.num_m
        self.file = a_file
        self.app = run_date.split('/')[0][:-4]
        self.date = datetime.strptime(run_date.split('/')[1], '%Y_%m_%d.%H:%M:%S')
        self.dir = a_file.replace('sim_param.log', '')
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

        idx = (np.abs(k-0.67*k_nyquist_par)).argmin() / 7

        supp_large.append(np.mean(P_diff[0:idx]))
        supp_large_std.append(np.std(P_diff[0:idx]))

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
    if a_sim_info.rerun(rerun, key, skip, zs):
        print 'Plotting power spectrum...'
        plot.plot_pwr_spec(files, zs, a_sim_info, pwr_spec_files_extrap=files_extrap, pwr_spec_files_emu=files_emu)
        a_sim_info.done(key)
    del zs, files, files_extrap, files_emu
    
    # Power spectrum difference
    key = "pwr_spec_diff"
    zs, files = try_get_zs_files(a_sim_info, 'pwr_diff/')
    if a_sim_info.rerun(rerun, key, skip, zs):
        print 'Plotting power spectrum difference...'
        plot.plot_pwr_spec_diff(files, zs, a_sim_info)
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
    zs, files_emu = try_get_zs_files(a_sim_info, 'corr_func/', a_file='*gsl*emu*.dat')
    zs, files_lin = try_get_zs_files(a_sim_info, 'corr_func/', a_file='*gsl*lin*.dat')
    zs, files = try_get_zs_files(a_sim_info, 'corr_func/', a_file='*gsl*par*.dat')
    if a_sim_info.rerun(rerun, key, skip, zs):
        print 'Plotting correlation function...'
        plot.plot_corr_func(files, zs, a_sim_info, corr_func_files_lin=files_lin, corr_func_files_emu=files_emu)
        a_sim_info.done(key)
    del zs, files, files_lin

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
    for args in files:
        sim_infos.append(SimInfo(*args))

    if only is not None:
        sim_infos = sim_infos[only]

    for a_sim_info in sim_infos:
        print 'Analyzing run %s' % a_sim_info.info_tr()
        analyze_run(a_sim_info, rerun=rerun, skip=skip)
    print 'All runs analyzed!'

if __name__ == 'main':
    out_dir = sys.argv[1]
    analyze_all(out_dir)
