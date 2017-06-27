from datetime import datetime
import sys
import gc
import matplotlib.pyplot as plt

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
            self.date = ''
            self.dir = ''

    def info(self):
        info = ''
        info += '%s:\n' % self.app
        info += '$N_m = %i^3$\n' % self.num_m
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


        self.app = run_date.split('/')[0][:-4]
        self.date = datetime.strptime(run_date.split('/')[1], '%Y_%m_%d.%H:%M:%S')
        self.dir = a_file.replace('sim_param.log', '')

    def done(self):
        with open(self.dir + 'sim_param.log', 'r+') as f:
            content = f.read()
            if content.startswith('Results: Done'):
                return None
            f.seek(0, 0)
            f.write('Results: Done\n' + content)

def analyze_run(a_sim_info, rerun=False, skip_ani=False):
    if rerun or not a_sim_info.results:
        out_dir = a_sim_info.dir + 'results/'
        create_dir(out_dir)
        all_good = True
        # Power spectrum
        try:
            zs, files = sort_get_fl_get_z(a_sim_info, 'pwr_spec/')
        except ValueError:
            print "WARNING! Missing data in '%s'. Skipping step." % (a_sim_info.dir + 'pwr_spec/')
            all_good = False
        else:
            print 'Plotting power spectrum...'
            plot.plot_pwr_spec(files, zs, a_sim_info, out_dir)

        # Power spectrum difference
        try:
            zs, files = sort_get_fl_get_z(a_sim_info, 'pwr_diff/')
        except ValueError:
            print "WARNING! Missing data in '%s'. Skipping step." % (a_sim_info.dir + 'pwr_diff/')
            all_good = False
        else:
            print 'Plotting power spectrum difference...'
            plot.plot_pwr_spec_diff(files, zs, a_sim_info, out_dir)

        # Power spectrum suppresion
            print 'Plotting power spectrum suppresion...'
            plot.plot_supp([a_sim_info], out_dir)

        # Density distribution
        try:
            zs, files = slice_zs_files(*sort_get_fl_get_z(a_sim_info, 'rho_bin/'))
        except TypeError:
            print "WARNING! Missing data in '%s'. Skipping step." % (a_sim_info.dir + 'rho_bin/')
            all_good = False
        else:
            print 'Plotting density distribution...'
            plot.plot_dens_histo(files, zs, a_sim_info, out_dir)

        if not skip_ani:
            # Particles evolution
            try:
                zs, files = sort_get_fl_get_z(a_sim_info, 'par_cut/', a_file='par*.dat')
                zs_t, files_t = sort_get_fl_get_z(a_sim_info, 'par_cut/', a_file='track*.dat')
            except ValueError:
                print "WARNING! Missing data in '%s'. Skipping step." % (a_sim_info.dir + 'par_cut/')
                all_good = False
            else:
                if zs != zs_t: print "ERROR! 'par_cut' files differ from 'track_par_pos' files. Skipping step."
                else:
                    print 'Plotting particles evolution...'
                    plot.plot_par_evol(files, files_t, zs, a_sim_info, out_dir)

            # Density evolution
            try:
                zs, files = sort_get_fl_get_z(a_sim_info, 'rho_map/', a_file='*.dat')
            except ValueError:
                print "WARNING! Missing data in '%s'. Skipping step." % (a_sim_info.dir + 'rho_map/')
                all_good = False
            else:
                print 'Plotting density evolution...'
                plot.plot_dens_evol(files, zs, a_sim_info, out_dir)

        # Update sim_param.log
        if all_good:
            print "Updating 'sim_param.log'..."
            a_sim_info.done()

        # Clean
        plt.close("all")
        gc.collect()
    else:
        print 'Run already analyzed!'

def analyze_all(out_dir='/home/vrastil/Documents/Adhesion-Approximation/output/'):
    files = get_files_in_traverse_dir(out_dir, 'sim_param.log')
    sim_infos = []
    for args in files:
        sim_infos.append(SimInfo(*args))

    for a_sim_info in sim_infos:
        print 'Analyzing run %s' % a_sim_info.info_tr()
        analyze_run(a_sim_info)
    print 'All runs analyzed!'

if __name__ == 'main':
    out_dir = sys.argv[1]
    analyze_all(out_dir)