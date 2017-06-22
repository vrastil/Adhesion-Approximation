from datetime import datetime

from . import *
from . import plot

class SimInfo(object):
    def __init__(self, *args):
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
            self.results = False

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

def analyze_run(a_sim_info, rerun=False):
    if rerun or not a_sim_info.results:
        out_dir = a_sim_info.dir + 'results/'
        create_dir(out_dir)

        # Power spectrum
        print 'Plotting power spectrum...'
        zs, files = sort_get_fl_get_z(a_sim_info, 'pwr_spec/')
        plot.plot_pwr_spec(files, zs, a_sim_info, out_dir)

        # Power spectrum difference
        print 'Plotting power spectrum difference...'
        zs, files = sort_get_fl_get_z(a_sim_info, 'pwr_diff/')
        plot.plot_pwr_spec_diff(files, zs, a_sim_info, out_dir)

        # Power spectrum suppresion
        print 'Plotting power spectrum suppresion...'
        plot.plot_supp([a_sim_info], out_dir)

        # Density distribution
        print 'Plotting density distribution...'
        zs, files = slice_zs_files(*sort_get_fl_get_z(a_sim_info, 'rho_bin/'))
        plot.plot_dens_histo(files, zs, a_sim_info, out_dir)

        # Particles evolution
        print 'Plotting particles evolution...'
        zs, files = sort_get_fl_get_z(a_sim_info, 'par_cut/', a_file='par*.dat')
        zs_t, files_t = sort_get_fl_get_z(a_sim_info, 'par_cut/', a_file='track*.dat')
        if zs != zs_t: print "ERROR! 'par_cut' files differ from 'track_par_pos' files"
        ani_par = plot.plot_par_evol(files, files_t, zs, a_sim_info, out_dir)

        # Density evolution
        print 'Plotting density evolution...'
        zs, files = sort_get_fl_get_z(a_sim_info, 'rho_map/', a_file='*.dat')
        ani_evol = plot.plot_dens_evol(files, zs, a_sim_info, out_dir)

        # Update sim_param.log
        print "Updating 'sim_param.log'..."
        a_sim_info.done()
    else:
        print 'Run already analyzed!'

    #return ani_par, ani_evol


def analyze_all(out_dir='/home/vrastil/Documents/Adhesion-Approximation/output/'):
    files = get_files_in_traverse_dir(out_dir, 'sim_param.log')
    sim_infos = []
    for args in files:
        sim_infos.append(SimInfo(*args))

    for a_sim_info in sim_infos:
        print 'Analyzing run %s' % a_sim_info.info_tr()
        analyze_run(a_sim_info)

