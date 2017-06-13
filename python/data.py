import os, fnmatch
from datetime import datetime

from . import plot

def get_files_in_traverse_dir(a_dir, a_file):
    """ return list of all files in directory and its subdirectories \
    which matches 'a_file' and its subdirectory path, support Unix \
    filename pattern matching ('*', '?', [seq], [!seq]) """

    ls_file = []
    for root, dirs, files in os.walk(a_dir):
        for name in files:
            if fnmatch.fnmatch(name, a_file):
                subdir = root.replace(a_dir, '')
                ls_file.append((os.path.join(root, name), subdir))
    return ls_file

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

    def info(self):
        info = ''
        info += '%s:\n' % self.app
        info += '$N_m = %i^3$\n' % self.num_m
        info += '$N_p = %i^3$\n' % self.num_p
        info += '$L = %i$ Mpc/h\n' % self.box
        if self.app == 'AA': info += '$\nu = %.1f$\n' % self.nu
        if self.app == 'FP_pp': info += '$r_s = %.1f$\n' % self.rs
        return info

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

        self.app = run_date.split('/')[0][:-4]
        self.date = datetime.strptime(run_date.split('/')[1], '%Y_%m_%d.%H:%M:%S')
        self.dir = a_file.replace('sim_param.log', '')

def sort_get_z(files, a_sim_info):
    zs = []
    for a_file in files:
        if a_sim_info.app + '_z' in a_file:
            z = float(a_file[a_file.index(a_sim_info.app + '_z') + len(a_sim_info.app+'_z'):-4])
        elif a_sim_info.app + '_init' in a_file:
            z = 'init'
        else:
            print "WARNING! Skipping file '%s', unknown format." % a_file
        zs.append(z)
    return zip(*sorted(zip(zs, files), reverse=True))

def sort_get_fl_get_z(a_sim_info, subdir):
    files = [x[0] for x in get_files_in_traverse_dir(a_sim_info.dir + subdir, '*.dat')]
    return sort_get_z(files, a_sim_info)

def analyze_run(a_sim_info):
    # Power spectrum
    zs, files = sort_get_fl_get_z(a_sim_info, 'pwr_spec/')

    plot.plot_pwr_spec(files, zs, a_sim_info)
    # Power spectrum difference
    zs, files = sort_get_fl_get_z(a_sim_info, 'pwr_diff/')
    plot.plot_pwr_spec_diff(files, zs, a_sim_info)

    # Power spectrum suppresion



def analyze_all(out_dir='/home/vrastil/Documents/Adhesion-Approximation/output/'):
    files = get_files_in_traverse_dir(out_dir, 'sim_param.log')
    sim_infos = []
    for args in files:
        sim_infos.append(SimInfo(*args))

    for a_sim_info in sim_infos:
        analyze_run(a_sim_info)

