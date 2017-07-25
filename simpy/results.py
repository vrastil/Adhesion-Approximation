import subprocess
from IPython.display import Image, display

from . import *
from .data import SimInfo, load_k_supp
from .plot import plot_supp

class Results(object):
    """
    Class to handle all the results. Find, sort, show.
    """
    def __init__(self, out_dir=''):
        if out_dir != '':
            self.load(out_dir)
        else: sim_infos = []

    def sort(self):
        self.sim_infos = sorted(self.sim_infos, key=lambda si: si.nu)
        self.sim_infos = sorted(self.sim_infos, key=lambda si: si.num_m)
        self.sim_infos = sorted(self.sim_infos, key=lambda si: si.app)

    def load(self, out_dir):
        files = get_files_in_traverse_dir(out_dir, 'sim_param.log')
        self.sim_infos = []
        for args in files:
            self.sim_infos.append(SimInfo(*args))
        self.sort()

    def get_subfiles(self, Nm=0, Np=0, L=0, nu=0, rs=0, app=''):
        subfiles = []
        for a_sim_info in self.sim_infos:
            if (
                    (Nm == 0 or a_sim_info.num_m == Nm) and
                    (Np == 0 or a_sim_info.num_p == Np) and
                    (L == 0 or a_sim_info.box == L) and
                    (nu == 0 or a_sim_info.nu == nu) and
                    (rs == 0 or a_sim_info.rs == rs) and
                    (app == '' or a_sim_info.app == app)
                ):
                subfiles.append(a_sim_info)

        return subfiles

    def info(self, Nm=0, Np=0, L=0, nu=0, rs=0, app=''):
        for a_sim_info in self.get_subfiles(Nm=Nm, Np=Np, L=L, nu=nu, rs=rs, app=app):
            print a_sim_info.info_tr()

    def show_folder(self, a_sim_info):
        subprocess.Popen(["xdg-open", a_sim_info.dir + 'results/'])

    def show_results(self, Nm=0, Np=0, L=0, nu=0, rs=0, app='',
                     pwr_spec=True, pwr_spec_diff=True, supp=True, dens_histo=True,
                     dens_last=True, par_evol_last=True):
        for a_sim_info in self.get_subfiles(Nm=Nm, Np=Np, L=L, nu=nu, rs=rs, app=app):
            results_dir = a_sim_info.dir + 'results/'
            if pwr_spec: display(Image(filename=results_dir + "pwr_spec.png"))
            if pwr_spec_diff: display(Image(filename=results_dir + "pwr_spec_diff.png"))
            if supp: display(Image(filename=results_dir + "supp.png"))
            if dens_histo: display(Image(filename=results_dir + "dens_histo.png"))
            if dens_last: display(Image(filename=results_dir + "dens_z0.00.png"))
            if par_evol_last: display(Image(filename=results_dir + "par_evol_last.png"))

    def load_k_supp(self, a_sim_info):
        if not hasattr(a_sim_info, "supp"):
            zs, files = try_get_zs_files(a_sim_info, 'pwr_diff/')
            if zs is not None:
                a_sim_info.supp = load_k_supp(files)
                a_sim_info.a = [1./(z+1) for z in zs]

    def plot_supp_compare(self, out_dir='/home/vrastil/Documents/GIT/Adhesion-Approximation/output/supp_comparison/',
                          Nm=0, Np=0, L=0, nu=0, rs=0, app='', scale=['small', 'medium', 'large'], show_k_lms=False):
        subfiles = self.get_subfiles(Nm=Nm, Np=Np, L=L, nu=nu, rs=rs, app=app)
        for a_sim_info in subfiles:
            self.load_k_supp(a_sim_info)
        for sc in scale:
            plot_supp(subfiles, out_dir+sc, suptitle=' on %s scales' % sc, save=False, show=True, scale=sc, show_k_lms=show_k_lms)
                    