"""
'struct.py' module defines all necessary data containers for working with results

 - Sim_Info : basic class storing all information about one particular simulation run
 - Stack_Info : subclass of Sim_Info, collects and stores information about all runs
   with the same parameters
 - Results : collects results, compare approximation, show plots
"""

from __future__ import print_function

import json
import os
import fnmatch
import subprocess
from IPython.display import Image, display
from .fastsim import Sim_Param

RESULTS_ALL = {
    "ani" : ["par_ani", "dens_ani"],
    "corr" : ["corr_func", "bao", "sigma_R"],
    "dens" : ["dens_hist", "dens_slice"],
    "chi" : ["pwr_spec_chi", "chi_pwr_diff", "chi_pwr_spec_supp", "chi_pwr_spec_supp_map"],
    "par" : ["par_slice"],
    "pwr" : ["pwr_spec", "pwr_diff", "pwr_diff_i", "pwr_diff_h", "pwr_spec_supp", "pwr_spec_supp_map", "pwr_slope"],
    "vel" : ["vel_pwr_spec", "vel_pwr_diff", "vel_pwr_spec_supp"],
    "eff_time" : ["Pk", "sigma_R"],
    "files" : ["corr_files", "sigma_files", "pwr_spec_files", "pwr_diff_files",
               "pwr_diff_files_i", "pwr_diff_files_h", "pwr_spec_chi_files"]
    }

def _is_key_strval(key, strval):
    """ check <strval> (rerun or skip) against <key>, true either for
        1) exact match
        2) <strval> stands for group of keys one of which is 'key'
        3) <strval> is 'all'
        """
    if (strval == key) or (
        strval in RESULTS_ALL and key in RESULTS_ALL[strval]) or (
        strval == 'all'):
        return True
    else:
        return False

def _is_key_val(key, val):
    """ check val (rerun or skip) against key and return bool """
    # single string
    if isinstance(val, str):
        return _is_key_strval(key, val)
    # list of strings
    elif isinstance(val, list):
        for strval in val:
            if _is_key_strval(key, strval):
                return True
    # None or unknown type
    else:
        return False

def get_files_in_traverse_dir(a_dir, patterns):
    # type: (str, str) -> List[str]
    """ return list of all files in directory which matches 'patterns'
    support Unix filename pattern matching ('*', '?', [seq], [!seq])
    and multiple option in 'patterns' (space delimetered) """

    return list(set([ # throw away duplicate files
        os.path.join(root, name) # full file name
        for root, _, files in os.walk(a_dir) # go through all subdirectores
        for pattern in patterns.split() # if multiple patterns given
        for name in fnmatch.filter(files, pattern) # pattern matching
        ]))

def create_dir(out_dir):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

class SimInfo(object):
    """ basic class storing all information about one particular simulation run """
    def __init__(self, a_file, verbose=True, **kwargs):
        # type: (str, dict) -> None
        """ load information in 'a_file', replace parameters in 'self.cosmo'
        by any additional parameters in 'kwargs' (optional) """
        # attributes to load from json file: string & dictionaries
        self.app = ""
        self.app_opt, self.box_opt, self.cosmo, self.integ_opt, self.k_nyquist, self.out_opt, self.results, self.run_opt, self.chi_opt = ({} for i in range(9))
        self.verbose = verbose

        if a_file.endswith('.json'):
            self.load_file(a_file)
        else:
            raise IOError("Invalid simulation parameters file '%s'." % a_file)

        # rewrite values in json info if new cosmo param passed
        for key, value in kwargs.iteritems():
            if self.verbose: print("Using new value for parameter '%s' = '%s'" % (key, value))
            self.cosmo[key] = value
        # other data attributes
        self.data = {}
        self._sim = None

    @property
    def sim(self):
        """ create C++ 'Sim_Param(self.file)' object
        may take a while to initialize CCL power spectra """
        if self._sim is None:
            self._sim = Sim_Param(self.file)
        return self._sim

    def info(self):
        # type: () -> str
        """ return string with basic info about run, to be used in plots """
        info = ''
        info += '%s:\n' % self.app
        info += '$N_m = %i$\n' % self.box_opt["mesh_num"]
        info += '$N_M = %i$\n' % self.box_opt["mesh_num_pwr"]
        info += '$N_p = %i^3$\n' % self.box_opt["par_num"]
        info += '$L = %i$ Mpc/h\n' % self.box_opt["box_size"]
        if self.app == 'AA': info += r'$\nu = %.1f$ (Mpc/h)$^2$' % self.app_opt["viscosity"]
        if self.app == 'FP_pp': info += r'$r_s = %.1f$' % self.app_opt["cut_radius"]
        if self.app == 'CHI':
            info += r'$\phi_s = %.1e$' % self.chi_opt["phi"]
            info += '  (lin)' if self.chi_opt["linear"] else '  (nl) '
        return info

    def info_tr(self):
        # type: () -> str
        """ return 'self.info()' as one-line string """
        info = self.info()
        return info.replace('\n', '  ')

    def info_supp(self):
        # type: () -> str
        """ return string with run resolution, viscozity and cut-off radius, to be used in plots """
        info = '%s: $res = %.2f$' % (self.app, float(self.box_opt["mesh_num"]) / self.box_opt["box_size"])
        if self.app == 'AA': info += r', $\nu = %.1f$' % self.app_opt["viscosity"]
        if self.app == 'FP_pp': info += r', $r_s = %.1f$' % self.app_opt["cut_radius"]
        return info

    def check_new_keys(self, key):
        if key not in self.results and key != 'files':
            self.results[key] = False

    def load_file(self, a_file):
        # type: (str) -> None
        """ load information in 'a_file', create results directory """
        with open(a_file) as data_file:
            data = json.loads(data_file.read())

        for key in data:
            setattr(self, key, data[key])

        if self.results is None:
            self.results = {}
        
        for key in sum(RESULTS_ALL.values(), []):
            self.check_new_keys(key)

        self.file = a_file
        self.dir = a_file.replace(a_file.split("/")[-1], '')
        self.res_dir = self.dir + 'results/'
        create_dir(self.res_dir)

    def rerun(self, rerun, key, skip, zs):
        """ steps to rerun or skip """
        # new analysis key?
        self.check_new_keys(key)

        # missing data
        if zs is None:
            if self.verbose: print("[Skipped]  (missing data)")
            return False
        # manually selected steps to rerun
        # check before skip-step in case of skip == 'all'
        elif _is_key_val(key, rerun):
            return True
        # manually selected steps to skip
        elif _is_key_val(key, skip):
            if self.verbose: print("[Skipped]")
            return False
        # step not done yet and not skipped
        elif not self.results[key]:
            return True
        # step already done
        else:
            if self.verbose: print("[Skipped]  (already done)")
            return False

    def done(self, key):
        self.results[key] = True
        with open(self.file) as data_file:
            data = json.loads(data_file.read())
        if data["results"] is None:
            data["results"] = {}
        data["results"][key] = True
        with open(self.file, 'w') as outfile:
            json.dump(data, outfile, indent=2)
        if self.verbose: print("[Done]")

class StackInfo(SimInfo):
    def __getitem__(self, key):
        return self.group_sim_infos[key]

    def __iter__(self):
        for x in self.group_sim_infos:
            yield x

    def __init__(self, group_sim_infos=None, stack_info_file=None, **kwargs):
        if group_sim_infos is not None:
            # use last SimInfo of SORTED list, i.e. the last run (if need to add
            # info, e.g. coorelation data)
            SimInfo.__init__(self, group_sim_infos[-1].file, **kwargs)
            self.load_group_sim_infos(group_sim_infos, **kwargs)
        elif stack_info_file is not None:
            SimInfo.__init__(self, stack_info_file, **kwargs)
        else:
            raise KeyError("Constructor 'StackInfo' called without arguments.")

        self.save() # need to save new cosmo param for C++ to load modified parameters
        self.data = {}
        for key in sum(RESULTS_ALL.values(), []):
            if key not in self.results:
                self.results[key] = False

    def load_group_sim_infos(self, group_sim_infos, **kwargs):
        self.group_sim_infos = group_sim_infos 

        # directory name
        self.dir = self.dir.replace(self.dir.split("/")[-2] + "/", "")
        self.dir += "STACK_%im_%ip_%iM_%ib" % (
            self.box_opt["mesh_num"], self.box_opt["par_num"],
            self.box_opt["mesh_num_pwr"], self.box_opt["box_size"])
        if self.app == "CHI":
            self.dir += "_%.1fn_%.0eY" % (self.chi_opt["n"], self.chi_opt["phi"])
            if self.chi_opt["linear"]:
                self.dir += "_lin"
        self.dir += "/"
        self.file = self.dir + 'stack_info.json'
        self.res_dir = self.dir + 'results/'

        create_dir(self.dir)
        create_dir(self.dir + "pwr_spec/")
        create_dir(self.dir + "pwr_diff/")
        create_dir(self.res_dir)

        self.seeds = [SI.run_opt["seed"] for SI in self]
        self.num_run = len(self.seeds)

        if os.path.isfile(self.file):  # there is already save StackInfo
            with open(self.file) as data_file:
                data = json.loads(data_file.read())
                self.results = data["results"]
                if self.num_run != data["num_run"]:
                    print("\tFound stack info but number of files does not seem right. Disregarding any saved data.")
                    self.results = {}
        else:  # save new StackInfo
            self.results = {}     

    def save(self):
        data = self.__dict__.copy()
        for key in ('ccl_cosmo', 'run_opt', 'out_dir', 'res_dir', 'data', 'sim', '_sim', 'group_sim_infos'):
            data.pop(key, None)
        # overriding
        with open(self.file, 'w') as outfile:
            json.dump(data, outfile, indent=2)


def compare(sim_info_1, sim_info_2):
    attributes = [a for a in dir(sim_info_1) if
                  not a.startswith('__') and a != "results" 
                  and a != "out_dir" and a != "dir" and a != "file"
                  and a != "run_opt" and a != "out_opt" and a != "res_dir"
                  and a != "data" and a != "sim"  and a != "_sim"
                  and not callable(getattr(sim_info_1, a))] #< last to avoid call to sim
    for att in attributes:
        att1 = getattr(sim_info_1, att)
        att2 = getattr(sim_info_2, att)

        if isinstance(att1, dict):
            for key in att1:
                if key not in att2 or att1[key] != att2[key]:
                    return False
        else:
            if att1 != att2:
                return False

    return True


def insert(a_sim_info, sep_files):
    for sep_sim_infos in sep_files:
        if compare(a_sim_info, sep_sim_infos[0]):
            sep_sim_infos.append(a_sim_info)
            return
    sep_files.append([a_sim_info])

class Results(object):
    """
    Class to handle all the results. Find, sort, show.
    """
    def __init__(self, out_dir=''):
        if out_dir != '':
            self.load(out_dir)
        else: self.sim_infos = []

    def sort(self):
        self.sim_infos = sorted(self.sim_infos, key=lambda si: si.app_opt["viscosity"])
        self.sim_infos = sorted(self.sim_infos, key=lambda si: si.box_opt["mesh_num_pwr"])
        self.sim_infos = sorted(self.sim_infos, key=lambda si: si.app)

    def load(self, out_dir):
        self.sim_infos = []
        files = get_files_in_traverse_dir(out_dir, 'stack_info.json')
        for a_file in files:
            self.sim_infos.append(StackInfo(stack_info_file=a_file))
        self.sort()

    def get_subfiles(self, Nm=0, NM=0, Np=0, L=0, nu=0, rs=0, phi=0, n=0, chi_lin=0, app='', app_not=''):
        subfiles = []
        for a_sim_info in self.sim_infos:
            check = bool (
                    (Nm == 0 or a_sim_info.box_opt["mesh_num"] == Nm) and
                    (NM == 0 or a_sim_info.box_opt["mesh_num_pwr"] == NM) and
                    (Np == 0 or a_sim_info.box_opt["par_num"]  == Np) and
                    (L == 0 or a_sim_info.box_opt["box_size"]  == L) and
                    (nu == 0 or a_sim_info.app_opt["viscosity"] == nu) and
                    (rs == 0 or a_sim_info.app_opt["cut_radius"] == rs) and
                    (app == '' or a_sim_info.app == app) and
                    (app_not == '' or a_sim_info.app != app_not)
                )
            if a_sim_info.chi_opt:
                check = check and (phi == 0 or a_sim_info.chi_opt["phi"] == phi)
                check = check and (n == 0 or a_sim_info.chi_opt["n"] == n)
                check = check and (chi_lin == 0 or a_sim_info.chi_opt["linear"] == chi_lin)
            if check:
                subfiles.append(a_sim_info)

        return subfiles

    def info(self, Nm=0, Np=0, L=0, nu=0, rs=0, phi=0, app=''):
        for a_sim_info in self.get_subfiles(Nm=Nm, Np=Np, L=L, nu=nu, rs=rs, phi=phi, app=app):
            info = a_sim_info.info_tr().replace('$', '')
            if hasattr(a_sim_info, 'num_run'):
                info += "\tnum runs = %i" % a_sim_info.num_run
            print(info)

    def show_folder(self, a_sim_info):
        subprocess.Popen(["xdg-open", a_sim_info.dir + 'results/'])

    def show_results(self, Nm=0, Np=0, L=0, nu=0, rs=0, phi=0, app='', plots=None):
        if plots is None:
            return
        elif plots == "all":
            return

        for a_sim_info in self.get_subfiles(Nm=Nm, Np=Np, L=L, nu=nu, rs=rs, phi=phi, app=app):
            results_dir = a_sim_info.dir + 'results/'
            for plot in plots:
                try:
                    filename=results_dir + plot + ".png"
                    display(Image(filename=filename))
                except IOError:
                    pass
    
    # def load_k_supp(self, a_sim_info):
    #     if not hasattr(a_sim_info, "supp"):
    #         zs, files = try_get_zs_files(a_sim_info, 'pwr_diff/')
    #         if zs is not None:
    #             a_sim_info.supp = load_k_supp(files)
    #             a_sim_info.a = [1./(z+1) for z in zs]

    # def plot_supp_compare(self, out_dir='/home/michal/Documents/GIT/Adhesion-Approximation/output/supp_comparison/',
    #                       Nm=0, Np=0, L=0, nu=0, rs=0, app='', scale=['small', 'medium', 'large'], show_k_lms=False, res=None):
    #     subfiles = self.get_subfiles(Nm=Nm, Np=Np, L=L, nu=nu, rs=rs, app=app)

    #     for a_sim_info in subfiles:
    #         self.load_k_supp(a_sim_info)
    #     for sc in scale:
    #         plot_supp(subfiles, out_dir+sc, suptitle=' on %s scales' % sc, save=True, show=True, scale=sc, show_k_lms=show_k_lms, res=res)

class Map(dict):
    """
    Example:
    m = Map({'first_name': 'Eduardo'}, last_name='Pool', age=24, sports=['Soccer'])
    """
    def __init__(self, *args, **kwargs):
        super(Map, self).__init__(*args, **kwargs)
        for arg in args:
            if isinstance(arg, dict):
                for k, v in arg.iteritems():
                    self[k] = v

        if kwargs:
            for k, v in kwargs.iteritems():
                self[k] = v

    def __getattr__(self, attr):
        return self.get(attr)

    def __setattr__(self, key, value):
        self.__setitem__(key, value)

    def __setitem__(self, key, value):
        super(Map, self).__setitem__(key, value)
        self.__dict__.update({key: value})

    def __delattr__(self, item):
        self.__delitem__(item)

    def __delitem__(self, key):
        super(Map, self).__delitem__(key)
        del self.__dict__[key]
