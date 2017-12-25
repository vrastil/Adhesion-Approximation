"""
'struct.py' module defines all necessary data containers for working with results

 - Sim_Info : basic class storing all information about one particular simulation run
 - Stack_Info : subclass of Sim_Info, collects and stores information about all runs
   with the same parameters
"""

import json
import os
from .fastsim import Sim_Param

RESULTS_KEYS = ["pwr_spec", "pwr_diff", "pwr_diff_i", "pwr_diff_h",
                "pwr_spec_supp", "pwr_spec_supp_map",
                "dens_hist", "vel_pwr_spec", "vel_pwr_diff", "vel_pwr_spec_supp",
                "par_slice", "par_ani", "dens_slice", "dens_ani", "corr_func"]

RESULTS_KEYS_FILES = ["corr_files", "pwr_spec_files", "pwr_diff_files",
                      "pwr_diff_files_i", "pwr_diff_files_h"]

RESULTS_KEYS_STACK = RESULTS_KEYS + RESULTS_KEYS_FILES


def create_dir(out_dir):
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

class SimInfo(object):
    """ basic class storing all information about one particular simulation run """
    def __init__(self, a_file, **kwargs):
        # type: (str, dict) -> None
        """ load information in 'a_file', replace parameters in 'self.cosmo'
        by any additional parameters in 'kwargs' (optional) """
        # attributes to load from json file: string & dictionaries
        self.app = ""
        self.app_opt = self.box_opt = self.cosmo = self.integ_opt = self.k_nyquist = self.out_opt = self.results = self.run_opt = {}

        if a_file.endswith('.json'):
            self.load_file(a_file)
        else:
            raise IOError("Invalid simulation parameters file '%s'." % a_file)

        # rewrite values in json info if new cosmo param passed
        for key, value in kwargs.iteritems():
            print "Using new value for parameter '%s' = '%s'" % (key, value)
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

    def load_file(self, a_file):
        # type: (str) -> None
        """ load information in 'a_file', create results directory """
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
        self.results[key] = True
        with open(self.file) as data_file:
            data = json.loads(data_file.read())
        if data["results"] is None:
            data["results"] = {}
        data["results"][key] = True
        with open(self.file, 'w') as outfile:
            json.dump(data, outfile, indent=2)

class StackInfo(SimInfo):
    def __getitem__(self, key):
        return self.group_sim_infos[key]

    def __iter__(self):
        for x in self.group_sim_infos:
            yield x

    def __init__(self, group_sim_infos=None, stack_info_file=None, **kwargs):
        if group_sim_infos is not None:
            self.load_group_sim_infos(group_sim_infos, **kwargs)
        elif stack_info_file is not None:
            SimInfo.__init__(self, stack_info_file, **kwargs)
        else:
            raise KeyError("Constructor 'StackInfo' called without arguments.")

        self.data = {}
        for key in RESULTS_KEYS_STACK:
            if key not in self.results:
                self.results[key] = False

    def load_group_sim_infos(self, group_sim_infos, **kwargs):
        # use last SimInfo of SORTED list, i.e. the last run (if need to add
        # info, e.g. coorelation data)
        self.group_sim_infos = group_sim_infos
        SimInfo.__init__(self, self[-1].file, **kwargs)
        self.dir = self.dir.replace(self.dir.split("/")[-2] + "/", "")
        self.dir += "STACK_%im_%ip_%iM_%ib/" % (
            self.box_opt["mesh_num"], self.box_opt["par_num"],
            self.box_opt["mesh_num_pwr"], self.box_opt["box_size"])
        self.file = self.dir + 'stack_info.json'
        self.res_dir = self.dir + 'results/'

        create_dir(self.dir)
        create_dir(self.dir + "pwr_spec/")
        create_dir(self.dir + "pwr_diff/")
        create_dir(self.res_dir)

        self.seeds = [SI.run_opt["seed"] for SI in self]

        if os.path.isfile(self.file):  # there is already save StackInfo
            with open(self.file) as data_file:
                data = json.loads(data_file.read())
                self.results = data["results"]
                if len(self.seeds) != len(data["seeds"]):
                    print "\tFound stack info but number of files does not seem right. Disregarding any saved data."
                    self.results = {}
                    self.save()
        else:  # save new StackInfo
            self.results = {}
            self.save()     

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
                if att1[key] != att2[key]:
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