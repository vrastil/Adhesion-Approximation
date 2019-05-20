"""
'struct.py' module defines all necessary data containers for working with results

 - Sim_Info : basic class storing all information about one particular simulation run
 - Stack_Info : subclass of Sim_Info, collects and stores information about all runs
   with the same parameters
 - Results : collects results, compare approximation, show plots
"""

from __future__ import print_function

import json
from bson.binary import Binary
from bson.son import SON
from bson.objectid import ObjectId
import pickle
import fnmatch
import datetime
import numpy as np
from . import utils as ut
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
               "pwr_diff_files_i", "pwr_diff_files_h", "pwr_spec_chi_files",
               "vel_pwr_spec_files", "vel_pwr_diff_spec_files"]
    }

RESULTS_DIRS = {
    'pwr_spec' : 'pwr_spec',
    'pwr_spec_chi' : 'pwr_spec',
    'vel_pwr_spec' : 'vel_pwr_spec',
    'pwr_slope' : 'pwr_spec',
    'pwr_diff' : 'pwr_diff',
    'pwr_diff_h' : 'pwr_diff',
    'pwr_diff_i' : 'pwr_diff',
    'vel_pwr_diff' : 'vel_pwr_diff',
    'chi_pwr_diff' : 'pwr_spec',
    'pwr_spec_supp' : 'pwr_diff',
    'pwr_spec_supp_map' : 'pwr_diff',
    'vel_pwr_spec_supp' : 'vel_pwr_diff',
    'chi_pwr_spec_supp' : 'pwr_spec',
    'chi_pwr_spec_supp_map' : 'pwr_spec',
    'corr_func' : 'pwr_spec',
    'bao' : 'pwr_spec',
    'sigma_R' : 'pwr_spec',
    'dens_hist' : 'rho_bin',
    'par_slice' : 'par_cut',
    'par_ani' : 'par_cut',
    'dens_slice' : 'rho_map',
    'dens_ani' : 'rho_map',
    'eff_time' : 'pwr_spec',
    'pwr_spec_files' : 'pwr_spec',
    'pwr_spec_chi_files' : 'pwr_spec',
    'pwr_diff_files' : 'pwr_diff',
    'pwr_diff_files_h' : 'pwr_diff',
    'pwr_diff_files_i' : 'pwr_diff',
    "vel_pwr_spec_files" : "vel_pwr_spec",
    "vel_pwr_diff_spec_files" : "vel_pwr_diff"
}

# save all initialized infos here for later retrieve, key by ObjectId
ALL_STACK_INFOS = {}

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

class SimInfo(object):
    """ basic class storing all information about one particular simulation run """

    def __init__(self, db, doc_id, collection='data', verbose=True, **kwargs):
        """ load information from database, replace parameters in 'self.cosmo'
        by any additional parameters in 'kwargs' (optional) """

        if isinstance(doc_id, dict) and '_id' in doc_id:
            pass # default option
        elif isinstance(doc_id, ObjectId):
            doc_id = {'_id' : doc_id}
        else:
            raise ValueError("Wrong type of argument: '%s'" % type(doc_id))
        
        # load data from the database, omit files
        doc = db[collection].find_one(doc_id, {'data' : 0})

        # attributes
        self.doc_id = doc_id
        self.db = db
        self.collection = db[collection]
        self.app = doc['app']
        self.cosmo = doc['cosmo']
        self.k_nyquist = doc['k_nyquist']
        self.integ_opt = doc['integ_opt']
        self.box_opt = doc['box_opt']
        self.app_opt = doc['app_opt']
        self.run_opt = doc['run_opt']
        self.out_opt = doc['out_opt']
        self.results = doc.get('results', {})
        self.chi_opt = doc.get('chi_opt', None)
        self.verbose = verbose
        self.file = doc['out_opt']['file']
        self.dir = doc['out_opt']['dir']
        self.res_dir = self.dir + 'results/'
        ut.create_dir(self.res_dir)
        
        # check if we added new result keys
        for key in sum(RESULTS_ALL.values(), []):
            self.check_new_keys(key)

        # rewrite values in database info if new cosmo param passed
        for key, value in kwargs.items():
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
            self._sim = Sim_Param(str(self.file)) # !!! VERY IMPORTANT -- do not use unicode string
        return self._sim

    def info(self, math_mode=False):
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
            info += '   $n = %.1f$' % self.chi_opt["n"]
            info += '  (lin)' if self.chi_opt["linear"] else '  (nl) '
        
        # for use in display(Math(...))
        if math_mode:
            info = info.replace('$', r'\ ')
        
        return info

    def info_tr(self, math_mode=False):
        # type: () -> str
        """ return 'self.info()' as one-line string """
        info = self.info(math_mode=math_mode)
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

    def rerun(self, rerun, key, skip, zs):
        """ steps to rerun or skip """
        # new analysis key?
        self.check_new_keys(key)

        # missing data
        if zs is None:
            if self.verbose: ut.print_skipped_miss()
            return False
        # manually selected steps to rerun
        # check before skip-step in case of skip == 'all'
        elif _is_key_val(key, rerun):
            return True
        # manually selected steps to skip
        elif _is_key_val(key, skip):
            if self.verbose: ut.print_skipped()
            return False
        # step not done yet and not skipped
        elif not self.results[key]:
            return True
        # step already done
        else:
            if self.verbose: ut.print_skipped(done=True)
            return False

    def update_db(self, upd_doc):
        self.collection.find_one_and_update(
            self.doc_id, # find self
            {'$set': upd_doc} # update
        )

    def done(self, key):
        self.results[key] = True
        self.update_db({'results': self.results})
        if self.verbose: ut.print_done()

    def save_data_to_db(self, keys=None):
        # if not specified, save everything
        if keys is None:
            keys = self.data.keys()
        
        err_keys = []
        # save in binary format -- allow saving numpy arrays (not C++ objects)
        for key in keys:
            try:
                binary_data = pickle.dumps(self.data[key])
            except TypeError: # can't pickle object
                err_keys.append(key)
            else:
                self.update_db({'data.processed.%s' % key : Binary(binary_data)})

        return err_keys

    def get_data_from_db(self, keys=None):
        # get doc with info about processed data
        doc = self.collection.find_one(self.doc_id, {'data.processed'})

        # check we have loaded any data
        if doc is None or 'data' not in doc or 'processed' not in doc['data']:
            return
        
        # if not specified, load everything
        if keys is None:
            keys = doc['data']['processed'].keys()

        # load from binary format
        for key in keys:
            binary_data = doc['data']['processed'][key]
            self.data[key] = pickle.loads(binary_data)     

    def get_zs_data(self, key, patterns, doc_id=None):
        # if doc_id is specified, do not load self (for StackInfo)
        if doc_id is None:
            doc_id = self.doc_id

        # get all data for given key
        subdir = RESULTS_DIRS[key]
        all_data = self.collection.find_one(doc_id, {'data.files.%s' % subdir}) # find data
        all_data = all_data['data']['files'].get(subdir, []) # get data or empty list if no data available

        # save zs, data
        zs, data = [], []

        for record in all_data: # go through all data
            match = False
            for pattern in patterns.split(): # if multiple patterns given
                match = match | fnmatch.fnmatch(record['file'], pattern) # pattern matching
            if match:
                zs.append(record['z'])
                binary_data = record['data']
                data.append(pickle.loads(binary_data))
        
        # check for empty list (no file matched)
        if zs:
            zs, data = ut.sort_lists(zs, data)
            return list(zs), list(data)
        else:
            return None, None

    def save_zs_data(self, key, zs, data_list, fname):
        data_dir = RESULTS_DIRS[key]
        db_key = "data.files.%s" % data_dir

        for z, data in zip(zs, data_list):
            z_str = 'init.dat' if z == 'init' else 'z%.2f.dat' % z
            fname_ = self.dir + RESULTS_DIRS[key] + '/' + fname + "_%s_%s" % (self.app, z_str)

            # load data and save them in binary
            binary_data = pickle.dumps(np.array(data))
            upd_doc = SON({
                'file' : fname_,
                'z' : z,
                'data' : Binary(binary_data)
            })
        
            # save into the database
            self.collection.find_one_and_update(
                self.doc_id, # find self
                {'$addToSet': {db_key : upd_doc}} # update
            )

class StackInfo(SimInfo):
    def __getitem__(self, key):
        return self.sim_ids[key]

    def __iter__(self):
        for x in self.sim_ids:
            yield x

    def __init__(self, db, sep_id, collection='data', **kwargs):
        # we can pass single sep_id of StackInfo or group of ids of SimInfo
        if isinstance(sep_id, list):
            SimInfo.__init__(self, db, sep_id[0], collection=collection, **kwargs)
            self.get_dir_name()
            self.sim_ids = sep_id
            self.num_run = len(self.sim_ids)
            self.find_data_in_db()
        else:
            SimInfo.__init__(self, db, sep_id, collection=collection, **kwargs)
            # load additional variables
            self.get_data_from_db()
            doc = self.collection.find_one(self.doc_id, {'database' : 1, 'results' : 1})
            self.sim_ids = doc['database']['sim_ids']
            self.num_run = len(self.sim_ids)
            self.results = doc['results']
        
        # need to save new cosmo param for C++ to load modified parameters
        self.save()

        # try to load data if we already save something during this run of program
        self.load_temp()

    def save_temp(self):
        glob_key = self.doc_id['_id']
        global ALL_STACK_INFOS

        keys = self.data.keys()

        if glob_key not in ALL_STACK_INFOS:
            ALL_STACK_INFOS[glob_key] = {}

        for key in keys:
            if key not in ALL_STACK_INFOS[glob_key]:
                ALL_STACK_INFOS[glob_key][key] = self.data[key]
        # need to use the sam SimParam as in Extra_Pk, etc.
        ALL_STACK_INFOS[glob_key]['_sim'] = self._sim

    def load_temp(self):
        glob_key = self.doc_id['_id']
        global ALL_STACK_INFOS

        if glob_key in ALL_STACK_INFOS:
            keys = ALL_STACK_INFOS[glob_key].keys()
            for key in keys:
                self.data[key] = ALL_STACK_INFOS[glob_key][key]
            # need to use the sam SimParam as in Extra_Pk, etc.
            self._sim = ALL_STACK_INFOS[glob_key]['_sim']

    def find_data_in_db(self):
        # get document with run information
        doc = self.collection.find_one(self.sim_ids[0], {
            '_id' : 0,
            'app_opt' : 0,
            'data' : 0,
            'database' : 0,
            'out_opt' : 0,
            'results' : 0,
            'run_opt' : 0,
        })

        doc['type'] = 'stack_info'
        doc_in_db = self.collection.find_one(doc, {'database' : 1, 'results' : 1})

        if doc_in_db is not None:
            self.doc_id = {'_id' : doc_in_db['_id']}
            if doc_in_db['database']['sim_ids'] == self.sim_ids:
                # same files, load data
                self.get_data_from_db()
                self.results = doc_in_db['results']

            else:
                # wrong number of files
                doc_in_db['database'] = {
                    'sim_ids' : self.sim_ids,
                    'timestamp' : str(datetime.datetime.utcnow())
                }
                doc_in_db['results'] = self.results
                self.update_db(doc_in_db)
        else:
            doc['database'] = {
                'sim_ids' : self.sim_ids,
                'timestamp' : str(datetime.datetime.utcnow())
            }
            doc['run_opt'] = self.run_opt
            doc['app_opt'] = self.app_opt
            doc['out_opt'] = self.out_opt
            doc['out_opt']['file'] = self.file
            doc['out_opt']['dir'] = self.dir
            doc['results'] = self.results
            doc['chi_opt'] = self.chi_opt
            self.doc_id = {'_id' : self.collection.insert(doc)}

    def get_dir_name(self):
        """get directory name for stack info from sim info, create structure"""
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

        ut.create_dir(self.dir)
        ut.create_dir(self.dir + "pwr_spec/")
        ut.create_dir(self.dir + "pwr_diff/")
        ut.create_dir(self.dir + "vel_pwr_spec/")
        ut.create_dir(self.dir + "vel_pwr_diff/")
        ut.create_dir(self.res_dir)         

    def save(self):
        # do not modify self, do a copy
        data = self.__dict__.copy()

        # get rid of unnecessary data
        for key in ('ccl_cosmo', 'run_opt', 'out_dir', 'res_dir', 'verbose', 'results',
            'data', 'sim', '_sim', 'collection', 'db', 'file', 'dir', 'doc_id', 'sim_ids'):
            data.pop(key, None)
        if data['chi_opt'] is None:
            data.pop('chi_opt')

        # convert ObjectId into string and save
        data['database'] = {
            "timestamp" : str(datetime.datetime.utcnow()),
            "_id" : str(self.doc_id['_id']),
            "sim_ids" : []
        }
        
        for doc_id in self.sim_ids:
            data['database']['sim_ids'].append(str(doc_id['_id']))

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
