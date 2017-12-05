import os.path
import sys
import traceback
import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import brentq
import json

from . import *
from . import plot
from . import fastsim as fs
from .data import SimInfo, RESULTS_KEYS, load_k_supp_from_data

RESULTS_KEYS_FILES = ["corr_files", "pwr_spec_files",
                      "pwr_spec_diff_files"]
RESULTS_KEYS_STACK = RESULTS_KEYS + RESULTS_KEYS_FILES


class StackInfo(SimInfo):
    def __getitem__(self, key):
        return self.group_sim_infos[key]

    def __iter__(self):
        for x in self.group_sim_infos:
            yield x

    def __init__(self, group_sim_infos=None, stack_info_file=None):
        if group_sim_infos is not None:
            self.load_group_sim_infos(group_sim_infos)
        elif stack_info_file is not None:
            SimInfo.__init__(self, stack_info_file)
            self.sim = fs.Sim_Param(stack_info_file)
        else:
            raise KeyError("Constructor 'StackInfo' called without arguments.")

        self.data = {}
        for key in RESULTS_KEYS_STACK:
            if key not in self.results:
                self.results[key] = False

    def load_group_sim_infos(self, group_sim_infos):
        # use last SimInfo of SORTED list, i.e. the last run (if need to add
        # info, e.g. coorelation data)
        self.group_sim_infos = group_sim_infos
        SimInfo.__init__(self, self[-1].file)
        self.sim = fs.Sim_Param(self[-1].file)
        self.dir = self.dir.replace(self.dir.split("/")[-2] + "/", "")
        self.dir += "STACK_%im_%ip_%iM_%ib/" % (
            self.box_opt["mesh_num"], self.box_opt["par_num"],
            self.box_opt["mesh_num_pwr"], self.box_opt["box_size"])
        self.file = self.dir + 'stack_info.json'
        self.pwr_dir = self.dir + "pwr_spec/"
        self.pwr_diff_dir = self.dir + "pwr_diff/"
        self.corr_dir = self.dir + "corr_func/"
        self.res_dir = self.dir + 'results/'

        create_dir(self.dir)
        create_dir(self.pwr_dir)
        create_dir(self.pwr_diff_dir)
        create_dir(self.corr_dir)
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
        for key in ('ccl_cosmo', 'run_opt', 'out_dir', 'sim', 'group_sim_infos'):
            data.pop(key, None)
        # overriding
        with open(self.file, 'w') as outfile:
            json.dump(data, outfile, indent=2)


def compare(sim_info_1, sim_info_2):
    attributes = [a for a in dir(sim_info_1) if not a.startswith('__') and not callable(getattr(sim_info_1, a))
                  and a != "results" and a != "out_dir" and a != "dir" and a != "file"
                  and a != "res_dir" and a != "run_opt" and a != "ccl_cosmo"]
    for att in attributes:
        att1 = getattr(sim_info_1, att)
        att2 = getattr(sim_info_2, att)

        if isinstance(att1, dict):
            for key in att1:
                if att1[key] != att2[key]:
                    # print "Runs differ in %s[%s] = %s != %s\n" % (att, key,
                    # str(att1[key]), str(att2[key]))
                    return False
        else:
            if att1 != att2:
                # print "Runs differ in %s = %s != %s\n" % (att, str(att1),
                # str(att2))
                return False

    return True


def insert(a_sim_info, sep_files):
    for sep_sim_infos in sep_files:
        if compare(a_sim_info, sep_sim_infos[0]):
            sep_sim_infos.append(a_sim_info)
            return
    sep_files.append([a_sim_info])


def stack_files(stack_info, subdir, a_file):
    """ assume data in files are x = data[:, 0], y = data[:, 1]
    return tuple (zs, data_list) where each entry in 'data_list' corresponds to entry in 'zs' list
    data[:, 0] = x, data[:, 1] = mean(y), data[:, 2] = std(y)
    """
    # load everything
    all_data_k = None
    all_data_Pk = None
    for a_sim_info in stack_info:
        # load files for ONE run
        zs, files = try_get_zs_files(a_sim_info, subdir, a_file=a_file)

        # create lists, only for the first run (the rest are assumed to have the same output)
        if all_data_k is None:
            all_data_k = [[] for x in xrange(len(zs))]
        if all_data_Pk is None:
            all_data_Pk = [[] for x in xrange(len(zs))]

        # load k, Pk
        for i, data_file in enumerate(files):
            data = np.loadtxt(data_file)
            k, P_k = data[:, 0], data[:, 1]
            all_data_k[i].append(k.tolist())
            all_data_Pk[i].append(P_k.tolist())

    # chceck lengths of lists, delete excess (in case outputs of simulations encountered any errors)
    for i, data_k in enumerate(all_data_k):
        del_num = 0
        j = 0
        while True:
            k_row = [k_vec[j] for k_vec in  data_k]
            k_max = np.max(k_row)
            for ik, k in  enumerate(k_row):
                k_ = k
                while not np.isclose(k_, k_max, rtol=1.e-5, atol=1.e-5) and k_ < k_max :
                    # remove k, Pk if not the same, look for first close
                    if j == len(all_data_k[i][ik]):
                        break
                    del_num += 1
                    del all_data_k[i][ik][j]
                    del all_data_Pk[i][ik][j]
                    # look at the next k
                    k_ = all_data_k[i][ik][j]
            
            j += 1
            # check if j is length of ALL arrays
            for x in all_data_k[i]:
                if len(x) != j:
                    break
            else:
                break
        if del_num:
            print "\tSome values in files %s for z = %s did not match each other. Deleted %i excess values." % (subdir+a_file, str(zs[i]), del_num)

    # compute means, std
    data_list = []
    for i in xrange(len(all_data_k)):
        data_list.append([])

        data_list[i].append(np.mean(all_data_k[i], axis=0))
        data_list[i].append(np.mean(all_data_Pk[i], axis=0))
        data_list[i].append(np.std(all_data_Pk[i], axis=0))

    return zs, data_list

def save_pwr(zs, data_list, out_dir, app):
    for i in xrange(len(zs)):
        z_str = 'init.dat' if zs[i] == 'init' else 'z%.2f.dat' % zs[i]
        fname = out_dir + \
            "pwr_spec_par_%s_%s" % (app, z_str)
        header = ("This file contains power spectrum P(k) in units [(Mpc/h)^3] "
                  "depending on wavenumber k in units [h/Mpc] with standard deviation in units [h/Mpc].\n"
                  "k [h/Mpc]\tP(k) [(Mpc/h)^3]\tstd [(Mpc/h)^3]")
        np.savetxt(fname, np.transpose(
            data_list[i]), fmt='%.6e', header=header)


def save_pwr_diff(zs, data_list_diff, diff_type, out_dir, app):
    for i in xrange(len(zs)):
        z_str = 'init.dat' if zs[i] == 'init' else 'z%.2f.dat' % zs[i]
        fname = out_dir + \
            "pwr_spec_diff_%s_%s_%s" % (diff_type, app, z_str)

        header = ("This file contains relative difference between power spectrum P(k)\n"
                  "and lineary extrapolated ")
        if diff_type == "par":
            header += "power spectrum of initial particle position\n"
        elif diff_type == "input":
            header += "'hybrid' power spectrum\n"
        elif diff_type == "hybrid":
            header += "input power spectrum\n"

        header += ("depending on wavenumber k in units [h/Mpc] with standard deviation in units [h/Mpc].\n"
                   "k [h/Mpc]\tP(k) [(Mpc/h)^3]\tstd [(Mpc/h)^3]")
        np.savetxt(fname, np.transpose(
            data_list_diff[diff_type][i]), fmt='%.6e', header=header)


def load_pwr(stack_info, subdir='pwr_spec/', a_file='pwr_spec_par*.dat'):
    zs, files = try_get_zs_files(
        stack_info, subdir, a_file=a_file)
    data_list = []
    for a_file in files:
        data_list.append(np.transpose(np.loadtxt(a_file)))
    return zs, data_list

def load_pwr_par(stack_info):
    key = "pwr_spec_files"
    if not stack_info.results[key]:
        # stack all pwr_spec files, get mean of k and Pk, compute std of Pk
        zs, data_list = stack_files(stack_info, 'pwr_spec/', '*par*.dat')
        zs_in, data_in = stack_files(
            stack_info, 'pwr_spec/', '*init*.dat')
        zs.insert(0, zs_in[0])  # insert input data at the beggining
        data_list.insert(0, data_in[0])  # there should be only one init file
        save_pwr(zs, data_list, stack_info.pwr_dir, stack_info.app)
        stack_info.done(key)
        return zs, data_list
    else:
        return load_pwr(stack_info)

def load_pwr_diff(stack_info):
    key = "pwr_spec_diff_files"
    data_list_diff = {}
    for diff_type in ("par", "input", "hybrid"):
        if not stack_info.results[key]:
            # stack all pwr_spec_diff files, get mean of k and Pk, compute std
            # of Pk
            zs_diff, data_list_diff[diff_type] = stack_files(
                stack_info, 'pwr_diff/', '*%s*.dat' % diff_type)
            save_pwr_diff(zs_diff, data_list_diff, diff_type,
                          stack_info.pwr_diff_dir, stack_info.app)
        else:
            zs_diff, data_list_diff[diff_type] = load_pwr(
                stack_info, subdir='pwr_diff/', a_file='*%s*.dat' % diff_type)
    stack_info.done(key)
    return zs_diff, data_list_diff

def save_corr(zs, xi_list, out_dir, app, corr_type=''):
    for i in xrange(len(zs)):
        z_str = 'init.dat' if zs[i] == 'init' else 'z%.2f.dat' % zs[i]
        fname = out_dir + \
            "corr_func_scipy_qawf_%s_%s_%s" % (corr_type, app, z_str)
        header = ("This file contains correlation function depending on distance r in units [Mpc/h]."
                  "r [Mpc/h]\txsi(r)")
        np.savetxt(fname, np.transpose(xi_list[i]), fmt='%.6e', header=header)


def load_corr(stack_info, zs, Pk_list, Pk_nl_list):
    key = "corr_files"
    xi_all = {}
    if not stack_info.results[key]:
        # for each Pk in Pk_list, Pk_nl_list compute correlation function
        xi_all['par'] = [corr_func(stack_info.sim, Pk=Pk) for Pk in Pk_list]
        xi_all['emu'] = [corr_func(stack_info.sim, Pk=Pk) for Pk in Pk_nl_list]
        xi_all['lin'] = [corr_func(stack_info.sim, z=z) for z in zs]
        for corr_type in ('par', 'lin', 'emu'):
            zs_ = zs if corr_type != 'emu' else [z for z in zs if z <= 2.02 or z == 'init']
            save_corr(zs_, xi_all[corr_type], stack_info.corr_dir, stack_info.app, corr_type=corr_type)
        stack_info.done(key)
    else:
        for corr_type in ('par', 'lin', 'emu'):
            files = try_get_zs_files(
                stack_info, 'corr_func/', a_file='corr_func_scipy_qawf_%s_*.dat' % corr_type)[1]
            xi_all[corr_type] = []
            for a_file in files:
                xi_all[corr_type].append(np.transpose(np.loadtxt(a_file)))

    return xi_all


def get_z_from_A(cosmo, A):
    # 'f = 0' <=> A = D^2 (linear power grows as D^2)
    f = lambda a : A - fs.growth_factor(a, cosmo)**2
    a = brentq(f, 0, 1)
    return 1./a - 1


def get_Extrap_pk_nl(sim, z):
    """ return Extrap_Pk object for non-linear power spectra"""
    if z == 'init' or z < 0: z = 0.0
    elif z > 2.02: z = 2.02
    return fs.Extrap_Pk(fs.init_emu(sim, z), sim, 0, 10, 341, 351)


def get_ndarray(Data_Vec):
    xx = [x for x in Data_Vec[0]]
    yy = [y for y in Data_Vec[1]]
    return np.array([xx, yy])


def get_Data_vec(data):
    size = len(data[0])
    Data_Vec = fs.Data_d3(size)
    for i in range(size):
        Data_Vec[0][i] = data[0][i]
        Data_Vec[1][i] = data[1][i]
        Data_Vec[2][i] = data[2][i] 
    return Data_Vec

def pwr_spec_non_linear(sim, kk, z):
    """ return ndarray of power spectra for every k in kk (iterable) """
    Pk_nl = get_Extrap_pk_nl(sim, z)
    return np.array([Pk_nl(k) for k in kk])

def pwr_spec_linear(sim, kk, z):
    D = fs.growth_factor(1./(1.+z), sim.cosmo)
    return np.array([D*D*fs.lin_pow_spec(k, sim.cosmo) for k in kk])

def pwr_spec_hybrid(sim, kk, z, A):
    if z > 2.02:
        # A = 0, TBD: extrapolate emulator
        return pwr_spec_linear(sim, kk, z)
    else:
        return (1-A)*pwr_spec_linear(sim, kk, z) + A*pwr_spec_non_linear(sim, kk, z)


def get_hybrid_pwr_spec(sim, data, k_nyquist_par, z=1):
    """ given data [k, Pk, std] and initial guess of redshift:
    return fit of hybrid power spectrum: (1-A)*P_lin(k) + A*P_nl(k) """
    # define functions which will be used in fitting
    pk_nl_func = lambda kk, z: pwr_spec_non_linear(sim, kk, z)
    pk_lin_func = lambda kk, z: pwr_spec_linear(sim, kk, z)
    pk_hyb_func = lambda kk, z, A: pwr_spec_hybrid(sim, kk, z, A)

    # extract data
    kk, Pk, std = data
    
    # get proper slice of data -- last decade befor half of particle nyquist
    idx = (np.abs(kk-0.5*k_nyquist_par)).argmin()
    idx = slice(idx - sim.out_opt.bins_per_decade, idx)

    # fit data
    popt, pcov = curve_fit(pk_hyb_func, kk[idx], Pk[idx], sigma=std[idx], p0=(z, 0.3))

    # get hybrid Extrap
    if popt[0] > 2.02:
        Pk_NL = fs.Extrap_Pk(get_Data_vec(data), sim)
    else:
        Pk_NL = fs.Extrap_Pk_Nl(get_Data_vec(data), sim, popt[1], popt[0])

    # return all info in dict
    return {"Pk_NL" : Pk_NL, "popt" : popt, "pcov" : pcov}


def corr_func(sim, Pk=None, z=None):
    corr = fs.Data_d2()
    if Pk is not None: # compute correlation function from given continuous power spectrum
        fs.gen_corr_func_binned_gsl_qawf(sim, Pk, corr)
    elif z is not None: # compute linear correlation function at given redshift
        a = 1./(1.+z) if z != 'init' else 1.0
        fs.gen_corr_func_binned_gsl_qawf_lin(sim, a, corr)
    else:
        raise KeyError("Function 'corr_func' called without arguments.")
    return get_ndarray(corr)


def stack_group(group_sim_infos=None, stack_info_file=None, plot_all=True):
    # load & save all info about stack
    stack_info = StackInfo(group_sim_infos=group_sim_infos, stack_info_file=stack_info_file)

    print "\tLoad & save data to stack..."
    # load data from stacked files if available, otherwise do stacking over all runs
    zs, data_list = load_pwr_par(stack_info)
    zs_diff, data_list_diff = load_pwr_diff(stack_info)

    # for each entry in data_list compute its continuous version
    # get non-linear power spectra for z in emulator range, z == 'init' transformed into z = 0.0
    Pk_hyb_lis = [get_hybrid_pwr_spec(stack_info.sim, data, stack_info.k_nyquist["particle"]) for data in data_list]
    #Pk_list = [fs.Extrap_Pk(get_Data_vec(data), stack_info.sim) for data in data_list]
    Pk_list = [x["Pk_NL"] for x in Pk_hyb_lis]
    Pk_nl_list = [get_Extrap_pk_nl(stack_info.sim, z) for z in zs if z <= 2.02 or z == 'init']

    # compute correlation function for simulated, linear and non-linear power spectra
    # return dictionary, keys = 'par', 'lin', 'emu'
    xi_all = load_corr(stack_info, zs, Pk_list, Pk_nl_list)

    # save all data to stack_info for later analysis
    stack_info.data["zs"] = zs
    stack_info.data["data_list"] = data_list
    stack_info.data["zs_diff"] = zs_diff
    stack_info.data["data_list_diff"] = data_list_diff
    stack_info.data["Pk_hyb_lis"] = Pk_hyb_lis
    stack_info.data["Pk_list"] = Pk_list
    stack_info.data["Pk_nl_list"] = Pk_nl_list
    stack_info.data["xi_all"] = xi_all

    if plot_all:
        print '\tPlotting power spectrum...'
        plot.plot_pwr_spec_stacked(data_list, zs, stack_info, Pk_list, Pk_nl_list)

        print '\tPlotting power spectrum difference...'
        for diff_type, data in data_list_diff.iteritems():
            plot.plot_pwr_spec_diff_from_data(data, zs_diff, stack_info, ext_title=diff_type)

        print '\tPlotting power spectrum suppression...'
        a = [1./(z+1) for z in zs_diff]
        supp_lms, supp_std_lms, k_lms = load_k_supp_from_data(data_list_diff["par"], stack_info.k_nyquist["particle"])
        plot.plot_supp_lms(supp_lms, a, stack_info, k_lms, supp_std_lms)

        print '\tPlotting power spectrum suppression (map)...'
        # WARNING! data in data_list_diff are potentially modified -- aligning values to the same counts
        plot.plot_pwr_spec_diff_map_from_data(data_list_diff["par"], zs_diff, stack_info, ext_title="par")

        print '\tPlotting correlation function...'
        plot.plot_corr_func_from_data(xi_all, zs, stack_info)

    return None
    return stack_info


def stack_all(in_dir='/home/vrastil/Documents/GIT/Adhesion-Approximation/output/', plot_all=True):
    # get all runs
    files = get_files_in_traverse_dir(in_dir, 'sim_param.json')
    sim_infos = [SimInfo(a_file[0]) for a_file in files]
    del files

    # separate files according to run parameters
    sep_files = []
    for a_sim_info in sim_infos:
        insert(a_sim_info, sep_files)
    del sim_infos

    # print info about separated files
    num_all_runs = num_all_sep_runs = num_sep_runs = 0
    for sep_sim_infos in sep_files:
        num_all_runs += len(sep_sim_infos)
        if len(sep_sim_infos) > 1:
            num_sep_runs += 1
            num_all_sep_runs += len(sep_sim_infos)

    # remove 1-length sep_sim_infos
    sep_files[:] = [x for x in sep_files if len(x) != 1]

    # sort sim_infos
    for sep_sim_infos in sep_files:
        sep_sim_infos.sort(key=lambda x: x.dir)

    print "There are in total %i different runs, from which %i share the same parameters, constituting %i group(s) eligible for stacking:" % (
        num_all_runs, num_all_sep_runs, num_sep_runs)
    for sep_sim_infos in sep_files:
        print "\n" + sep_sim_infos[0].info_tr()
        for i, a_sim_info in enumerate(sep_sim_infos):
            if i == 10:
                print "\t...and %i more runs" % (len(sep_sim_infos) - 10)
                break
            print "\t" + a_sim_info.dir

    stack_infos = []
    for sep_sim_infos in sep_files:
        print "\nStacking group %s" % sep_sim_infos[0].info_tr()
        try:
            stack_infos.append(stack_group(sep_sim_infos, plot_all=plot_all))
        except:
            print "=" * 110
            traceback.print_exc(file=sys.stdout)
            print "=" * 110
    return stack_infos
