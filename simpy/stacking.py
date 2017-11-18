from datetime import datetime
import os.path
import sys
import gc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import curve_fit
from scipy.integrate import quad
import json

from . import *
from . import plot
from .data import SimInfo
from .data import RESULTS_KEYS

RESULTS_KEYS_FILES = ["corr_files", "pwr_spec_files",
                      "pwr_spec_extra_files", "pwr_spec_diff_files"]
RESULTS_KEYS_STACK = RESULTS_KEYS + RESULTS_KEYS_FILES


class StackInfo(SimInfo):
    def __init__(self, group_sim_infos):
        # use last SimInfo of SORTED list, i.e. the last run (if need to add
        # info, e.g. coorelation data)
        SimInfo.__init__(self, group_sim_infos[-1].file)
        self.last = group_sim_infos[-1]
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

        self.seeds = [SI.run_opt["seed"] for SI in group_sim_infos]

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

        for key in RESULTS_KEYS_STACK:
            if key not in self.results:
                self.results[key] = False

    def save(self):
        data = self.__dict__.copy()
        for key in ('ccl_cosmo', 'run_opt', 'out_dir', 'last'):
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


def stack_files(group_sim_infos, subdir, a_file):
    """ assume data in files are x = data[:, 0], y = data[:, 1]
    return tuple (zs, data_list) where each entry in 'data_list' corresponds to entry in 'zs' list
    data[:, 0] = x, data[:, 1] = mean(y), data[:, 2] = std(y)
    """
    # load everything
    all_data_k = None
    all_data_Pk = None
    for a_sim_info in group_sim_infos:
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
        print "\n\tz = ", zs[i]
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


class Interp_obj(object):
    def __init__(self, data):
        self.interpolate = PchipInterpolator(data[0], data[1])


class Extrap_Pk(Interp_obj):
    def __init__(self, data, cosmo, ccl_cosmo=None, k_max=None):
        Interp_obj.__init__(self, data)

        k = data[0]
        Pk = data[1]
        dPk = data[2]
        # number of points per decade
        decade = int(len(k) / np.log10(k[-1] / k[0]))

        self.cosmo = cosmo
        self.ccl_cosmo = ccl_cosmo

        def f(k_, A_): return A_ * \
            plot.pwr_spec(k_, self.cosmo, self.ccl_cosmo)
        popt, pcov = curve_fit(
            f, k[0:decade], Pk[0:decade], sigma=dPk[0:decade])  # fit over a decade
        self.A_lower = popt
        self.k_lower = k[decade / 2]
        ind = (np.abs(k - 0.5 * k_max)).argmin() if k_max is not None else -1
        popt, pcov = curve_fit(f, k[ind - decade / 2:ind], Pk[ind - decade / 2:ind],
                               sigma=dPk[ind - decade / 2:ind])  # fit over a half of a decade
        self.A_upper = popt
        self.k_upper = k[ind]

    def eval(self, k):
        if k < self.k_lower:
            return self.A_lower * plot.pwr_spec(k, self.cosmo, self.ccl_cosmo)
        elif k < self.k_upper:
            return self.interpolate(k)
        else:
            return self.A_upper * plot.pwr_spec(k, self.cosmo, self.ccl_cosmo)

    def eval_list(self, k_list):
        Pk_list = []
        for k in k_list:
            Pk_list.append(self.eval(k))
        return Pk_list


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


def save_extra_pwr(zs, data_list, Pk_list, out_dir, app):
    for i in xrange(len(zs)):
        z_str = 'init.dat' if zs[i] == 'init' else 'z%.2f.dat' % zs[i]
        fname = out_dir + \
            "pwr_spec_extrap_%s_%s" % (app, z_str)
        header = ("This file contains extra/intra-polated power spectrum P(k) in units [(Mpc/h)^3] "
                  "depending on wavenumber k in units [h/Mpc].\n"
                  "k [h/Mpc]\tP(k) [(Mpc/h)^3]")
        np.savetxt(fname, np.transpose([
            data_list[i][0], Pk_list[i].eval_list(data_list[i][0])]), fmt='%.6e', header=header)


def save_corr(zs, data_list, xi_list, out_dir, app):
    for i in xrange(len(zs)):
        z_str = 'init.dat' if zs[i] == 'init' else 'z%.2f.dat' % zs[i]
        fname = out_dir + \
            "corr_func_scipy_qawf_par_%s_%s" % (app, z_str)
        header = ("This file contains correlation function depending on distance r in units [Mpc/h]."
                  "r [Mpc/h]\txsi(r)")
        np.savetxt(fname, np.transpose(xi_list[i]), fmt='%.6e', header=header)


def load_corr(stack_info):
    zs, files = try_get_zs_files(
        stack_info, 'corr_func/', a_file='corr_func_scipy_qawf_par_*.dat')
    xi_list = []
    for a_file in files:
        xi_list.append(np.transpose(np.loadtxt(a_file)))
    return xi_list


def corr_func_r(r, Pk):
    def f(k): return 1. / (2. * np.pi**2) * k / r * Pk.eval(k)
    sys.stdout.write('\tr = %f\r' % r)
    sys.stdout.flush()
    xi, err = quad(f, 0, np.inf, weight='sin', wvar=r)
    return xi


def corr_func(Pk, r_min=1, r_max=200, num=50):
    r_range = np.linspace(r_min, r_max, num=num)
    print "\tComputing correlation function..."
    corr = [corr_func_r(r, Pk) for r in r_range]
    print "\n\tDone!"
    return [r_range, corr]


def stack_group(group_sim_infos):
    # load & save all info about stack
    stack_info = StackInfo(group_sim_infos)
    print "\tLoad & save data to stack..."

    key = "pwr_spec_files"
    if not stack_info.results[key]:
        # stack all pwr_spec files, get mean of k and Pk, compute std of Pk
        zs, data_list = stack_files(group_sim_infos, 'pwr_spec/', '*par*.dat')
        zs_in, data_in = stack_files(
            group_sim_infos, 'pwr_spec/', '*init*.dat')
        zs.insert(0, zs_in[0])  # insert input data at the beggining
        data_list.insert(0, data_in[0])  # there should be only one init file
        save_pwr(zs, data_list, stack_info.pwr_dir, stack_info.app)
        stack_info.done(key)
    else:
        zs, data_list = load_pwr(stack_info)

    key = "pwr_spec_diff_files"
    data_list_diff = {}
    for diff_type in ("par", "input", "hybrid"):
        if not stack_info.results[key]:
            # stack all pwr_spec_diff files, get mean of k and Pk, compute std
            # of Pk
            zs_diff, data_list_diff[diff_type] = stack_files(
                group_sim_infos, 'pwr_diff/', '*%s*.dat' % diff_type)
            save_pwr_diff(zs_diff, data_list_diff, diff_type,
                          stack_info.pwr_diff_dir, stack_info.app)
        else:
            zs_diff, data_list_diff[diff_type] = load_pwr(
                stack_info, subdir='pwr_diff/', a_file='*%s*.dat' % diff_type)
    stack_info.done(key)

    # for each entry in data_list compute its extra/inter-polated version
    Pk_list = [Extrap_Pk(data, stack_info.cosmo, stack_info.ccl_cosmo,
                         stack_info.k_nyquist["particle"]) for data in data_list]
    key = "pwr_spec_extra_files"
    if not stack_info.results[key]:
        save_extra_pwr(zs, data_list, Pk_list,
                       stack_info.pwr_dir, stack_info.app)
        stack_info.done(key)

    # key = "corr_files"
    # if not stack_info.results[key]:
    #     # for each Pk in Pk_list compute correlation function
    #     xi_list = [corr_func(Pk) for Pk in Pk_list]
    #     save_corr(zs, data_list, xi_list, stack_info.corr_dir, stack_info.app)
    #     stack_info.done(key)
    # else:
    #     xi_list = load_corr(stack_info)

    print '\tPlotting power spectrum...'
    plot.plot_pwr_spec_stacked(
        data_list, zs, stack_info, Pk_list)

    print '\tPlotting power spectrum difference...'
    for diff_type in ("par", "input", "hybrid"):
        plot.plot_pwr_spec_diff_from_data(
            data_list_diff[diff_type], zs_diff, stack_info, ext_title=diff_type)

    # print '\tPlotting correlation function...'
    # zs_emu, files_emu = try_get_zs_files(
    #     stack_info.last, 'corr_func/', a_file='*gsl*emu*.dat')
    # files_lin = try_get_zs_files(
    #     stack_info.last, 'corr_func/', a_file='*gsl*lin*.dat')[1]
    # plot.plot_corr_func_from_data(
    #     xi_list, zs, stack_info, files_lin,
    #     files_emu, zs_emu)


def stack_all(in_dir='/home/vrastil/Documents/GIT/Adhesion-Approximation/output/'):
    # get all runs
    files = get_files_in_traverse_dir(in_dir, 'sim_param.json')
    sim_infos = []
    for a_file, sub in files:
        sim_infos.append(SimInfo(a_file))

    # separate files according to run parameters
    sep_files = []
    for a_sim_info in sim_infos:
        insert(a_sim_info, sep_files)

    # print info about separated files
    num_all_runs = 0
    num_all_sep_runs = 0
    num_sep_runs = 0
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

    for sep_sim_infos in sep_files:
        print "\nStacking group %s" % sep_sim_infos[0].info_tr()
        try:
            stack_group(sep_sim_infos)
        except:
            print("Unexpected error:", sys.exc_info()[0])
            print("Continuing with next group")
