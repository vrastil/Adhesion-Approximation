from datetime import datetime
import sys
import gc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.optimize import curve_fit
import json

from . import *
from . import plot
from .data import SimInfo


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
    all_data_k = None
    all_data_Pk = None
    for a_sim_info in group_sim_infos:
        zs, files = try_get_zs_files(a_sim_info, subdir, a_file=a_file)
        if all_data_k is None:
            all_data_k = [[] for x in xrange(len(zs))]
        if all_data_Pk is None:
            all_data_Pk = [[] for x in xrange(len(zs))]
        for i, data_file in enumerate(files):
            data = np.loadtxt(data_file)
            k, P_k = data[:, 0], data[:, 1]
            all_data_k[i].append(k)
            all_data_Pk[i].append(P_k)

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

        f = lambda k_, A_: A_ * plot.pwr_spec(k_, self.cosmo, self.ccl_cosmo)
        popt, pcov = curve_fit(
            f, k[0:decade], Pk[0:decade], sigma=dPk[0:decade])  # fit over a decade
        self.A_lower = popt
        self.k_lower = k[decade / 2]
        ind = (np.abs(k - 0.5 * k_max)).argmin() if k_max is not None else -1
        popt, pcov = curve_fit(f, k[ind - decade / 2:ind], Pk[ind - decade / 2:ind],
                               sigma=dPk[ind - decade / 2:ind])  # fit over a half of a decade
        self.A_upper = popt
        self.k_upper = k[ind]
        print "k < %e, A = %e; k > %e, A = %e" % (self.k_lower, self.A_lower, self.k_upper, self.A_upper)
        del k, Pk, dPk, f

    def eval(self, k):
        if k < self.k_lower:
            return self.A_lower*plot.pwr_spec(k, self.cosmo, self.ccl_cosmo)
        elif k < self.k_upper:
            return self.interpolate(k)
        else:
            return self.A_upper*plot.pwr_spec(k, self.cosmo, self.ccl_cosmo)

    def eval_list(self, k_list):
        Pk_list = []
        for k in k_list:
            Pk_list.append(self.eval(k))
        return Pk_list


def stack_group(group_sim_infos):
    # Power spectrum
    a_sim_info = group_sim_infos[0]
    out_dir = a_sim_info.dir.replace(a_sim_info.dir.split("/")[-2] + "/", "")
    out_dir += "STACK_%im_%ip_%iM_%ib/" % (
        a_sim_info.box_opt["mesh_num"], a_sim_info.box_opt["par_num"],
        a_sim_info.box_opt["mesh_num_pwr"], a_sim_info.box_opt["box_size"])
    create_dir(out_dir)
    create_dir(out_dir + "pwr_spec/")
    create_dir(out_dir + "results/")

    print "\tLoading data to stack..."
    zs, data_list = stack_files(group_sim_infos, 'pwr_spec/', '*par*.dat')
    print "\tSaving stacked data..."
    for i in xrange(len(zs)):
        fname = out_dir + \
            "pwr_spec/pwr_spec_par_%s_z%.2f" % (a_sim_info.app, zs[i])
        header = ("# This file contains power spectrum P(k) in units [(Mpc/h)^3] "
                  "depending on wavenumber k in units [h/Mpc] with standard deviation in units [h/Mpc].\n"
                  "# k [h/Mpc]\tP(k) [(Mpc/h)^3]\tstd [(Mpc/h)^3]")
        np.savetxt(fname, np.transpose(
            data_list[i]), fmt='%.6e', header=header)
        
        Pk = Extrap_Pk(data_list[i], a_sim_info.cosmo,
                       a_sim_info.ccl_cosmo, a_sim_info.k_nyquist["particle"])
        fname = out_dir + \
            "pwr_spec/pwr_spec_extrap_%s_z%.2f" % (a_sim_info.app, zs[i])
        header = ("# This file contains extra/intra-polated power spectrum P(k) in units [(Mpc/h)^3] "
                  "depending on wavenumber k in units [h/Mpc].\n"
                  "# k [h/Mpc]\tP(k) [(Mpc/h)^3]")
        np.savetxt(fname, np.transpose([
            data_list[i][0], Pk.eval_list(data_list[i][0])]), fmt='%.6e', header=header)
    
        for k in np.logspace(-3, 1, num=10):
            print "A = [%e, %e]\tk = %e\tP(k) = %e\tP_lin(k) = %e" % (Pk.A_lower, Pk.A_upper, k, Pk.eval(k),
                plot.pwr_spec(k, a_sim_info.cosmo, a_sim_info.ccl_cosmo, z=zs[i]))
        del Pk

    print '\tPlotting power spectrum...'
    plot.plot_pwr_spec_stacked(
        data_list, zs, a_sim_info, out_dir=out_dir + "results/")


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
        stack_group(sep_sim_infos)
