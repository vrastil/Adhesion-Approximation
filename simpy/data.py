"""
'data.py' module serves for loading and saving all necessary data
run analysis of runs
"""

import os
import sys
import traceback
import fnmatch
import numpy as np
from scipy.optimize import curve_fit

from . import plot
from .fastsim import Extrap_Pk_Nl
from .struct import SimInfo
from .power import hybrid_pow_spec, get_Data_vec, corr_func

# *******************
# FIND, SORT, SLICE *
# *******************

def del_duplicate(seq):
    # type: (Seq) -> List
    """remove duplicate elements in sequence while preserving order, return list"""
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]

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

def sort_get_z(files, a_sim_info):
    # type: (List[str], SimInfo) -> List[str], List[TypeVar(str, float)]
    zs = []
    for a_file in files:
        if a_sim_info.app + '_z' in a_file:
            zs.append(float(a_file[a_file.index(a_sim_info.app + '_z') + len(a_sim_info.app+'_z'):-4]))
        elif a_sim_info.app + '_init' in a_file:
            zs.append('init')
        else:
            print "WARNING! Skipping file '%s', unknown format." % a_file
    return zip(*sorted(zip(zs, files), reverse=True))

def sort_get_fl_get_z(a_sim_info, subdir, patterns='*.dat'):
    return sort_get_z(get_files_in_traverse_dir(a_sim_info.dir + subdir, patterns), a_sim_info)

def try_get_zs_files(a_sim_info, subdir, patterns):
    try:
        zs, files = sort_get_fl_get_z(a_sim_info, subdir, patterns=patterns)
        return list(zs), list(files)
    except ValueError:
        print "WARNING! Missing data in '%s'. Skipping step." % (a_sim_info.dir + subdir)
        return None, None

# ***************************
# LOAD DATA FROM SINGLE RUN *
# ***************************

def load_k_supp_from_data(supp_list, k_nyquist_par):
    """
    Divide available k-values into 7 subinterval from
    k_min = 2*PI / L to k_max = 50% k_nyquist_par
    large scale :: k = 1st subinterval
    medium scale :: k = 4rd subinterval
    small scale :: k = 7th subinterval """

    supp = [[[] for x in xrange(3)] for y in xrange(3)]
    for i, j in enumerate([0, 3, 6]):
        for data in supp_list:
            k = data[0]
            P_diff = data[1]
            idx = (np.abs(k-0.5*k_nyquist_par)).argmin() / 7
            supp[i][0].append(np.mean(P_diff[j*idx:(j+1)*idx]))
            supp[i][1].append(np.std(P_diff[j*idx:(j+1)*idx]))
        supp[i][2] = [k[j*idx], k[(j+1)*idx]]
    return supp

def load_k_supp(files, k_nyquist_par):
    supp_list = [np.transpose(np.loadtxt(a_file)) for a_file in files]
    return load_k_supp_from_data(supp_list, k_nyquist_par)

def get_extrap_pk(a_sim_info, files):
    if "extrap_pk" in a_sim_info.data:
        return
    a_sim_info.data["extrap_pk"] = [
        get_hybrid_pow_spec_amp(
            a_sim_info.sim, np.transpose(np.loadtxt(x)),
            a_sim_info.k_nyquist["particle"])
        for x in files]
    a_sim_info.data["pk_list"] = [x["Pk_NL"] for x in a_sim_info.data["extrap_pk"]]

def load_plot_pwr(files, zs, a_sim_info, **kwargs):
    data_list = [np.transpose(np.loadtxt(x)) for x in files]
    get_extrap_pk(a_sim_info, files)
    plot.plot_pwr_spec(data_list, zs, a_sim_info, a_sim_info.data["pk_list"], **kwargs)

def get_plot_corr(files, zs, a_sim_info, load=False):
    if "corr_func" not in a_sim_info.data:
        a_sim_info.data["corr_func"] = {}
        if load:
            a_sim_info.data["corr_func"]["par"] = [np.transpose(np.loadtxt(a_file)) for a_file in files]
        else:
            get_extrap_pk(a_sim_info, files)
            a_sim_info.data["corr_func"]["par"] = [corr_func(a_sim_info.sim, Pk=a_sim_info.data["pk_list"]) for z in zs]
        a_sim_info.data["corr_func"]["lin"] = [corr_func(a_sim_info.sim, z=z) for z in zs]
        a_sim_info.data["corr_func"]["nl"] = [corr_func(a_sim_info.sim, z=z, non_lin=True) for z in zs]

    plot.plot_corr_func(a_sim_info.data["corr_func"], zs, a_sim_info)

def get_plot_supp(files, zs, a_sim_info, pk_type='dens'):
    supp = load_k_supp(files, a_sim_info.k_nyquist["particle"])
    a = [1./(z+1.) for z in zs]
    plot.plot_supp_lms(supp, a, a_sim_info, pk_type=pk_type)

def load_plot_pwr_spec_diff(files, zs, a_sim_info, **kwargs):
    data_list = [np.transpose(np.loadtxt(a_file))
                 for a_file in files]
    plot.plot_pwr_spec_diff_from_data(data_list, zs, a_sim_info, **kwargs)

def split_particle_files(files, zs):
    files_t = [x for x in files if 'track' in x.split('/')[-1]]
    files = [x for x in files if 'par_cut' in x.split('/')[-1]]
    zs = del_duplicate(zs)
    if len(files_t) != len(files):
        raise IndexError("Different number of 'particle' and 'track' files.")
    else:
        return files, files_t, zs

def split_plot_par_slice(files, zs, a_sim_info):
    files, files_t, zs = split_particle_files(files, zs)
    plot.plot_par_last_slice(files, files_t, zs, a_sim_info)

def split_plot_par_evol(files, zs, a_sim_info):
    files, files_t, zs = split_particle_files(files, zs)
    plot.plot_par_evol(files, files_t, zs, a_sim_info)

def load_plot_dens_histo(files, zs, a_sim_info):
    data_list = [np.transpose(np.loadtxt(a_file))
                 for a_file in files]
    plot.plot_dens_histo(data_list, zs, a_sim_info)

def load_check_plot(a_sim_info, key, patterns, # type: SimInfo, str, str,
                    rerun, skip, info_str,     # type: List[str], List[str], str
                    plot_func, subdir=None,    # type: Callable[List[str], List[str], kwargs], str (optional)
                    **kwargs                   # type: Dict (optional)
                   ):
    """bring structure of all steps together:
    1) load redshifts 'zs' and data files 'files'
        - subdirectory is key + '/' if not specified in arguments
    2) check for succes in loading files and if the step should be performed (rerun / skip)
    3) print step name and 'info_str', not printing when None
    4) plot -- need to pass Callable function with arguments: files, zs, a_sim_info, kwargs
    5) write info about done step into a_sim_info
    """
    subdir = key + '/' if subdir is None else subdir
    zs, files = try_get_zs_files(a_sim_info, subdir, patterns)
    if a_sim_info.rerun(rerun, key, skip, zs):
        if info_str is not None:
            print 'step: %s %s' % (key, info_str)
        plot_func(files, zs, a_sim_info, **kwargs)
        a_sim_info.done(key)

# ****************************
# RUN ANALYSIS -- SINGLE RUN *
# ****************************

def analyze_run(a_sim_info, rerun=None, skip=None):
    if skip is None:
        skip = []
    elif skip == "ani":
        skip = ["par_ani", "dens_ani"]
    if rerun is None:
        rerun = []
    
    # Steps to perform -- each entry represents full information needed to perform one step
    # type: Tuple[step_key, data_file_patterns, plot_func, dict_kwargs]
    all_steps = [
        # Power spectrum
        ("pwr_spec", '*par*.dat *init*.dat', '', load_plot_pwr, {}),
        # Power spectrum difference -- input, hybrid, particle
        ("pwr_diff", '*par*', '(particle)', load_plot_pwr_spec_diff, {'ext_title' : 'par'}),
        ("pwr_diff", '*hybrid*', '(hybrid)', load_plot_pwr_spec_diff, {'ext_title' : 'hybrid'}),
        ("pwr_diff", '*input*', '(input)', load_plot_pwr_spec_diff, {'ext_title' : 'input'}),
        # Power spectrum suppression
        ("pwr_spec_supp", '*par*', '', get_plot_supp, {'subdir' : 'pwr_diff/'}),
        # Velocity power spectrum
        ("vel_pwr_spec", '*.dat', '', load_plot_pwr, {'pk_type' : 'vel'}),
        # Velocity power spectrum difference
        ("vel_pwr_diff", '*.dat', '', load_plot_pwr_spec_diff, {'pk_type' : 'vel'}),
        # Velocity power spectrum suppression
        ("vel_pwr_spec_supp", '*.dat', '', get_plot_supp, {'subdir' : 'vel_pwr_diff/', 'pk_type' : 'vel'}),
        # Correlation function
        ("corr_func", '*par*.dat *init*.dat', '', get_plot_corr, {'subdir' : 'pwr_spec/'}),
        # Density distribution
        ("dens_hist", '*.dat', '', load_plot_dens_histo, {'subdir' : 'rho_bin/'}),
        # Particles -- last slice, evolution
        ("par_slice", 'par*.dat track*.dat', '', split_plot_par_slice, {'subdir' : 'par_cut/'}),
        ("par_ani", 'par*.dat track*.dat', '', split_plot_par_evol, {'subdir' : 'par_cut/'}),
        # Density -- two slices, evolution
        ("dens_slice", '*.dat', '', plot.plot_dens_two_slices, {'subdir' : 'rho_map/'}),
        ("dens_ani", '*.dat', '', plot.plot_dens_evol, {'subdir' : 'rho_map/'})
    ]

    # perform all steps
    for key, patterns, info_str, plot_func, kwargs  in all_steps:
        try:
            load_check_plot(a_sim_info, key, patterns, rerun,
                            skip, info_str, plot_func, **kwargs)
        except KeyboardInterrupt:
            raise
        except Exception:
            print "=" * 110
            traceback.print_exc(file=sys.stdout)
            print "=" * 110

def analyze_all(out_dir='/home/vrastil/Documents/GIT/Adhesion-Approximation/output/',
                rerun=None, skip=None, only=None):
    files = get_files_in_traverse_dir(out_dir, 'sim_param.json')
    sim_infos = []
    for a_file in files:
        sim_infos.append(SimInfo(a_file))

    if only is not None:
        sim_infos = sim_infos[only]

    for a_sim_info in sim_infos:
        print 'Analyzing run %s' % a_sim_info.info_tr()
        try:
            analyze_run(a_sim_info, rerun=rerun, skip=skip)
        except KeyboardInterrupt:
            return
        except Exception:
            print "=" * 110
            traceback.print_exc(file=sys.stdout)
            print "=" * 110
    print 'All runs analyzed!'

# ******************************
# LOAD DATA FROM MULTIPLE RUNS *
# ******************************

def load_data_for_stack(stack_info, subdir, a_file):
    # load everything
    all_data_k = None
    all_data_Pk = None
    for a_sim_info in stack_info:
        # load files for ONE run
        zs, files = try_get_zs_files(a_sim_info, subdir, patterns=a_file)

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
    return zs, all_data_k, all_data_Pk

def check_data_consistency(all_data_k, all_data_Pk):
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

def check_data_consistency_diff(data_list):
    """check if all data have the same lengths"""
    try:
        data_array = np.array(data_list)
    except ValueError:
        print '\t\tData in data_list have different shapes. Trying to cut...'
        # convert to list so elements can be deleted
        data_list = [np.array(data).tolist() for data in data_list]
        del_num = 0
        j = 0
        while True:
            k_row = [data[0][j] for data in  data_list]
            k_max = np.max(k_row)
            for ik, k in  enumerate(k_row):
                k_ = k
                while not np.isclose(k_, k_max, rtol=1.e-5, atol=1.e-5) and k_ < k_max :
                    # remove k, Pk if not the same, look for first close
                    if j == len(data_list[ik][0]):
                        break
                    del_num += 1
                    del data_list[ik][0][j]
                    del data_list[ik][1][j]
                    try:
                        del data_list[ik][2][j]
                    except IndexError:
                        pass
                    # look at the next k
                    k_ = data_list[ik][0][j]
            
            j += 1
            # check if j is length of ALL arrays
            for x in data_list:
                if len(x[0]) != j:
                    break
            else:
                break
        if del_num:
            print "\t\tDeleted %i excess values." % (del_num)
        data_array = np.array(data_list)
    return data_array

def stack_files(stack_info, subdir, a_file):
    """ assume data in files are x = data[:, 0], y = data[:, 1]
    return tuple (zs, data_list) where each entry in 'data_list' corresponds to entry in 'zs' list
    data[:, 0] = x, data[:, 1] = mean(y), data[:, 2] = std(y)
    """
    # load everything
    zs, all_data_k, all_data_Pk = load_data_for_stack(stack_info, subdir, a_file)

    # chceck lengths of lists, delete excess (in case outputs of simulations encountered any errors)
    check_data_consistency(all_data_k, all_data_Pk)

    # compute means, std
    data_list = []
    for i in xrange(len(all_data_k)):
        data_list.append([])

        data_list[i].append(np.mean(all_data_k[i], axis=0))
        data_list[i].append(np.mean(all_data_Pk[i], axis=0))
        data_list[i].append(np.std(all_data_Pk[i], axis=0))

    return zs, data_list

def get_perr_pcor(pcov):
    """convert covariance matrix into standard errors and correlation matrix"""
    perr = np.sqrt(np.diag(pcov))
    inv = np.diag([1/x if x else 0 for x in perr]) # cut out values with 0 variance
    pcor = np.dot(np.dot(inv, pcov), inv)
    return perr, pcor

def get_hybrid_pow_spec_amp(sim, data, k_nyquist_par):
    """ fit data [k, Pk, std] to hybrid power spectrum (1-A)*P_lin(k) + A*P_nl(k)
    return dictionary with C++ class Extrap_Pk_Nl, fit values and covariance """
    # extract data
    kk, Pk = data[0], data[1]
    
    # get proper slice of data -- last decade befor half of particle nyquist
    idx = (np.abs(kk-0.5*k_nyquist_par)).argmin()
    idx = slice(idx - sim.out_opt.bins_per_decade, idx)

    # define functions which will be used in fitting
    pk_hyb_func = lambda k, a, A: hybrid_pow_spec(a, k, A, sim.cosmo)
    
    # fit data, a = <0, 1>, A = <0, 1>
    sigma = data[2][idx] if len(data) > 2 else None
    bounds = ([0, 0], [1, 1])
    popt, pcov = curve_fit(pk_hyb_func, kk[idx], Pk[idx], bounds=bounds, sigma=sigma)
    perr, pcor = get_perr_pcor(pcov)

    # get hybrid Extrap
    Pk_par = Extrap_Pk_Nl(get_Data_vec(data), sim, popt[1], popt[0])

    # return all info in dict
    return {"Pk_par" : Pk_par, "popt" : popt, "pcov" : pcov, 'perr' : perr, 'pcor' : pcor}

# ***************************
# SAVE ALREADY STACKED DATA *
# ***************************

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

def save_corr(zs, xi_list, out_dir, app, corr_type=''):
    for i in xrange(len(zs)):
        z_str = 'init.dat' if zs[i] == 'init' else 'z%.2f.dat' % zs[i]
        fname = out_dir + \
            "corr_func_scipy_qawf_%s_%s_%s" % (corr_type, app, z_str)
        header = ("This file contains correlation function depending on distance r in units [Mpc/h]."
                  "r [Mpc/h]\txsi(r)")
        np.savetxt(fname, np.transpose(xi_list[i]), fmt='%.6e', header=header)

# ***************************
# LOAD ALREADY STACKED DATA *
# ***************************

def load_pwr(stack_info, subdir='pwr_spec/', a_file='pwr_spec_par*.dat'):
    zs, files = try_get_zs_files(
        stack_info, subdir, patterns=a_file)
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
            zs_diff, data_list_diff[diff_type] = stack_files(
                stack_info, 'pwr_diff/', '*%s*.dat' % diff_type)
            save_pwr_diff(zs_diff, data_list_diff, diff_type,
                          stack_info.pwr_diff_dir, stack_info.app)
        else:
            zs_diff, data_list_diff[diff_type] = load_pwr(
                stack_info, subdir='pwr_diff/', a_file='*%s*.dat' % diff_type)
    stack_info.done(key)
    return zs_diff, data_list_diff

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

# **********************************
# RUN ANALYSIS -- STACKING OF RUNS *
# **********************************

def stack_group(group_sim_infos=None, stack_info_file=None, plot_all=True, **kwargs):
    # load & save all info about stack
    stack_info = StackInfo(group_sim_infos=group_sim_infos, stack_info_file=stack_info_file, **kwargs)

    print "\tLoad & save data to stack..."
    # load data from stacked files if available, otherwise do stacking over all runs
    zs, data_list = load_pwr_par(stack_info)
    zs_diff, data_list_diff = load_pwr_diff(stack_info)

    # for each entry in data_list compute its continuous version
    # get non-linear power spectra z == 'init' transformed into z = 0.0
    Pk_hyb_lis = [get_hybrid_pow_spec_amp(stack_info.sim, data, stack_info.k_nyquist["particle"]) for data in data_list]
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
        supp = load_k_supp_from_data(data_list_diff["par"], stack_info.k_nyquist["particle"])
        plot.plot_supp_lms(supp, a, stack_info)

        print '\tPlotting power spectrum suppression (map)...'
        # WARNING! data in data_list_diff are potentially modified -- aligning values to the same counts
        data_array = check_data_consistency_diff(data_list_diff["par"])
        plot.plot_pwr_spec_diff_map_from_data(data_array, zs_diff, stack_info, ext_title="par")

        print '\tPlotting correlation function...'
        plot.plot_corr_func_from_data(xi_all, zs, stack_info)

    return stack_info

def stack_all(in_dir='/home/vrastil/Documents/GIT/Adhesion-Approximation/output/', plot_all=True, **kwargs):
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
            stack_infos.append(stack_group(sep_sim_infos, plot_all=plot_all, **kwargs))
        except:
            print "=" * 110
            traceback.print_exc(file=sys.stdout)
            print "=" * 110
    return stack_infos
