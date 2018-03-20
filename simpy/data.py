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
from .fastsim import Extrap_Pk_Nl_2, Extrap_Pk_Nl_3
from .struct import SimInfo, StackInfo, insert
from .power import hybrid_pow_spec, get_Data_vec, corr_func, chi_trans_to_supp

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
        print "\tWARNING! Missing data. Skipping step."
        return None, None

def sort_chi_files(files, zs):
    """ separate chi_files (*_chi_*) from given files, return 4 lists """
    files_chi, zs_chi = zip(*[x for x in zip(files, zs) if 'chi' in x[0]])
    files, zs = zip(*[x for x in zip(files, zs) if 'chi' not in x[0]])
    return map(list, [files, zs, files_chi, zs_chi])

# ***************************
# LOAD DATA FROM SINGLE RUN *
# ***************************

def load_k_supp(files, k_nyquist_par, a_sim_info=None, a=None, pk_type='dens'):
    """
    Divide available k-values into 7 subinterval from
    k_min = 2*PI / L to k_max = 50% k_nyquist_par
    large scale :: k = 1st subinterval
    medium scale :: k = 4rd subinterval
    small scale :: k = 7th subinterval 
    
    For chameleon field divide P(k) values by linear prediction """

    supp = [[[] for x in xrange(3)] for y in xrange(3)]
    for i, j in enumerate([0, 3, 6]):
        for ii, a_file in enumerate(files):
            data = np.transpose(np.loadtxt(a_file))
            k = data[0]
            P_diff = data[1]
            if pk_type == 'chi':
                P_diff = chi_trans_to_supp(a[ii], k, P_diff, a_sim_info.sim.cosmo, a_sim_info.chi_opt)
            idx = (np.abs(k-0.5*k_nyquist_par)).argmin() / 7
            supp[i][0].append(np.mean(P_diff[j*idx:(j+1)*idx]))
            supp[i][1].append(np.std(P_diff[j*idx:(j+1)*idx]))
        supp[i][2] = [k[j*idx], k[(j+1)*idx]]
    return supp

def get_extrap_pk(a_sim_info, files):
    if "extrap_pk" in a_sim_info.data:
        return
    a_sim_info.data["extrap_pk"] = [
        get_hybrid_pow_spec_amp(
            a_sim_info.sim, np.transpose(np.loadtxt(x)),
            a_sim_info.k_nyquist["particle"])
        for x in files]
    a_sim_info.data["pk_list"] = [x["Pk_par"] for x in a_sim_info.data["extrap_pk"]]

def load_plot_pwr(files, zs, a_sim_info, **kwargs):
    data_list = [np.transpose(np.loadtxt(x)) for x in files]
    get_extrap_pk(a_sim_info, files)
    plot.plot_pwr_spec(data_list, zs, a_sim_info, a_sim_info.data["pk_list"], **kwargs)

def load_plot_chi_pwr(files, zs, a_sim_info, **kwargs):
    data_list = [np.transpose(np.loadtxt(x)) for x in files]
    plot.plot_chi_pwr_spec(data_list, zs, a_sim_info, **kwargs)

def get_plot_corr(files, zs, a_sim_info, load=False):
    if "corr_func" not in a_sim_info.data:
        a_sim_info.data["corr_func"] = {}
        if load:
            a_sim_info.data["corr_func"]["par"] = [np.transpose(np.loadtxt(a_file)) for a_file in files]
        else:
            get_extrap_pk(a_sim_info, files)
            a_sim_info.data["corr_func"]["par"] = [corr_func(a_sim_info.sim, Pk=Pk) for Pk in a_sim_info.data["pk_list"]]
        a_sim_info.data["corr_func"]["lin"] = [corr_func(a_sim_info.sim, z=z) for z in zs]
        a_sim_info.data["corr_func"]["nl"] = [corr_func(a_sim_info.sim, z=z, non_lin=True) for z in zs]

    plot.plot_corr_func(a_sim_info.data["corr_func"], zs, a_sim_info)

def get_plot_supp(files, zs, a_sim_info, pk_type='dens'):
    a = [1./(z+1.) for z in zs]
    supp = load_k_supp(files, a_sim_info.k_nyquist["particle"], a_sim_info=a_sim_info, a=a, pk_type=pk_type)
    plot.plot_supp_lms(supp, a, a_sim_info, pk_type=pk_type)

def get_plot_supp_map(files, zs, a_sim_info, pk_type='dens'):
    data_array = [np.transpose(np.loadtxt(a_file)) for a_file in files]
    data_array = check_data_consistency_diff(data_array)
    plot.plot_pwr_spec_diff_map_from_data(data_array, zs, a_sim_info, ext_title="par")

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
                    rerun, skip, plot_func,    # type: List[str], List[str], Callable[List[str], List[str], kwargs],
                    info_str='', subdir=None,  # type: str, str
                    **kwargs                   # type: Dict
                   ):
    """bring structure of all steps together:
    1) load redshifts 'zs' and data files 'files'
        - subdirectory is key + '/' if not specified in arguments
    2) check for success in loading files and if the step should be performed (rerun / skip)
    3) print step name and 'info_str', not printing when None
    4) plot -- need to pass Callable function with arguments: files, zs, a_sim_info, kwargs
    5) write info about done step into a_sim_info
    """
    print 'step: %s %s' % (key, info_str)
    subdir = key + '/' if subdir is None else subdir
    zs, files = try_get_zs_files(a_sim_info, subdir, patterns)
    if a_sim_info.rerun(rerun, key, skip, zs):
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
    # type: Tuple[step_key, data_file_patterns, plot_func, opt_kwargs]
    all_steps = [
        # Power spectrum -- particle, velocity, chameleon
        ("pwr_spec", '*par*.dat *init*.dat', load_plot_pwr, {}),
        ("pwr_spec_chi", '*chi*.dat*', load_plot_chi_pwr, {'subdir' : 'pwr_spec/'}),
        ("vel_pwr_spec", '*.dat', load_plot_pwr, {'pk_type' : 'vel'}),
        # Power spectrum difference -- input, hybrid, particle, velocity, chameleon
        ("pwr_diff", '*par*', load_plot_pwr_spec_diff,
            {'info_str' : '(particle)', 'ext_title' : 'par'}),
        ("pwr_diff_h", '*hybrid*', load_plot_pwr_spec_diff,
            {'subdir' : 'pwr_diff/', 'info_str' : '(hybrid)', 'ext_title' : 'hybrid'}),
        ("pwr_diff_i", '*input*', load_plot_pwr_spec_diff,
            {'subdir' : 'pwr_diff/', 'info_str' : '(input)', 'ext_title' : 'input'}),
        ("vel_pwr_diff", '*.dat', load_plot_pwr_spec_diff, {'pk_type' : 'vel'}),
        ("chi_pwr_diff", '*chi*.dat*', load_plot_pwr_spec_diff,
            {'subdir' : 'pwr_spec/', 'pk_type' : 'chi'}),
        # Power spectrum suppression (includes maps) -- particle, velocity, chameleon,
        ("pwr_spec_supp", '*par*', get_plot_supp, {'subdir' : 'pwr_diff/'}),
        ("pwr_spec_supp_map", '*par*', get_plot_supp_map, {'subdir' : 'pwr_diff/'}),
        ("vel_pwr_spec_supp", '*.dat', get_plot_supp, {'subdir' : 'vel_pwr_diff/', 'pk_type' : 'vel'}),
        ("chi_pwr_spec_supp", '*chi*.dat*', get_plot_supp, {'subdir' : 'pwr_spec/', 'pk_type' : 'chi'}),
        # Correlation function
        ("corr_func", '*par*.dat *init*.dat', get_plot_corr, {'subdir' : 'pwr_spec/'}),
        # Density distribution
        ("dens_hist", '*.dat', load_plot_dens_histo, {'subdir' : 'rho_bin/'}),
        # Particles -- last slice, evolution
        ("par_slice", 'par*.dat track*.dat', split_plot_par_slice, {'subdir' : 'par_cut/'}),
        ("par_ani", 'par*.dat track*.dat', split_plot_par_evol, {'subdir' : 'par_cut/'}),
        # Density -- two slices, evolution
        ("dens_slice", '*.dat', plot.plot_dens_two_slices, {'subdir' : 'rho_map/'}),
        ("dens_ani", '*.dat', plot.plot_dens_evol, {'subdir' : 'rho_map/'})
    ]

    # perform all steps, skip step if Exception occurs
    for key, patterns, plot_func, kwargs in all_steps:
        try:
            load_check_plot(a_sim_info, key, patterns, rerun,
                            skip, plot_func, **kwargs)
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
            print 'Exiting...'
            return
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
    """ assume data in files are k = data[:, 0], Pk = data[:, 1]
    return tuple (zs, data_list) where each entry in 'data_list' corresponds to entry in 'zs' list
    data[:, 0] = k, data[:, 1] = mean(Pk), data[:, 2] = std(Pk)
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
    Pk_par = Extrap_Pk_Nl_3 if sigma else Extrap_Pk_Nl_2
    Pk_par = Pk_par(get_Data_vec(data), sim, popt[1], popt[0])

    # return all info in dict
    return {"Pk_par" : Pk_par, "popt" : popt, "pcov" : pcov, 'perr' : perr, 'pcor' : pcor}

# **************************
# LOAD & SAVE STACKED DATA *
# **************************

HEADER_PWR = ("This file contains power spectrum P(k) in units [(Mpc/h)^3] "
              "depending on wavenumber k in units [h/Mpc] with standard deviation in units [h/Mpc].\n"
              "k [h/Mpc]\tP(k) [(Mpc/h)^3]\tstd [(Mpc/h)^3]")
HEADER_PWR_DIFF_1 = ("This file contains relative difference between power spectrum P(k)\n"
                     "and lineary extrapolated ")
HEADER_PWR_DIFF_2 = ("depending on wavenumber k in units [h/Mpc] with standard deviation in units [h/Mpc].\n"
                     "k [h/Mpc]\tP(k) [(Mpc/h)^3]\tstd [(Mpc/h)^3]")
HEADER_PWR_DIFF_PAR = HEADER_PWR_DIFF_1 + "power spectrum of initial particle position\n" + HEADER_PWR_DIFF_2
HEADER_PWR_DIFF_INPUT = HEADER_PWR_DIFF_1 + "'hybrid' power spectrum\n" + HEADER_PWR_DIFF_2
HEADER_PWR_DIFF_HYBRID = HEADER_PWR_DIFF_1 + "input power spectrum\n" + HEADER_PWR_DIFF_2
HEADER_CORR = ("This file contains correlation function depending on distance r in units [Mpc/h]."
               "r [Mpc/h]\txsi(r)")

def load_stack_save(stack_info, key, patterns,  # type: StackInfo, str, str,
                    rerun, skip, header, fname, # type: List[str], List[str], str, str
                    info_str='', subdir=None    # type: str, str
                   ):
    """bring structure of stacking together:
    1) check if the step should be performed (rerun / skip)
    2) stack files through function 'stack_files'
       - subdirectory is key without '_files' + '/' if not specified in arguments
    3) save stacked files
       - out_dir = stack_info.dir + subdir
       - fname = fname + "_%s_%s" % (app, z_str)
    4) write info about done step into stack_info
    """
    if stack_info.rerun(rerun, key, skip, True):
        print 'step: %s %s' % (key, info_str) 
        subdir = key.replace('_files', '/') if subdir is None else subdir
        for z, data in zip(*stack_files(stack_info, subdir, patterns)):
            z_str = 'init.dat' if z == 'init' else 'z%.2f.dat' % z
            fname_ = stack_info.dir + subdir + fname + "_%s_%s" % (stack_info.app, z_str)
            np.savetxt(fname_, np.transpose(data), fmt='%.6e', header=header)
        stack_info.done(key)

# **********************************
# RUN ANALYSIS -- STACKING OF RUNS *
# **********************************

def stack_group(rerun=None, skip=None, **kwargs):
    if skip is None:
        skip = []
    elif skip == "ani":
        skip = ["par_ani", "dens_ani"]
    if rerun is None:
        rerun = []

    # load & save all info about stack
    stack_info = StackInfo(**kwargs)

    # load, stack & save all files
    all_files_steps = [
        # Power spectrum
        ("pwr_spec_files", '*par*.dat *init*.dat', HEADER_PWR, 'pwr_spec_par', {}),
        # Power spectrum difference -- input, hybrid, particle
        ("pwr_diff_files", '*par*', HEADER_PWR_DIFF_PAR, 'pwr_spec_diff_par', {'info_str' : '(particle)'}),
        ("pwr_diff_files_h", '*hybrid*', HEADER_PWR_DIFF_INPUT, 'pwr_spec_diff_hybrid',
            {'subdir' : 'pwr_diff/', 'info_str' : '(hybrid)'}),
        ("pwr_diff_files_i", '*input*', HEADER_PWR_DIFF_HYBRID, 'pwr_spec_diff_input',
            {'subdir' : 'pwr_diff/', 'info_str' : '(input)'})
    ]
    for key, patterns, header, fname, kwargs in all_files_steps:
        try:
            load_stack_save(stack_info, key, patterns, rerun,
                            skip, header, fname, **kwargs)
        except KeyboardInterrupt:
            raise
        except Exception:
            print "=" * 110
            traceback.print_exc(file=sys.stdout)
            print "=" * 110

    # load and plot files
    all_steps = [
        # Power spectrum
        ("pwr_spec", '*par*.dat *init*.dat', load_plot_pwr, {'err' : True}),
        # Power spectrum difference -- input, hybrid, particle
        # Power spectrum difference -- input, hybrid, particle
        ("pwr_diff", '*par*', load_plot_pwr_spec_diff,
            {'info_str' : '(particle)', 'ext_title' : 'par'}),
        ("pwr_diff_h", '*hybrid*', load_plot_pwr_spec_diff,
            {'subdir' : 'pwr_diff/', 'info_str' : '(hybrid)', 'ext_title' : 'hybrid'}),
        ("pwr_diff_i", '*input*', load_plot_pwr_spec_diff,
            {'subdir' : 'pwr_diff/', 'info_str' : '(input)', 'ext_title' : 'input'}),
        # Correlation function
        ("corr_func", '*par*.dat *init*.dat', get_plot_corr, {'subdir' : 'pwr_spec/'}),
        # Power spectrum suppression
        ("pwr_spec_supp", '*par*', get_plot_supp, {'subdir' : 'pwr_diff/'}),
        ("pwr_spec_supp_map", '*par*', get_plot_supp_map, {'subdir' : 'pwr_diff/'})
    ]

    # perform all steps, skip step if Exception occurs
    for key, patterns, plot_func, kwargs in all_steps:
        try:
            load_check_plot(stack_info, key, patterns, rerun,
                            skip, plot_func, **kwargs)
        except KeyboardInterrupt:
            raise
        except Exception:
            print "=" * 110
            traceback.print_exc(file=sys.stdout)
            print "=" * 110

def stack_all(in_dir='/home/vrastil/Documents/GIT/Adhesion-Approximation/output/', rerun=None, skip=None, **kwargs):
    # get all runs
    files = get_files_in_traverse_dir(in_dir, 'sim_param.json')
    sim_infos = [SimInfo(a_file) for a_file in files]
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

    for sep_sim_infos in sep_files:
        print "\nStacking group %s" % sep_sim_infos[0].info_tr()
        try:
            stack_group(group_sim_infos=sep_sim_infos, rerun=rerun, skip=skip, **kwargs)
        except KeyboardInterrupt:
            print 'Exiting...'
            return
    print 'All groups analyzed!'
