"""
'data.py' module serves for loading and saving all necessary data
run analysis of runs
"""

import os
import sys
import traceback
from pygments import highlight
from pygments.lexers import get_lexer_by_name
from pygments.formatters import TerminalFormatter
from itertools import izip

import fnmatch
import numpy as np
from scipy.optimize import curve_fit

from . import plot
from .fastsim import Extrap_Pk_Nl_2, Extrap_Pk_Nl_3
from . import struct
from . import power as pwr

def print_exception(file=sys.stdout):
    """ print catched exception with colors """
    tbtext = traceback.format_exc()
    lexer = get_lexer_by_name("pytb", stripall=True)
    formatter = TerminalFormatter()

    print "\n"
    print "=" * 110
    file.write(highlight(tbtext, lexer, formatter))
    print "=" * 110            

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

def sort_lists(*lists):
    return zip(*sorted(zip(*lists), reverse=True))

def sort_get_z(files, a_sim_info, skip_init=False):
    # type: (List[str], struct.SimInfo) -> List[str], List[TypeVar(str, float)]
    zs = []
    for a_file in files:
        if a_sim_info.app + '_z' in a_file:
            zs.append(float(a_file[a_file.index(a_sim_info.app + '_z') + len(a_sim_info.app+'_z'):-4]))
        elif a_sim_info.app + '_init' in a_file and not skip_init:
            zs.append('init')
        else:
            print "WARNING! Skipping file '%s', unknown format." % a_file
    return sort_lists(zs, files)

def sort_get_fl_get_z(a_sim_info, subdir, patterns='*.dat', skip_init=False):
    files = get_files_in_traverse_dir(a_sim_info.dir + subdir, patterns)
    return sort_get_z(files, a_sim_info, skip_init=skip_init)

def try_get_zs_files(a_sim_info, subdir, patterns, skip_init=False):
    try:
        zs, files = sort_get_fl_get_z(a_sim_info, subdir, patterns=patterns, skip_init=skip_init)
        return list(zs), list(files)
    except ValueError:
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
                P_diff = -1 + pwr.chi_trans_to_supp(a[ii], k, P_diff, a_sim_info.sim.cosmo, a_sim_info.chi_opt)
            idx = (np.abs(k-0.5*k_nyquist_par)).argmin() / 7
            supp[i][0].append(np.mean(P_diff[j*idx:(j+1)*idx]))
            supp[i][1].append(np.std(P_diff[j*idx:(j+1)*idx]))
        supp[i][2] = [k[j*idx], k[(j+1)*idx]]
    return supp

def get_single_hybrid_pow_spec_amp_w_z(sim, a_file, z, k_nyquist_par):
    data = np.transpose(np.loadtxt(a_file))
    x = get_hybrid_pow_spec_amp(sim, data, k_nyquist_par)
    x["z"] = z
    return x

def get_extrap_pk(a_sim_info, files, zs):
    if "extrap_pk" not in a_sim_info.data:
        # needed variables
        sim = a_sim_info.sim
        k_nyquist_par = a_sim_info.k_nyquist["particle"]
        func = lambda a_file, z : get_single_hybrid_pow_spec_amp_w_z(sim, a_file, z, k_nyquist_par)

        # get hybrid power spectrum with redshift for each file
        a_sim_info.data["extrap_pk"] = list(map(func, files, zs))

        # extract 'Pk_par' and 'z' for convenience
        a_sim_info.data["pk_list"] = [x["Pk_par"] for x in a_sim_info.data["extrap_pk"]]
        a_sim_info.data["zs"] = [x["z"] for x in a_sim_info.data["extrap_pk"]]

def load_plot_pwr(files, zs, a_sim_info, **kwargs):
    data_list = [np.transpose(np.loadtxt(x)) for x in files]
    get_extrap_pk(a_sim_info, files, zs)
    plot.plot_pwr_spec(data_list, zs, a_sim_info, a_sim_info.data["pk_list"], **kwargs)

def load_plot_chi_pwr(files, zs, a_sim_info, **kwargs):
    data_list = [np.transpose(np.loadtxt(x)) for x in files]
    plot.plot_chi_pwr_spec(data_list, zs, a_sim_info, **kwargs)

def load_plot_slope(files, zs, a_sim_info, **kwargs):
    data_list = [np.transpose(np.loadtxt(x)) for x in files]
    get_extrap_pk(a_sim_info, files, zs)
    plot.plot_slope(data_list, zs, a_sim_info, a_sim_info.data["pk_list"], **kwargs)

def get_key_func(files, zs, a_sim_info, key, load=False):
    """ key must be name of callable function placed in power with signature:
        key(sim, Pk=None, z=None, non_lin=False) """

    if key not in a_sim_info.data:
        a_sim_info.data[key] = {}
        gen_func = getattr(pwr, key)
        if load:
            a_sim_info.data[key]["par"] = [np.transpose(np.loadtxt(a_file)) for a_file in files]
        else:
            get_extrap_pk(a_sim_info, files, zs)
            a_sim_info.data[key]["par"] = [gen_func(a_sim_info.sim, Pk=Pk) for Pk in a_sim_info.data["pk_list"]]
        a_sim_info.data[key]["lin"] = [gen_func(a_sim_info.sim, z=z) for z in zs]
        a_sim_info.data[key]["nl"] = [gen_func(a_sim_info.sim, z=z, non_lin=True) for z in zs]
        a_sim_info.data[key]["zs"] = zs

def get_corr_func(files, zs, a_sim_info, load=False):
    get_key_func(files, zs, a_sim_info, "corr_func", load=load)

def get_plot_corr(files, zs, a_sim_info, load=False, **kwargs):
    get_corr_func(files, zs, a_sim_info, load=load)
    plot.plot_corr_func(a_sim_info.data["corr_func"], zs, a_sim_info, **kwargs)
    
def get_sigma_R(files, zs, a_sim_info, load=False):
    get_key_func(files, zs, a_sim_info, "sigma_R", load=load)

def get_plot_sigma(files, zs, a_sim_info, load=False):
    get_sigma_R(files, zs, a_sim_info, load=load)
    plot.plot_corr_func(a_sim_info.data["sigma_R"], zs, a_sim_info, is_sigma=True)

def find_nearest_idx(array, value, axis=None):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin(axis=axis)
    return idx

def get_plot_supp(files, zs, a_sim_info, pk_type='dens'):
    a = [1./(z+1.) for z in zs]
    supp = load_k_supp(files, a_sim_info.k_nyquist["particle"], a_sim_info=a_sim_info, a=a, pk_type=pk_type)
    plot.plot_supp_lms(supp, a, a_sim_info, pk_type=pk_type)

def get_plot_supp_map(files, zs, a_sim_info, pk_type='dens'):
    data_array = [np.transpose(np.loadtxt(a_file)) for a_file in files]
    data_array = check_data_consistency_diff(data_array)
    plot.plot_pwr_spec_diff_map_from_data(data_array, zs, a_sim_info, ext_title="par", pk_type=pk_type)

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

def load_check_plot(a_sim_info, key, patterns, # type: struct.SimInfo, str, str,
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
    print 'step: %-25s' % (key + ' ' + info_str),
    sys.stdout.flush()
    subdir = key + '/' if subdir is None else subdir
    zs, files = try_get_zs_files(a_sim_info, subdir, patterns)
    if a_sim_info.rerun(rerun, key, skip, zs):
        plot_func(files, zs, a_sim_info, **kwargs)
        a_sim_info.done(key)

def cut_zs_files(zs, files, z=None):
    if z is not None:
        if z == 'init':
            zs = zs[:1]
            files = files[:1]
        else:
            zs = zs[1:]
            files = files[1:]
            idx = find_nearest_idx(zs, z)
            zs = [zs[idx]]
            files = [files[idx]]
    return zs, files

def get_initialized_StackInfo(a_file, z=None, get_corr=False, get_sigma=False, get_data=False):
    # get struct.StackInfo
    a_sim_info = struct.StackInfo(stack_info_file=a_file)

    # get files for extrapolated power spectrum
    zs, files = try_get_zs_files(a_sim_info, subdir='pwr_spec/', patterns='*par*.dat')

    # if z is specified, find nearest-value file
    zs, files = cut_zs_files(zs, files, z)

    # get data
    get_extrap_pk(a_sim_info, files, zs)

    # get correlation function
    if get_corr:
        get_corr_func(files, zs, a_sim_info)

    # get amplitude of density fluctuations
    if get_sigma:
        get_sigma_R(files, zs, a_sim_info)

    # get data about power spectrum
    if get_data:
        a_sim_info.data["pk_data_par"] = [np.transpose(np.loadtxt(x)) for x in files]

    # return initialized object
    return a_sim_info

# ****************************
# RUN ANALYSIS -- SINGLE RUN *
# ****************************

def analyze_run(a_sim_info, rerun=None, skip=None):
    # Steps to perform -- each entry represents full information needed to perform one step
    # type: Tuple[step_key, data_file_patterns, plot_func, opt_kwargs]
    all_steps = [
        # Power spectrum -- particle, velocity, chameleon, slope
        ("pwr_spec", '*par*.dat *init*.dat', load_plot_pwr, {}),
        ("pwr_spec_chi", '*chi*.dat*', load_plot_chi_pwr, {'subdir' : 'pwr_spec/'}),
        ("vel_pwr_spec", '*.dat', load_plot_pwr, {'pk_type' : 'vel'}),
        ("pwr_slope", '*par*.dat*', load_plot_slope, {'subdir' : 'pwr_spec/'}),
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
        ("chi_pwr_spec_supp_map", '*chi*.dat*', get_plot_supp_map, {'subdir' : 'pwr_spec/', 'pk_type' : 'chi'}),
        # Correlation function, amplitude of density fluctuations
        ("corr_func", '*par*.dat *init*.dat', get_plot_corr, {'subdir' : 'pwr_spec/'}),
        ("sigma_R", '*par*.dat *init*.dat', get_plot_sigma, {'subdir' : 'pwr_spec/'}),
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
            print_exception()

def analyze_all(out_dir='/home/vrastil/Documents/GIT/Adhesion-Approximation/output/',
                rerun=None, skip=None, only=None):
    files = get_files_in_traverse_dir(out_dir, 'sim_param.json')
    sim_infos = []
    for a_file in files:
        sim_infos.append(struct.SimInfo(a_file))

    if only is not None:
        sim_infos = sim_infos[only]

    info = ''
    for a_sim_info in sim_infos:
        info = 'Analyzing run %s' % a_sim_info.info_tr()
        print '*'*len(info), '\n', info, '\n', '*'*len(info)
        try:
            analyze_run(a_sim_info, rerun=rerun, skip=skip)
        except KeyboardInterrupt:
            print 'Exiting...'
            return
    print '*'*len(info), '\n','All runs analyzed!'

# ******************************
# LOAD DATA FROM MULTIPLE RUNS *
# ******************************

def load_data_for_stack(stack_info, subdir, a_file):
    # load everything
    all_data_k = None
    all_data_Pk = None
    all_zs = None
    for a_sim_info in stack_info:
        # load files for ONE run
        zs, files = try_get_zs_files(a_sim_info, subdir, patterns=a_file)

        # if data are missing for this run
        if zs is None:
            continue
        else:
            all_zs = zs

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
    return all_zs, all_data_k, all_data_Pk

def check_data_consistency(all_data_k, all_data_Pk):
    # chceck lengths of lists, delete excess (in case outputs of simulations encountered any errors)
    for i, data_k in enumerate(all_data_k):
        j = 0
        while True:
            min_l = np.min([len(k_vec) for k_vec in  data_k])
            if min_l == j:
                # delete all rows beyond this j
                for ik in range(len(data_k)):
                    while len(all_data_k[i][ik]) > min_l:
                        del all_data_k[i][ik][j]
                        del all_data_Pk[i][ik][j]
                break
            k_row = [k_vec[j] for k_vec in  data_k]
            k_max = np.max(k_row)
            for ik, k in  enumerate(k_row):
                k_ = k
                while not np.isclose(k_, k_max, rtol=1.e-5, atol=1.e-5) and k_ < k_max :
                    # remove k, Pk if not the same, look for first close
                    if j == len(all_data_k[i][ik]):
                        break

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

    # if data are missing for all runs
    if zs is None:
        return None, None

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
    kk, Pk = np.array(data[0]), np.array(data[1])
    
    # get proper slice of data -- last decade befor half of particle nyquist
    idx = (np.abs(kk-0.5*k_nyquist_par)).argmin()
    idx = slice(idx - sim.out_opt.bins_per_decade, idx)

    # define functions which will be used in fitting
    pk_hyb_func = lambda k, a, A: pwr.hybrid_pow_spec(a, k, A, sim.cosmo)
    
    # fit data, a = <0, 1>, A = <0, 1>
    sigma = data[2][idx] if len(data) > 2 else None
    bounds = ([0, 0], [1, 1])
    popt, pcov = curve_fit(pk_hyb_func, kk[idx], Pk[idx], bounds=bounds, sigma=sigma)
    perr, pcor = get_perr_pcor(pcov)

    # get hybrid Extrap
    Pk_par = Extrap_Pk_Nl_2 if sigma is None else Extrap_Pk_Nl_3
    Pk_par = Pk_par(pwr.get_Data_vec(data), sim, popt[1], popt[0])

    # return all info in dict
    return {"Pk_par" : Pk_par, "popt" : popt, "pcov" : pcov, 'perr' : perr, 'pcor' : pcor}

def get_a_eff_from_dens_fluct(a_sim_info):
    """ fit data [r, sigma] to """
    # check zs
    zs = a_sim_info.data["sigma_R"]["zs"]
    idx = 1 if zs[0] == 'init' else 0

    # sigma_R_0
    data_par_0 = np.array(a_sim_info.data["sigma_R"]["par"])[idx,1]
    data_nl_0 = np.array(a_sim_info.data["sigma_R"]["nl"])[idx,1]

    # sigma_R
    data_par = np.array(a_sim_info.data["sigma_R"]["par"])[idx:,1]
    data_nl = np.array(a_sim_info.data["sigma_R"]["nl"])[idx:,1]

    # ratio and its std
    data_D_eff = np.sqrt(data_par / data_par_0 * data_nl_0 / data_nl)
    D_eff_ratio = np.mean(data_D_eff, axis=1)
    D_eff_std = np.std(data_D_eff, axis=1)

    # save data
    a_sim_info.data["sigma_R"]["D_eff_ratio"] = D_eff_ratio
    a_sim_info.data["sigma_R"]["D_eff_std"] = D_eff_std

def get_a_eff_from_Pk(stack_info):
    # go through all extrapolated Pk
    a, popt, pcov, perr = map(np.array, zip(*[ # extract back, store as np.array
        (1/(Pk['z'] + 1), Pk['popt'], Pk['pcov'], Pk['perr']) # a, popt, pcov, perr
        for Pk in stack_info.data["extrap_pk"] if Pk["z"] != 'init']))
    
    # derived variables
    a_eff = popt[:,0]
    A = popt[:,1]
    
    # store
    stack_info.data["eff_time"] = {'a' : a, 'popt' : popt, 'pcov' : pcov, 'a_eff' : a_eff, 'A' : A, 'perr' : perr}

# **************************
# LOAD & SAVE STACKED DATA *
# **************************

HEADER_PWR = ("This file contains power spectrum P(k) in units [(Mpc/h)^3] "
              "depending on wavenumber k in units [h/Mpc] with standard deviation in units [h/Mpc].\n"
              "k [h/Mpc]\tP(k) [(Mpc/h)^3]\tstd [(Mpc/h)^3]")
HEADER_PWR_CHI = HEADER_PWR.replace("power", "chameleon power").replace("units", "units of chi_a/(1-n)")
HEADER_PWR_DIFF = ("This file contains relative difference between power spectrum P(k)\n"
                   "and lineary extrapolated <STRING_TO_REPLACE>\n"
                   "depending on wavenumber k in units [h/Mpc] with standard deviation in units [h/Mpc].\n"
                   "k [h/Mpc]\tP(k) [(Mpc/h)^3]\tstd [(Mpc/h)^3]")
HEADER_PWR_DIFF_PAR = HEADER_PWR_DIFF.replace("<STRING_TO_REPLACE>", "power spectrum of initial particle position")
HEADER_PWR_DIFF_INPUT = HEADER_PWR_DIFF.replace("<STRING_TO_REPLACE>", "input power spectrum")
HEADER_PWR_DIFF_HYBRID = HEADER_PWR_DIFF.replace("<STRING_TO_REPLACE>", "'hybrid' power spectrum")
HEADER_CORR = ("This file contains correlation function depending on distance r in units [Mpc/h]."
               "r [Mpc/h]\txsi(r)")

def load_stack_save(stack_info, key, patterns,  # type: struct.StackInfo, str, str,
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
    print 'step: %-25s' % (key + ' ' + info_str),
    sys.stdout.flush()

    # check only rerun / key / skip -- no files loading
    if stack_info.rerun(rerun, key, skip, True):
        subdir = key.replace('_files', '/') if subdir is None else subdir
        zs, data_list = stack_files(stack_info, subdir, patterns)

        # check again if we loaded any data
        if stack_info.rerun(rerun, key, skip, zs):
            for z, data in zip(zs, data_list):
                z_str = 'init.dat' if z == 'init' else 'z%.2f.dat' % z
                fname_ = stack_info.dir + subdir + fname + "_%s_%s" % (stack_info.app, z_str)
                np.savetxt(fname_, np.transpose(data), fmt='%.6e', header=header)
            stack_info.done(key)    

# **********************************
# RUN ANALYSIS -- STACKING OF RUNS *
# **********************************

def stack_group(rerun=None, skip=None, **kwargs):
    # load & save all info about stack
    stack_info = struct.StackInfo(**kwargs)

    # load, stack & save all files
    all_files_steps = [
        # Power spectrum
        ("pwr_spec_files", '*par*.dat *init*.dat', HEADER_PWR, 'pwr_spec_par', {}),
        ("pwr_spec_chi_files", '*chi*.dat*', HEADER_PWR_CHI, 'pwr_spec_chi',
            {'subdir' : 'pwr_spec/'}),
        # Power spectrum difference -- input, hybrid, particle
        ("pwr_diff_files", '*par*', HEADER_PWR_DIFF_PAR, 'pwr_spec_diff_par', {'info_str' : '(particle)'}),
        ("pwr_diff_files_h", '*hybrid*', HEADER_PWR_DIFF_HYBRID, 'pwr_spec_diff_hybrid',
            {'subdir' : 'pwr_diff/', 'info_str' : '(hybrid)'}),
        ("pwr_diff_files_i", '*input*', HEADER_PWR_DIFF_INPUT, 'pwr_spec_diff_input',
            {'subdir' : 'pwr_diff/', 'info_str' : '(input)'})
    ]
    for key, patterns, header, fname, kwargs in all_files_steps:
        try:
            load_stack_save(stack_info, key, patterns, rerun,
                            skip, header, fname, **kwargs)
        except KeyboardInterrupt:
            raise
        except Exception:
            print_exception()

    # load and plot files
    all_steps = [
        # Power spectrum
        ("pwr_spec", '*par*.dat *init*.dat', load_plot_pwr, {'err' : True}),
        ("pwr_spec_chi", '*chi*.dat*', load_plot_chi_pwr,
            {'subdir' : 'pwr_spec/', 'err' : True}),
        ("pwr_slope", '*par*.dat*', load_plot_slope, {'subdir' : 'pwr_spec/'}),
        # Power spectrum difference -- input, hybrid, particle
        ("pwr_diff", '*par*', load_plot_pwr_spec_diff,
            {'info_str' : '(particle)', 'ext_title' : 'par'}),
        ("pwr_diff_h", '*hybrid*', load_plot_pwr_spec_diff,
            {'subdir' : 'pwr_diff/', 'info_str' : '(hybrid)', 'ext_title' : 'hybrid'}),
        ("pwr_diff_i", '*input*', load_plot_pwr_spec_diff,
            {'subdir' : 'pwr_diff/', 'info_str' : '(input)', 'ext_title' : 'input'}),
        ("chi_pwr_diff", '*chi*.dat*', load_plot_pwr_spec_diff,
            {'subdir' : 'pwr_spec/', 'pk_type' : 'chi'}),
        # Correlation function, amplitude of density fluctuations
        ("corr_func", '*par*.dat *init*.dat', get_plot_corr, {'subdir' : 'pwr_spec/'}),
        ("sigma_R", '*par*.dat *init*.dat', get_plot_sigma, {'subdir' : 'pwr_spec/'}),
        # Power spectrum suppression
        ("pwr_spec_supp", '*par*', get_plot_supp, {'subdir' : 'pwr_diff/'}),
        ("pwr_spec_supp_map", '*par*', get_plot_supp_map, {'subdir' : 'pwr_diff/'}),
        ("chi_pwr_spec_supp", '*chi*.dat*', get_plot_supp, {'subdir' : 'pwr_spec/', 'pk_type' : 'chi'}),
        ("chi_pwr_spec_supp_map", '*chi*.dat*', get_plot_supp_map, {'subdir' : 'pwr_spec/', 'pk_type' : 'chi'})
    ]

    # perform all steps, skip step if Exception occurs
    for key, patterns, plot_func, kwargs in all_steps:
        try:
            load_check_plot(stack_info, key, patterns, rerun,
                            skip, plot_func, **kwargs)
        except KeyboardInterrupt:
            raise
        except Exception:
            print_exception()

def get_runs_siminfo(in_dir):
    # get all runs
    files = get_files_in_traverse_dir(in_dir, 'sim_param.json')
    sim_infos =  [struct.SimInfo(a_file) for a_file in files]

    # separate files according to run parameters
    sep_infos = []
    for a_sim_info in sim_infos:
        struct.insert(a_sim_info, sep_infos)
    
    return sep_infos

def count_runs(sep_files):
    # count number of runs
    num_all_runs = num_all_sep_runs = num_sep_runs = 0
    for sep_sim_infos in sep_files:
        num_all_runs += len(sep_sim_infos)
        if len(sep_sim_infos) > 1:
            num_sep_runs += 1
            num_all_sep_runs += len(sep_sim_infos)

    return num_all_runs, num_all_sep_runs, num_sep_runs

def print_runs_info(sep_files, num_all_runs, num_all_sep_runs, num_sep_runs):
    print "There are in total %i different runs, from which %i share the same parameters, constituting %i group(s) eligible for stacking:" % (
        num_all_runs, num_all_sep_runs, num_sep_runs)
    for sep_sim_infos in sep_files:
        print "\n" + sep_sim_infos[0].info_tr()
        for i, a_sim_info in enumerate(sep_sim_infos):
            if i == 10:
                print "\t...and %i more runs" % (len(sep_sim_infos) - 10)
                break
            print "\t" + a_sim_info.dir

def stack_all(in_dir='/home/vrastil/Documents/GIT/Adhesion-Approximation/output/', rerun=None, skip=None, **kwargs):
    # get & count all runs
    sep_files = get_runs_siminfo(in_dir)
    num_all_runs, num_all_sep_runs, num_sep_runs = count_runs(sep_files)
    
    # remove 1-length sep_sim_infos
    sep_files[:] = [x for x in sep_files if len(x) != 1]

    # sort sim_infos
    for sep_sim_infos in sep_files:
        sep_sim_infos.sort(key=lambda x: x.dir)

    # print info about separated files
    print_runs_info(sep_files, num_all_runs, num_all_sep_runs, num_sep_runs)

    # analysis
    for sep_sim_infos in sep_files:
        info = "Stacking group %s" % sep_sim_infos[0].info_tr()
        print '*'*len(info), '\n', info, '\n', '*'*len(info)
        try:
            stack_group(group_sim_infos=sep_sim_infos, rerun=rerun, skip=skip, **kwargs)
        except KeyboardInterrupt:
            print 'Exiting...'
            return
    print '*'*len(info), '\n', 'All groups analyzed!'

# ********************************
# RUN ANALYSIS -- CHI COMPARISON *
# ********************************

def plot_chi_wave_pot(a_file="/home/vrastil/Documents/GIT/FastSim/jobs/output/CHI_run/STACK_512m_512p_1024M_2000b_1e-06Y/stack_info.json",
                      outdir = "/home/vrastil/Documents/GIT/FastSim/report/plots/"):
    a_sim_info = struct.SimInfo(a_file)
    zs = np.linspace(0,5)
    beta = a_sim_info.chi_opt["beta"]
    n = [0.1,0.5,0.7]
    phi = [10**(-5)]
    chi_opt = [{'beta' : beta, 'n' : n_, 'phi' : phi_} for phi_ in phi for n_ in n]
    plot.plot_chi_evol(zs, a_sim_info, chi_opt=chi_opt, out_dir=outdir, show=True)


def get_data_fp_chi_ratio(group):
    data_all = []
    
    data_fp = [np.transpose(np.loadtxt(a_file)) for a_file in group["FP_files"]]
    for files in group["CHI_files"]:
        data_chi = [np.transpose(np.loadtxt(a_file)) for a_file in files]
        data = []
        for data_fp_s, data_chi_s in izip(data_fp, data_chi):
            shape = data_fp_s.shape
            data_ = np.empty(shape)
            data_[0] = data_fp_s[0] # k
            data_[1] = data_chi_s[1] / data_fp_s[1] # Pk
            data_[2] = data_[1]*np.sqrt(np.square(data_chi_s[2] / data_chi_s[1]) + np.square(data_fp_s[2] / data_fp_s[1]))
            data.append(data_)
        data_all.append(data)

    return data_all

def rm_extra_zs(zs_unique, zs_list, other_list):
    i = 0
    for z in zs_list:
        if z not in zs_unique:
            del other_list[i]
        else:
            i += 1

def load_chi_fp_files(group, subdir, patterns):
    zs_fp, files = try_get_zs_files(group["FP"], subdir, patterns)
    group["FP_zs"] = zs_fp
    group["FP_files"] = files
    group["CHI_zs"] = []
    group["CHI_files"] = []
    for chi_info in group["CHI"]:
        zs, files = try_get_zs_files(chi_info, subdir, patterns)
        group["CHI_zs"].append(zs)
        group["CHI_files"].append(files)

    # create zs which are in all subgroups
    zs_unique = set(zs_fp)
    for zs in group["CHI_zs"]:
        zs_unique &= set(zs)
    zs_unique = list(sort_lists(zs_unique)[0])

    # remove extra files from all list
    rm_extra_zs(zs_unique, zs_fp, group["FP_files"])
    for zs_list, other_list in izip(group["CHI_zs"], group["CHI_files"]):
        rm_extra_zs(zs_unique, zs_list, other_list)

    group["FP_zs"] = zs_unique
    group["CHI_zs"] = None

def get_fp_chi_groups(in_dir):
    res = struct.Results(in_dir)
    groups = []
    for a_sim_info in res.get_subfiles(app='FP'):
        Nm = a_sim_info.box_opt["mesh_num"]
        NM = a_sim_info.box_opt["mesh_num_pwr"]
        Np = a_sim_info.box_opt["par_num"]
        L = a_sim_info.box_opt["box_size"]
        chi_infos = res.get_subfiles(Nm=Nm, NM=NM, Np=Np, L=L, app='CHI')
        groups.append({ "FP" : a_sim_info, "CHI" : chi_infos})

    for group in groups:
        load_chi_fp_files(group, 'pwr_spec', '*par*')
    
    return groups

def my_shape(data):
    try:
        data = np.array(data)
        return data.shape
    except ValueError:
        return len(data)
        
def compare_chi_fp(in_dir="/home/vrastil/Documents/GIT/Adhesion-Approximation/output/",
                   out_dir="/home/vrastil/Documents/GIT/Adhesion-Approximation/report/plots/",
                   use_group=None):
    groups = get_fp_chi_groups(in_dir)    
    if use_group is not None:
        groups = [groups[use_group]]
    data_all = [get_data_fp_chi_ratio(group) for group in groups]
    #return data_all
    
    for group, data_g in izip(groups, data_all): # per group
        # struct.SimInfo, zs should be the same for all FP / CHI
        a_sim_info = group["FP"]
        zs = group["FP_zs"]
        phi_s = [si.chi_opt["phi"] for si in group["CHI"]]
        
        plot.plot_chi_fp_map(data_g, zs, a_sim_info)

        data_g_tr = map(list, zip(*data_g)) # transpoose
        
        for lab, data_z in plot.iter_data(zs, [data_g_tr]): # per redshift
            suptitle = "Relative chameleon power spectrum, " + lab
            plot.plot_chi_fp_z(data_z, a_sim_info, phi_s, out_dir=out_dir ,suptitle=suptitle, show=True, save=True)

# *********************************
# RUN ANALYSIS -- CORR COMPARISON *
# *********************************

def load_get_corr(a_file, z=None):
    # get struct.StackInfo
    a_sim_info = struct.StackInfo(stack_info_file=a_file)
    
    # get files for correlation function
    patterns = '*par*.dat *init*.dat'
    subdir = 'pwr_spec/'
    zs, files = try_get_zs_files(a_sim_info, subdir, patterns)
    
    # if z is specified, find nearest-value file
    zs, files = cut_zs_files(zs, files, z)
    
    # get correlation function
    get_corr_func(files, zs, a_sim_info)
    
    # return struct.SimInfo with all loaded data and redshifts
    return a_sim_info, zs
    
def corr_func_ZA_FP(ZA_file="/home/vrastil/Documents/GIT/Adhesion-Approximation/output/FP_run/STACK_512m_512p_1024M_2000b/stack_info.json",
                    FP_file="/home/vrastil/Documents/GIT/Adhesion-Approximation/output/ZA_run/STACK_512m_512p_1024M_2000b/stack_info.json",
                    outdir = "/home/vrastil/Documents/GIT/Adhesion-Approximation/report/plots/",
                    z=1.):
    # load struct.SimInfo and get correlation data
    ZA_ST, z_ZA = load_get_corr(ZA_file, z=z)
    FP_ST, z_FP = load_get_corr(FP_file, z=z_ZA)
    
    # check redshifts
    if z_ZA != z_FP:
        raise IndexError("ZA and FP do not have the same redshift-slice.")
    
    # plot all (one for 'z=None') correlation plots
    r, xi = FP_ST.data["corr_func"]["par"][0]
    extra_data = [
        {'r' : r, 'xi' : xi, 'lab' : FP_ST.app, 'mlt' : 1}
    ]

    plot.plot_corr_func(ZA_ST.data["corr_func"], z_ZA, ZA_ST, out_dir=outdir, save=True, show=True, extra_data=extra_data)
    

def get_pk_broad_k(data_list, sim_infos):
    data_list_new = [[] for i in range(3)]
    
    # sort from lowest k
    data_list, sim_infos = zip(*sorted(zip(data_list, sim_infos), key=lambda x : x[0][0][0]))

    # go through sorted list from, do not go higher values than k_nyqist / 2
    k_last = 0
    for data, a_sim_info in zip(data_list, sim_infos):
        k, Pk, Pk_std = data
        k_max = a_sim_info.k_nyquist["particle"] / 10
        idx = (k >= k_last) & (k < k_max)
        data_list_new[0] += k[idx].tolist()
        data_list_new[1] += Pk[idx].tolist()
        data_list_new[2] += Pk_std[idx].tolist()
        k_last = k_max
        
    # last data get up to k_nyquist
    k_nq = a_sim_info.k_nyquist["particle"]
    idx = (k >= k_last) & (k <= k_nq)
    data_list_new[0] =  np.array(data_list_new[0] + k[idx].tolist())
    data_list_new[1] =  np.array(data_list_new[1] + Pk[idx].tolist())
    data_list_new[2] =  np.array(data_list_new[2] + Pk_std[idx].tolist())
        
    # get new extrapolated Pk, use last a_sim_info for k_nyquist and simulation param.
    sim = a_sim_info.sim
    extrap_pk = get_hybrid_pow_spec_amp(sim, data_list_new, k_nq)
    
    # plot data only to half nyquist frequency (previusly needed for extra_pk)
    idx = (data_list_new[0] <= k_nq/2)
    data_list_new[0] = data_list_new[0][idx]
    data_list_new[1] = data_list_new[1][idx]
    data_list_new[2] = data_list_new[2][idx]
    
    return np.array(data_list_new), extrap_pk

def get_check_pk_broad(stack_infos):
    data_list = []
    zs = []
    for a_sim_info in stack_infos:    
        data_list += a_sim_info.data["pk_data_par"]
        zs += a_sim_info.data["zs"]

    # check of consistency
    if len(set(zs)) != 1:
        raise IndexError("Redshift do not have same value.")
        
    APP = [x.app for x in stack_infos]
    if len(set(APP)) != 1:
        raise IndexError("Different simulations.")

    # get data
    data_list_new, extrap_pk = get_pk_broad_k(data_list, stack_infos)
    return data_list_new, extrap_pk

def get_plot_mlt_pk_broad(stack_infos):
    # sort according to app
    app_all = set([x.app for x in stack_infos])
    groups = {app : [] for app in app_all}
    for a_sim_info in stack_infos:
        app = a_sim_info.app
        groups[app].append(a_sim_info)

    zs, Pk_list_extrap, data_all = [], [], []
    for app in app_all:
        infos = groups[app]
        data_list_new, extrap_pk = get_check_pk_broad(infos)
        groups[app + "_data"] = data_list_new
        groups[app + "_pk"] = extrap_pk
        zs.append(groups[app][0].data["zs"][0])
        Pk_list_extrap.append(extrap_pk["Pk_par"])
        data_all.append(data_list_new)
    
    # plot
    plot.plot_pwr_spec_comparison(data_all, zs, app_all, stack_infos[0].sim.cosmo, save=True, show=True)