"""
'data.py' module serves for loading and saving all necessary data
run analysis of runs
"""

from __future__ import print_function

# system modules
import sys
try:
    from itertools import izip # python 2
except ImportError:
    izip = zip # python 3

# mathematical modules
import numpy as np

# simpy modules
from . import utils as ut
from . import plot
from . import struct
from . import power as pwr
from . import database
from plot import report_dir

###################################
# FIND, SORT, SLICE
###################################

def sort_get_fl_get_z(a_sim_info, subdir, patterns='*.dat'):
    files = ut.get_files_in_traverse_dir(a_sim_info.dir + subdir, patterns)
    return ut.sort_get_z(files, a_sim_info.app)

def try_get_zs_files(a_sim_info, subdir, patterns):
    try:
        zs, files = sort_get_fl_get_z(a_sim_info, subdir, patterns=patterns)
        return list(zs), list(files)
    except ValueError:
        return None, None

def sort_chi_files(files, zs):
    """ separate chi_files (*_chi_*) from given files, return 4 lists """
    files_chi, zs_chi = zip(*[x for x in zip(files, zs) if 'chi' in x[0]])
    files, zs = zip(*[x for x in zip(files, zs) if 'chi' not in x[0]])
    return map(list, [files, zs, files_chi, zs_chi])

###################################
# LOAD DATA FROM SINGLE RUN 
###################################

def has_app_lin_pwr(app):
    return True if app == 'ZA' or app == 'TZA' else False

def load_k_supp(data_list, k_nyquist_par, a_sim_info=None, a=None, pk_type='dens'):
    """
    Divide available k-values into 7 subinterval from
    k_min = 2*PI / L to k_max = 50% k_nyquist_par
    large scale :: k = 1st subinterval
    medium scale :: k = 4rd subinterval
    small scale :: k = 7th subinterval 
    
    For chameleon field divide P(k) values by linear prediction """

    supp = [[[] for x in range(3)] for y in range(3)]
    for i, j in enumerate([0, 3, 6]):
        for ii, data in enumerate(data_list):
            if a_sim_info is not None and a_sim_info.app == "TZA":
                correct_tza_single(a_sim_info, data)
            k = data[0]
            P_diff = data[1]
            if pk_type == 'chi':
                P_diff = -1 + pwr.chi_trans_to_supp(a[ii], k, P_diff, a_sim_info.sim.cosmo, a_sim_info.chi_opt)
            idx = (np.abs(k-0.5*k_nyquist_par)).argmin() / 7
            supp[i][0].append(np.mean(P_diff[j*idx:(j+1)*idx]))
            supp[i][1].append(np.std(P_diff[j*idx:(j+1)*idx]))
        supp[i][2] = [k[j*idx], k[(j+1)*idx]]
    return supp

def get_single_hybrid_pow_spec_amp_w_z(sim, data, z, k_nyquist_par, a=None, fit_lin=False):
    x = pwr.get_hybrid_pow_spec_amp(sim, data, k_nyquist_par, a=a, fit_lin=fit_lin)
    x["z"] = z
    return x

def get_extrap_pk(a_sim_info, data_list, zs):
    if "extrap_pk" not in a_sim_info.data:
        # needed variables
        sim = a_sim_info.sim
        k_nyquist_par = a_sim_info.k_nyquist["particle"]
        fit_lin = has_app_lin_pwr(a_sim_info.app)
        func = lambda data, z : get_single_hybrid_pow_spec_amp_w_z(sim, data, z, k_nyquist_par, fit_lin=fit_lin)

        # get hybrid power spectrum with redshift for each file
        a_sim_info.data["extrap_pk"] = list(map(func, data_list, zs))

        # extract 'Pk_par' and 'z' for convenience
        a_sim_info.data["pk_list"] = [x["Pk_par"] for x in a_sim_info.data["extrap_pk"]]
        a_sim_info.data["zs"] = [x["z"] for x in a_sim_info.data["extrap_pk"]]

def get_pk_nl_amp(a_sim_info):
    # initialize data
    init_data(a_sim_info, get_pk=True)
    load_a_eff(a_sim_info, use_z_eff='Pk')
        
    # needed variables
    sim = a_sim_info.sim
    k_nyquist_par = a_sim_info.k_nyquist["particle"]
    as_eff = 1/(1+a_sim_info.data['eff_time']['Pk']['z_eff'])
    zs = a_sim_info.data["zs"]
    data_all = a_sim_info.data["pk_data_par"]

    # correct for z = 'init'
    if zs[0] == 'init':
        zs = zs[1:]
        data_all = data_all[1:]

    # check lengths
    if not (len(data_all) == len(zs) == len(as_eff)):
        raise IndexError("Data have wrong lengths!")
    
    # get amplitude of non-linear power spectrum with redshift for a_eff
    # fit_lin = has_app_lin_pwr(a_sim_info.app)
    func = lambda a_eff, z, data : get_single_hybrid_pow_spec_amp_w_z(sim, data, z, k_nyquist_par, a=a_eff)
    data_w_amp = list(map(func, as_eff, zs, data_all))        

    # extract amplitude and redshift
    a_sim_info.data["pk_nl_amp"] = {
        'z' : [x['z'] for x in data_w_amp],
        'A' : [x['popt'][0] for x in data_w_amp],
        'A_err' : [x['perr'][0] for x in data_w_amp],
    }

def get_z_eff(a_sim_info, files=None, zs=None, **kwargs):
    use_z_eff = kwargs.get('use_z_eff', None)
    if use_z_eff is None or not get_a_eff(a_sim_info, files, zs, use_z_eff=use_z_eff):
        return zs
    else:
        return a_sim_info.data["eff_time"][use_z_eff] ["z_eff"]

def transform_supp_data_to_z_eff(a_sim_info, use_z_eff='Pk'):
    # get all data
    key = 'pk_supp_input'
    get_supp_map(a_sim_info, key='input')
    data_array = a_sim_info.data[key]['supp']
    zs = a_sim_info.data[key]['zs']
    zs_eff = get_z_eff(a_sim_info, use_z_eff=use_z_eff)
    data_array_new = []
    cosmo = a_sim_info.sim.cosmo

    # transform all data to effective redshift
    for z, z_eff, data in izip(zs, zs_eff, data_array):
        a, a_eff, k, Pk, Pk_std = 1/(1.+z), 1/(1.+z_eff), data[0], data[1], data[2]
        ratio = pwr.growth_factor(a, cosmo) / pwr.growth_factor(a_eff, cosmo)
        Pk = ((Pk + 1) * ratio**2) - 1
        data_array_new.append([k, Pk, Pk_std])

    # correction from initial power spectrum
    data_array_new = np.array(data_array_new)
    Pk_init = data_array_new[0, 1]
    data_array_new[:, 1] -= Pk_init

    return zs_eff, data_array_new
        
def load_plot_pwr(data_list, zs, a_sim_info, **kwargs):
    get_extrap_pk(a_sim_info, data_list, zs)
    plot.plot_pwr_spec(data_list, zs, a_sim_info, a_sim_info.data["pk_list"], **kwargs)

def load_plot_chi_pwr(data_list, zs, a_sim_info, **kwargs):
    plot.plot_chi_pwr_spec(data_list, zs, a_sim_info, **kwargs)

def load_plot_slope(data_list, zs, a_sim_info, **kwargs):
    get_extrap_pk(a_sim_info, data_list, zs)
    plot.plot_slope(data_list, zs, a_sim_info, a_sim_info.data["pk_list"], **kwargs)

def get_key_func(data_list, zs, a_sim_info, key, load=False):
    """ key must be name of callable function placed in power with signature:
        key(sim, Pk=None, z=None, non_lin=False) """

    if key not in a_sim_info.data:
        a_sim_info.data[key] = {
            'par' : [],
            'lin' : [],
            'nl' : [],
            'zs' : []
        }
        gen_func = getattr(pwr, key)
        if load:
            a_sim_info.data[key]["par"] = data_list
            a_sim_info.data[key]["zs"] = zs
        else:
            get_extrap_pk(a_sim_info, data_list, zs)
            zs = a_sim_info.data["zs"]
            for Pk, z in izip(a_sim_info.data["pk_list"], zs):
                data = gen_func(a_sim_info.sim, Pk=Pk)
                if data is not None:
                    a_sim_info.data[key]["par"].append(data)
                    a_sim_info.data[key]["zs"].append(z)

        a_sim_info.data[key]["lin"] = [gen_func(a_sim_info.sim, z=z) for z in a_sim_info.data[key]["zs"]]
        a_sim_info.data[key]["nl"] = [gen_func(a_sim_info.sim, z=z, non_lin=True) for z in a_sim_info.data[key]["zs"]]
        

def get_corr_func(data_list, zs, a_sim_info, load=False):
    get_key_func(data_list, zs, a_sim_info, "corr_func", load=load)

def get_plot_corr(data_list, zs, a_sim_info, load=False, **kwargs):
    get_corr_func(data_list, zs, a_sim_info, load=load)
    plot.plot_corr_func(a_sim_info.data["corr_func"], zs, a_sim_info, **kwargs)

def get_corr_peak(a_sim_info):
    # for each type of measured correlation function
    for key in ('par', 'lin', 'nl'):
        # init data
        init_data(a_sim_info, get_corr=True)

        # for each redshift
        a_sim_info.data["corr_func"][key + '_peak'] = []
        for i, corr in enumerate(a_sim_info.data["corr_func"][key]):
            if corr is not None:
                bao_peak = pwr.get_bao_peak(corr)
                bao_peak["z"] = a_sim_info.data["corr_func"]["zs"][i]
                a_sim_info.data["corr_func"][key + '_peak'].append(bao_peak)

def get_plot_corr_peak(data_list, zs, a_sim_info, load=False, **kwargs):
    get_corr_func(data_list, zs, a_sim_info, load=load)
    get_corr_peak(a_sim_info)
    plot.plot_corr_peak([a_sim_info], **kwargs)
    
def get_sigma_R(data_list, zs, a_sim_info, load=False):
    get_key_func(data_list, zs, a_sim_info, "sigma_R", load=load)

def get_plot_sigma(data_list, zs, a_sim_info, load=False, **kwargs):
    get_sigma_R(data_list, zs, a_sim_info, load=load)
    plot.plot_corr_func(a_sim_info.data["sigma_R"], zs, a_sim_info, is_sigma=True, **kwargs)

def find_nearest_idx(array, value, axis=None):
    # special case for redshifts
    if value == 'init':
        idx = array.index(value)
    else:
        array_ = list(array) # copy
        idx_init = len(array_)
        if 'init' in array:
            idx_init = array.index('init')
            array_ = array_[:idx_init] + array_[idx_init+1:]
        array_ = np.asarray(array_)
        idx = (np.abs(array_ - value)).argmin(axis=axis)
        if idx >= idx_init:
            idx += 1
    return idx

def get_plot_supp(data_list, zs, a_sim_info, pk_type='dens', **kwargs):
    a = [1./(z+1.) for z in zs]
    supp = load_k_supp(data_list, a_sim_info.k_nyquist["particle"], a_sim_info=a_sim_info, a=a, pk_type=pk_type)
    plot.plot_supp_lms(supp, a, a_sim_info, pk_type=pk_type, **kwargs)

def get_supp_map(a_sim_info, key='input'):
    if 'pk_supp_%s' % key not in a_sim_info.data:
        zs, data_array = a_sim_info.get_zs_data('pwr_diff', '*%s*' % key)
        
        a_sim_info.data['pk_supp_%s' % key] = {
            'zs' : zs,
            'supp' : check_data_consistency_diff(data_array)
        }

def correct_tza_single(a_sim_info, data):
    k, P_k = data[0], data[1]
    truncation = pwr.get_truncation(k, a_sim_info.cosmo)
    data[1] = (P_k + 1) * truncation - 1

def correct_tza(a_sim_info, data_array):
    if a_sim_info.app == 'TZA':
        for data in data_array:
            correct_tza_single(a_sim_info, data)

def get_plot_supp_map(data_list, zs, a_sim_info, pk_type='dens', **kwargs):
    data_array = check_data_consistency_diff(data_list)
    plot.plot_pwr_spec_diff_map_from_data(data_array, zs, a_sim_info, ext_title="par", pk_type=pk_type, **kwargs)

def load_plot_pwr_spec_diff(data_list, zs, a_sim_info, **kwargs):
    data_array = check_data_consistency_diff(data_list)
    correct_tza(a_sim_info, data_array)
    plot.plot_pwr_spec_diff_from_data(data_array, zs, a_sim_info, **kwargs)

#
# !!! TO-DO
#
# def split_particle_files(files, zs):
#     files_t = [x for x in files if 'track' in x.split('/')[-1]]
#     files = [x for x in files if 'par_cut' in x.split('/')[-1]]
#     zs = ut.del_duplicate(zs)
#     if len(files_t) != len(files):
#         raise IndexError("Different number of 'particle' and 'track' files.")
#     else:
#         return files, files_t, zs

# def split_plot_par_slice(files, zs, a_sim_info):
#     files, files_t, zs = split_particle_files(files, zs)
#     plot.plot_par_last_slice(files, files_t, zs, a_sim_info)

# def split_plot_par_evol(files, zs, a_sim_info):
#     files, files_t, zs = split_particle_files(files, zs)
#     plot.plot_par_evol(files, files_t, zs, a_sim_info)

def load_plot_dens_histo(data_list, zs, a_sim_info):
    plot.plot_dens_histo(data_list, zs, a_sim_info)

def load_a_eff(a_sim_info, data_list=None, zs=None, use_z_eff='all'):
    # create structure in data
    create_a_eff_struct(a_sim_info)

    # load files
    if data_list is None or zs is None:
        zs, data_list = a_sim_info.get_zs_data('pwr_spec', '*.dat')

    # effective time from power spectrum and density fluctuations
    get_a_eff(a_sim_info, data_list, zs, use_z_eff=use_z_eff)


def load_plot_a_eff(data_list, zs, a_sim_info, **kwargs):
    # load a_eff
    load_a_eff(a_sim_info, data_list=data_list, zs=zs, use_z_eff="Pk")
    # load_a_eff(a_sim_info, data_list=data_list, zs=zs, use_z_eff="sigma_R")
    
    # plot
    # plot.plot_eff_time([a_sim_info], a_eff_type="sigma_R", **kwargs)
    plot.plot_eff_time([a_sim_info], a_eff_type="Pk", **kwargs)

def load_check_plot(a_sim_info, key, patterns, # type: struct.SimInfo, str, str,
                    rerun, skip, plot_func,    # type: List[str], List[str], Callable[List[str], List[str], kwargs],
                    info_str='',               # type: str, str
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
    if a_sim_info.verbose: print('step: %-25s' % (key + ' ' + info_str), end='')
    sys.stdout.flush()
    
    zs, data_list = a_sim_info.get_zs_data(key, patterns)
    if a_sim_info.rerun(rerun, key, skip, zs):
        plot_func(data_list, zs, a_sim_info, **kwargs)
        a_sim_info.done(key)

def cut_zs_list(zs, a_list, z=None):
    if z is not None:
        idx = find_nearest_idx(zs, z)
        zs = [zs[idx]]
        a_list = [a_list[idx]]
    return zs, a_list

def init_data(a_sim_info, z=None, get_pk=False, get_corr=False, get_sigma=False):
    # look for already initialized data
    a_sim_info.load_temp()

    # get files for data
    zs, data_list = a_sim_info.get_zs_data('pwr_spec', '*par*.dat *init*.dat')

    # if z is specified, find nearest-value file
    zs, data_list = cut_zs_list(zs, data_list, z)

    # get extrapolated power spectrum
    if get_pk:
        get_extrap_pk(a_sim_info, data_list, zs)
        a_sim_info.data["pk_data_par"] = data_list

    # get correlation function
    if get_corr:
        get_corr_func(data_list, zs, a_sim_info)

    # get amplitude of density fluctuations
    if get_sigma:
        get_sigma_R(data_list, zs, a_sim_info)

    # save initialized data
    a_sim_info.save_temp()

def get_initialized_StackInfo(a_file, z=None, get_pk=False, get_corr=False, get_sigma=False):
    # get struct.StackInfo
    a_sim_info = struct.StackInfo(stack_info_file=a_file)

    # get data
    init_data(a_sim_info, z=z, get_pk=True, get_corr=get_corr, get_sigma=get_sigma)

    # return initialized object
    return a_sim_info

def reinit_data(sim_infos, get_pk=True, get_corr=True, get_sigma=False):
    for si in sim_infos:
        if get_pk: si.data.pop("extrap_pk", None)
        if get_corr: si.data.pop("corr_func", None)
        if get_sigma: si.data.pop("sigma_R", None)
        init_data(si, get_pk=get_pk, get_corr=get_corr, get_sigma=get_sigma)

# ****************************
# RUN ANALYSIS -- SINGLE RUN *
# ****************************

def analyze_run(a_sim_info, rerun=None, skip=None):
    # try to load data, if exists
    a_sim_info.get_data_from_db()
    a_sim_info.load_temp()

    # Steps to perform -- each entry represents full information needed to perform one step
    # type: Tuple[step_key, data_file_patterns, plot_func, opt_kwargs]
    all_steps = [
        # Power spectrum -- particle, velocity, chameleon, slope
        ("pwr_spec", '*par*.dat *init*.dat', load_plot_pwr, {}),
        ("pwr_spec_chi", '*chi*.dat*', load_plot_chi_pwr, {}),
        ("vel_pwr_spec", '*.dat', load_plot_pwr, {'pk_type' : 'vel'}),
        ("pwr_slope", '*par*.dat*', load_plot_slope, {}),
        # Power spectrum difference -- input, hybrid, particle, velocity, chameleon
        ("pwr_diff", '*par*', load_plot_pwr_spec_diff,
            {'info_str' : '(particle)', 'ext_title' : 'par'}),
        ("pwr_diff_h", '*hybrid*', load_plot_pwr_spec_diff,
            {'info_str' : '(hybrid)', 'ext_title' : 'hybrid'}),
        ("pwr_diff_i", '*input*', load_plot_pwr_spec_diff,
            {'info_str' : '(input)', 'ext_title' : 'input'}),
        ("vel_pwr_diff", '*.dat', load_plot_pwr_spec_diff, {'pk_type' : 'vel'}),
        ("chi_pwr_diff", '*chi*.dat*', load_plot_pwr_spec_diff,
            {'pk_type' : 'chi'}),
        # Power spectrum suppression (includes maps) -- particle, velocity, chameleon,
        ("pwr_spec_supp", '*input*', get_plot_supp, {}),
        ("pwr_spec_supp_map", '*input*', get_plot_supp_map, {}),
        ("vel_pwr_spec_supp", '*.dat', get_plot_supp, {'pk_type' : 'vel'}),
        ("chi_pwr_spec_supp", '*chi*.dat*', get_plot_supp, {'pk_type' : 'chi'}),
        ("chi_pwr_spec_supp_map", '*chi*.dat*', get_plot_supp_map, {'pk_type' : 'chi'}),
        # Correlation function, BAO peak, amplitude of density fluctuations
        ("corr_func", '*par*.dat *init*.dat', get_plot_corr, {}),
        ("bao", '*par*.dat *init*.dat', get_plot_corr_peak, {}),
        ("sigma_R", '*par*.dat *init*.dat', get_plot_sigma, {}),
        # Density distribution
        ("dens_hist", '*.dat', load_plot_dens_histo, {}),
        # Particles -- last slice, evolution
        # !!! TO-DO
        # ("par_slice", 'par*.dat track*.dat', split_plot_par_slice, {}),
        # ("par_ani", 'par*.dat track*.dat', split_plot_par_evol, {}),
        # Density -- two slices, evolution
        ("dens_slice", '*.dat', plot.plot_dens_two_slices, {}),
        ("dens_ani", '*.dat', plot.plot_dens_evol, {}),
        # Effective time
        ("eff_time", '*.dat', load_plot_a_eff, {})
    ]

    # perform all steps, skip step if Exception occurs
    for key, patterns, plot_func, kwargs in all_steps:
        try:
            load_check_plot(a_sim_info, key, patterns, rerun,
                            skip, plot_func, **kwargs)
        except KeyboardInterrupt:
            raise
        except Exception:
            ut.print_exception()

    # save all processed data
    if a_sim_info.verbose: print('step: %-25s' % 'save', end='')
    err_keys = a_sim_info.save_data_to_db()
    if a_sim_info.verbose:
        if err_keys:
            ut.print_warning("[Done]     (some data could not be tranformed into binary)")
        else:
            ut.print_done()
    a_sim_info.save_temp()

def analyze_all(db, collection='data', query=None, rerun=None, skip=None, only=None):
    # filter database
    if query is None:
        query = {}
    query["type"] = 'sim_info'

    # get all ids and sim infos
    sim_infos = []
    for doc_id in db[collection].find(query, {'_id' : 1}):
        sim_infos.append(struct.SimInfo(db, doc_id, collection=collection))

    if only is not None:
        sim_infos = sim_infos[only]

    for a_sim_info in sim_infos:
        ut.print_info('Analyzing run: ' , math_mode=a_sim_info.info_tr(math_mode=True))
        try:
            analyze_run(a_sim_info, rerun=rerun, skip=skip)
        except KeyboardInterrupt:
            print('Exiting...')
            return
    ut.print_info_end()

# ******************************
# LOAD DATA FROM MULTIPLE RUNS *
# ******************************

def load_data_for_stack(stack_info, key, patterns):
    # load everything
    all_data_k = None
    all_data_Pk = None
    all_zs = None
    for doc_id in stack_info:
        # load data for ONE run
        zs, data_list = stack_info.get_zs_data(key, patterns, doc_id=doc_id)

        # if data are missing for this run
        if zs is None:
            continue
        else:
            all_zs = zs

        # create lists, only for the first run (the rest are assumed to have the same output)
        if all_data_k is None:
            all_data_k = [[] for x in range(len(zs))]
        if all_data_Pk is None:
            all_data_Pk = [[] for x in range(len(zs))]

        # load k, Pk
        for i, data in enumerate(data_list):
            k, P_k = data[0], data[1]
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
        print('\t\tData in data_list have different shapes. Trying to cut...')
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
            print("\t\tDeleted %i excess values." % (del_num))
        data_array = np.array(data_list)
    return data_array

def stack_files(stack_info, key, a_file):
    """ assume data in files are k = data[:, 0], Pk = data[:, 1]
    return tuple (zs, data_list) where each entry in 'data_list' corresponds to entry in 'zs' list
    data[:, 0] = k, data[:, 1] = mean(Pk), data[:, 2] = std(Pk)
    """
    # load everything
    zs, all_data_k, all_data_Pk = load_data_for_stack(stack_info, key, a_file)

    # if data are missing for all runs
    if zs is None:
        return None, None

    # chceck lengths of lists, delete excess (in case outputs of simulations encountered any errors)
    check_data_consistency(all_data_k, all_data_Pk)

    # compute means, std
    data_list = []
    for k, Pk in izip(all_data_k, all_data_Pk):
        data_list.append([
            np.mean(k, axis=0),
            np.mean(Pk, axis=0),
            np.std(Pk, axis=0)
        ])

    return zs, data_list

def create_a_eff_struct(a_sim_info):
    if "eff_time" not in a_sim_info.data:
        a_sim_info.data["eff_time"] = {"Pk" : {}, "Pk_nl" : {}, "sigma_R" : {}}

def get_a_eff_from_dens_fluct(a_sim_info):
    """ fit data [r, sigma] to """
    # check zs
    zs = a_sim_info.data["sigma_R"]["zs"]
    idx = 1 if zs[0] == 'init' else 0
    a = [1./(1+z) for z in zs if z != 'init']

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

    # effective redshift
    D_eff = D_eff_ratio * pwr.growth_factor(a, a_sim_info.sim.cosmo)
    a_eff = pwr.get_a_from_growth(D_eff, a_sim_info.sim.cosmo)

    # create structure in data
    create_a_eff_struct(a_sim_info)

    # store
    a_sim_info.data["eff_time"]["sigma_R"] = {
        'a' : a,
        'z_eff' : 1./a_eff - 1,
        'a_err' : D_eff_std,
        'D_eff_ratio' : D_eff_ratio,
        'D_eff' : D_eff
    }


def get_a_eff_from_Pk(stack_info):
    # go through all extrapolated Pk
    a, A = map(np.array, zip(*[ # extract back, store as np.array
        (1/(Pk['z'] + 1), Pk['Pk_par'].A_low) # a, Extrap_Pk_Nl
        for Pk in stack_info.data["extrap_pk"] if Pk["z"] != 'init']))
    
    # get a_eff from amlitude of linear power spectrum
    a_eff = pwr.get_a_from_A(stack_info.sim.cosmo, A)
    
    # derived variables
    D = pwr.growth_factor(a, stack_info.sim.cosmo)
    D_eff = pwr.growth_factor(a_eff, stack_info.sim.cosmo)

    # create structure in data
    create_a_eff_struct(stack_info)

    # store
    stack_info.data["eff_time"]["Pk"] = {
        'a' : a,
        'z_eff' : 1./a_eff - 1,
        'a_err' : 0,
        'D_eff_ratio' : D_eff / D,
        'D_eff' : D_eff
    }

def get_a_eff_from_Pk_nl(stack_info):
    # go through all extrapolated Pk
    a, popt, perr = map(np.array, zip(*[ # extract back, store as np.array
        (1/(Pk['z'] + 1), Pk['popt'], Pk['perr']) # a, popt, perr
        for Pk in stack_info.data["extrap_pk"] if Pk["z"] != 'init']))
    

    # ZA and TZA do not have non-linear fit
    if None in popt:
        return

    # derived variables
    a_eff = popt[:,0]
    D = pwr.growth_factor(a, stack_info.sim.cosmo)
    D_eff = pwr.growth_factor(a_eff, stack_info.sim.cosmo)

    # create structure in data
    create_a_eff_struct(stack_info)

    # store
    stack_info.data["eff_time"]["Pk_nl"] = {
        'a' : a,
        'z_eff' : 1./a_eff - 1,
        'a_err' : perr[:,0],
        'D_eff_ratio' : D_eff / D
    }

def get_a_eff(a_sim_info, files, zs, use_z_eff='Pk'):
    success = False

    # effective time from power spectrum
    if use_z_eff == 'Pk' or use_z_eff == 'all':
        get_extrap_pk(a_sim_info, files, zs)
        get_a_eff_from_Pk(a_sim_info)
        success = True

    # effective time from non-linear power spectrum
    if use_z_eff == 'Pk_nl' or use_z_eff == 'all':
        get_extrap_pk(a_sim_info, files, zs)
        get_a_eff_from_Pk_nl(a_sim_info)
        success = True
    
    # effective time from density fluctuations
    if use_z_eff == 'sigma_R' or use_z_eff == 'all':
        get_sigma_R(files, zs, a_sim_info)
        get_a_eff_from_dens_fluct(a_sim_info)
        success = True

    return success

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
                    info_str=''                 # type: str
                   ):
    """bring structure of stacking together:
    1) check if the step should be performed (rerun / skip)
    2) stack files through function 'stack_files'
    3) save stacked files
       - out_dir = stack_info.dir + subdir
       - fname = fname + "_%s_%s" % (app, z_str)
    4) write info about done step into stack_info
    """
    if stack_info.verbose: print('step: %-25s' % (key + ' ' + info_str), end='')
    sys.stdout.flush()
    # check only rerun / key / skip -- no files loading
    if stack_info.rerun(rerun, key, skip, True):
        zs, data_list = stack_files(stack_info, key, patterns)

        # check again if we loaded any data
        if stack_info.rerun(rerun, key, skip, zs):
            # save all stacked data into the files
            for z, data in zip(zs, data_list):
                z_str = 'init.dat' if z == 'init' else 'z%.2f.dat' % z
                fname_ = stack_info.dir + struct.RESULTS_DIRS[key] + '/' + fname + "_%s_%s" % (stack_info.app, z_str)
                np.savetxt(fname_, np.transpose(data), fmt='%.6e', header=header)

            # save all stacked data into the database
            stack_info.save_zs_data(key, zs, data_list, fname)
            stack_info.done(key)

# **********************************
# RUN ANALYSIS -- STACKING OF RUNS *
# **********************************

def stack_group(db, sep_id, collection='data', rerun=None, skip=None, verbose=False, **kwargs):
    # load & save all info about stack
    stack_info = struct.StackInfo(db, sep_id, collection=collection, verbose=verbose, **kwargs)

    # load, stack & save all files
    all_files_steps = [
        # Power spectrum
        ("pwr_spec_files", '*par*.dat *init*.dat', HEADER_PWR, 'pwr_spec_par', {}),
        ("pwr_spec_chi_files", '*chi*.dat*', HEADER_PWR_CHI, 'pwr_spec_chi', {}),
        # Power spectrum difference -- input, hybrid, particle
        ("pwr_diff_files", '*par*', HEADER_PWR_DIFF_PAR, 'pwr_spec_diff_par', {'info_str' : '(particle)'}),
        ("pwr_diff_files_h", '*hybrid*', HEADER_PWR_DIFF_HYBRID, 'pwr_spec_diff_hybrid', {'info_str' : '(hybrid)'}),
        ("pwr_diff_files_i", '*input*', HEADER_PWR_DIFF_INPUT, 'pwr_spec_diff_input', {'info_str' : '(input)'})
    ]
    for key, patterns, header, fname, kwargs in all_files_steps:
        try:
            load_stack_save(stack_info, key, patterns, rerun,
                            skip, header, fname, **kwargs)
        except KeyboardInterrupt:
            raise
        except Exception:
            ut.print_exception()

    return stack_info

def get_runs_siminfo(in_dir):
    # get all runs
    files = ut.get_files_in_traverse_dir(in_dir, 'sim_param.json')
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
    print("There are in total %i different runs, from which %i share the same parameters, constituting %i group(s) eligible for stacking:" % (
        num_all_runs, num_all_sep_runs, num_sep_runs))
    for sep_sim_infos in sep_files:
        print("\n" + sep_sim_infos[0].info_tr())
        for i, a_sim_info in enumerate(sep_sim_infos):
            if i == 10:
                print("\t...and %i more runs" % (len(sep_sim_infos) - 10))
                break
            print("\t" + a_sim_info.dir)

def stack_all(db, collection='data', query=None, rerun=None, skip=None, verbose=False, **kwargs):
    # separate all runs
    all_sep_id = database.get_separated_ids(db, query=query)
    
    # remove 1-length sep_sim_infos
    all_sep_id[:] = [x for x in all_sep_id if len(x) != 1]

    # stack runs and performe analysis
    stack_infos = []
    for i, sep_id in enumerate(all_sep_id):
        if verbose:
            ut.print_info('Analyzing run ', math_mode=struct.SimInfo(db, sep_id[0]['_id']).info_tr(math_mode=True))
        else:
            print("\rStacking group %i/%i" % (i + 1, len(all_sep_id)), end="")
            sys.stdout.flush()
        try:
            stack_info = stack_group(db, sep_id, collection=collection, rerun=rerun, skip=skip, verbose=verbose, **kwargs)
            analyze_run(stack_info, rerun=rerun, skip=skip)
            stack_infos.append(stack_info)
        except KeyboardInterrupt:
            print('Exiting...')

    if verbose:
        ut.print_info_end()
    else:
        print('\n')

    return stack_infos

# ********************************
# RUN ANALYSIS -- CHI COMPARISON *
# ********************************

def plot_chi_wave_pot(db, collection='data', outdir=report_dir,
                      n=None, phi=None, zs=None, save=True, show=True, k_scr=False):
    doc_id = db.data.find_one({'app' : 'CHI'}, {'_id' : 1})
    a_sim_info = struct.SimInfo(db, doc_id)

    # parameters of the plot (if not given)
    if zs is None:
        zs = np.linspace(0,5)
    beta = a_sim_info.chi_opt["beta"]
    if n is None:
        n = [0.1,0.5,0.7]
    if phi is None:
        phi = [10**(-5)]

    # all chameleon parameters to go through
    chi_opt = [{'beta' : beta, 'n' : n_, 'phi' : phi_} for phi_ in phi for n_ in n]

    # plot
    plot.plot_chi_evol(zs, a_sim_info, chi_opt=chi_opt, out_dir=outdir, save=save, show=show, k_scr=k_scr)


def get_data_fp_chi_ratio(group, z=None):
    data_all = []

    if z is not None:
        zs = group["FP_zs"]
        idx = find_nearest_idx(zs, z, axis=None)
        cut = slice(idx, idx+1)
    else:
        cut = slice(None)

    data_fp = [np.transpose(np.loadtxt(a_file)) for a_file in group["FP_files"][cut]]
    for files in group["CHI_files"]:
        data_chi = [np.transpose(np.loadtxt(a_file)) for a_file in files[cut]]
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
    zs_unique = list(ut.sort_lists(zs_unique)[0])

    # remove extra files from all list
    rm_extra_zs(zs_unique, zs_fp, group["FP_files"])
    for zs_list, other_list in izip(group["CHI_zs"], group["CHI_files"]):
        rm_extra_zs(zs_unique, zs_list, other_list)

    group["FP_zs"] = zs_unique
    del group["CHI_zs"]

def get_fp_chi_groups(in_dir,Nm=0, NM=0, Np=0, L=0, n=0, phi=0):
    res = struct.Results(in_dir)
    groups = []
    for a_sim_info in res.get_subfiles(app='FP', Nm=Nm, NM=NM, Np=Np, L=L):
        Nm = a_sim_info.box_opt["mesh_num"]
        NM = a_sim_info.box_opt["mesh_num_pwr"]
        Np = a_sim_info.box_opt["par_num"]
        L = a_sim_info.box_opt["box_size"]
        chi_infos = res.get_subfiles(Nm=Nm, NM=NM, Np=Np, L=L, app='CHI', n=n, phi=phi)
        if chi_infos:
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
        


def compare_chi_fp(in_dir="/home/michal/Documents/GIT/FastSim/output/", out_dir=report_dir, use_group=None,
                   z=None, Nm=0, NM=0, Np=0, L=0, n=0, phi=0):
    groups = get_fp_chi_groups(in_dir, Nm=Nm, Np=Np, L=L, n=n, phi=phi)
    if use_group is not None:
        groups = [groups[use_group]]

    for group in groups: # per group
        # load data
        data_g = get_data_fp_chi_ratio(group, z=z)
        
        # struct.SimInfo, zs should be the same for all FP / CHI
        a_sim_info = group["FP"]
        zs = group["FP_zs"] if z is None else [z]

        # transpoose first and second dimension
        data_g = map(list, zip(*data_g))
        
        for lab, data_z in plot.iter_data(zs, [data_g]): # per redshift
            # sort from lowest screening potential and xhameleon exponent
            sim_infos, data_z = zip(*sorted(sorted(sorted(zip(group["CHI"], data_z),
                                key=lambda x : x[0].chi_opt['linear']),
                                key=lambda x : x[0].chi_opt['n']),
                                key=lambda x : x[0].chi_opt['phi'])
                            )

            # get labels
            labels = plot.get_chi_labels(sim_infos)
            suptitle = "Relative chameleon power spectrum, " + lab
            ut.print_info(suptitle)
            
            plot.plot_chi_fp_z(data_z, a_sim_info, labels, out_dir=out_dir ,suptitle=suptitle, show=True, save=True)

def compare_chi_fp_map(chi_info, fp_info, in_dir="/home/michal/Documents/GIT/FastSim/output/", out_dir=report_dir):
    # check we have same parameters
    for key in ("mesh_num", "mesh_num_pwr", "par_num", "box_size"):
        if chi_info.box_opt[key] != fp_info.box_opt[key]:
            raise IndexError("CHI and FPA files do not have the same parameters.")

    # load data
    group = { "FP" : fp_info, "CHI" : [chi_info]}
    load_chi_fp_files(group, 'pwr_spec', '*par*')
    data = get_data_fp_chi_ratio(group)[0] # we have only one set of chi-parameters

    # get effective redshift
    load_a_eff(fp_info, use_z_eff='Pk')
    zs_eff = fp_info.data["eff_time"]['Pk']["z_eff"]

    # get rid if 'init' redshift and different ones
    if group["FP_zs"][0] == 'init':
        data = np.array(data[1:])

    # check lengths
    if len(zs_eff) != len(data):
        raise IndexError("Files do not have the same time slices.")

    # plot map
    plot.plot_chi_fp_map(data, zs_eff, chi_info, out_dir=out_dir, save=True, show=True)


def compare_chi_res(in_dir, out_dir, n=0.5, phi=1e-5, z=0):
    # load all data
    groups = get_fp_chi_groups(in_dir, n=n, phi=phi)
    data_all = [get_data_fp_chi_ratio(group, z=z) for group in groups]
    
    # check we have same chameleon parameters
    phi, n = set(), set()
    for group in groups:
        for si in group["CHI"]:
            phi.add(si.chi_opt["phi"])
            n.add(si.chi_opt["n"])
    if len(phi) != len(n) != 1:
            raise IndexError("CHI files do not have the same chameleon parameters.")
    phi = phi.pop()
    n = n.pop()

    # # check we loaded only one redshift, remove needless dimension
    for i, data in enumerate(data_all):
        data = data_all[i] = np.array(data)
        if data.shape[1] != 1:
            raise IndexError("Data for different redshifts when expecting only one!")
        new_shape = data.shape[:1] + data.shape[2:]
        data_all[i] = data.reshape(new_shape)

    print('Phi = ', phi, "\tn = ", n)
    sim_infos = [group["CHI"] for group in groups]

    # sort from lowest to highest resolution
    sim_infos, data_all = zip(*sorted(zip(sim_infos, data_all), key=lambda x : x[0][0].k_nyquist["potential"]))

    plot.plot_chi_fp_res(data_all, sim_infos, out_dir=out_dir, show=True, save=True)

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
    zs, files = cut_zs_list(zs, files, z)
    
    # get correlation function
    get_corr_func(files, zs, a_sim_info)
    
    # return struct.SimInfo with all loaded data and redshifts
    return a_sim_info, zs
    
def corr_func_comp_plot(files=None, sim_infos=None, outdir=report_dir, z=1., bao_peak=True):

    extra_data = []
    zs = None

    # load struct.SimInfo and get correlation data
    if sim_infos is None:
        sim_infos = [load_get_corr(a_file)[0] for a_file in files]

    # get data, check redshift
    for sim_info in sim_infos:
        # init data
        init_data(sim_info, get_corr=True)
        load_a_eff(sim_info, use_z_eff='Pk')

        zs_ = sim_info.data["corr_func"]["zs"]
        idx = find_nearest_idx(zs_, z)

        # save needed values
        r, xi = sim_info.data["corr_func"]["par"][idx]
        idx_eff = idx - 1 if 'init' in zs_ else idx
        z_eff = sim_info.data["eff_time"]["Pk"]["z_eff"][idx_eff]
        extra_data.append({'r' : r, 'xi' : xi, 'lab' : sim_info.app, 'mlt' : 1, 'z_eff' : z_eff})
    
        # check redshifts
        if zs is None:
            zs = zs_[idx]
        elif zs != zs_[idx]:
            raise IndexError("Files do not have the same redshift-slices.")

    # plot non-linear BAO peak
    if bao_peak:
        get_corr_peak(sim_info)
        peak_loc = sim_info.data["corr_func"]["nl_peak"][idx]["popt"][1]
    else:
        peak_loc = None

    data = {
        'par' : [sim_info.data["corr_func"]["par"][idx]],
        'lin': [sim_info.data["corr_func"]["lin"][idx]],
        'nl' : [sim_info.data["corr_func"]["nl"][idx]]
    }
    use_z_eff = {
        'z' : z_eff,
        'sim' : sim_info.sim
    }

    # plot simple correlation function and ratio
    plot.plot_corr_func(data, [zs], sim_info, out_dir=outdir, save=True, show=True, extra_data=extra_data[:-1], peak_loc=peak_loc, use_z_eff=use_z_eff)

def corr_func_comp_plot_peak(files=None, sim_infos=None, outdir=report_dir, plot_all=True):
    # load struct.SimInfo and get correlation data
    if sim_infos is None:
        sim_infos = [load_get_corr(a_file)[0] for a_file in files]

    # sort from lowest screening potential and chameleon exponent
    non_chi_sim_infos = [x for x in sim_infos if x.app != 'CHI']
    chi_sim_infos = [x for x in sim_infos if x.app == 'CHI']
    chi_sim_infos = sorted(sorted(chi_sim_infos,
                        key=lambda x : x.chi_opt['n']),
                        key=lambda x : x.chi_opt['phi'])
    sim_infos = non_chi_sim_infos + chi_sim_infos

    # get data, check redshift
    use_z_eff = []
    for sim_info in sim_infos:
        init_data(sim_info, get_corr=True)
        get_a_eff_from_Pk(sim_info)
        get_corr_peak(sim_info)

        use_z_eff.append({
            'z' : sim_info.data["eff_time"]["Pk"]["z_eff"],
            'sim' : sim_info.sim
        })

    # plot bao peak and location
    if plot_all:
        plot.plot_corr_peak(sim_infos, out_dir=outdir, save=True, show=True, use_z_eff=use_z_eff)
    else:
        plot.plot_corr_peak(sim_infos, out_dir=outdir, save=True, show=True, use_z_eff=use_z_eff, plot_loc=True, plot_amp=False, plot_width=False)
        plot.plot_corr_peak(sim_infos, out_dir=outdir, save=True, show=True, use_z_eff=use_z_eff, plot_loc=False, plot_amp=True, plot_width=False)
        plot.plot_corr_peak(sim_infos, out_dir=outdir, save=True, show=True, use_z_eff=use_z_eff, plot_loc=False, plot_amp=False, plot_width=True)


def get_pk_broad_k(data_list, sim_infos, get_extrap_pk=True, cutoff_high=4.3, lim_kmax=None):
    data_list_new = [[] for _ in range(3)]
    
    # sort from lowest k
    data_list, sim_infos = zip(*sorted(zip(data_list, sim_infos), key=lambda x : x[0][0][0]))

    # go through sorted list, do not go higher values than k_nyqist / cutoff_high
    k_last = 0
    for data, a_sim_info in zip(data_list, sim_infos):
        k, Pk, Pk_std = data
        k_max = a_sim_info.k_nyquist["particle"] / cutoff_high
        if lim_kmax is not None:
            k_max = np.min((k_max, lim_kmax))
        idx = (k >= k_last) & (k < k_max)
        data_list_new[0] += k[idx].tolist()
        data_list_new[1] += Pk[idx].tolist()
        data_list_new[2] += Pk_std[idx].tolist()
        k_last = k_max
        
    # last data get up to k_nyquist
    k_nq = a_sim_info.k_nyquist["particle"]
    if lim_kmax is not None:
        k_nq = np.min((k_nq, lim_kmax))
    idx = (k >= k_last) & (k <= k_nq)
    data_list_new[0] =  np.array(data_list_new[0] + k[idx].tolist())
    data_list_new[1] =  np.array(data_list_new[1] + Pk[idx].tolist())
    data_list_new[2] =  np.array(data_list_new[2] + Pk_std[idx].tolist())
        
    # get new extrapolated Pk, use last a_sim_info for k_nyquist and simulation param.
    sim = a_sim_info.sim
    fit_lin = has_app_lin_pwr(sim_infos[0].app)
    extrap_pk = pwr.get_hybrid_pow_spec_amp(sim, data_list_new, k_nq, fit_lin=fit_lin) if get_extrap_pk else None
    
    # plot data only to half nyquist frequency (previusly needed for extra_pk)
    idx = (data_list_new[0] <= k_nq/2)
    data_list_new[0] = data_list_new[0][idx]
    data_list_new[1] = data_list_new[1][idx]
    data_list_new[2] = data_list_new[2][idx]
    
    return np.array(data_list_new), extrap_pk

def get_check_pk_broad(stack_infos, idx, data_key="pk", get_extrap_pk=True, cutoff_high=4.48, lim_kmax=None):
    data_list = []
    zs = []
    for a_sim_info in stack_infos:
        if data_key == 'pk':
            data = a_sim_info.data["pk_data_par"][idx]
            z = a_sim_info.data["zs"][idx]
        elif data_key == 'supp':
            data = a_sim_info.data["pk_supp_input"]["supp"][idx]
            z = a_sim_info.data["pk_supp_input"]["zs"][idx]
        elif data_key == 'supp_z_eff':
            _, data_list_tmp = transform_supp_data_to_z_eff(a_sim_info)
            data = data_list_tmp[idx]
            z = a_sim_info.data["pk_supp_input"]["zs"][idx]

        data_list.append(data)
        zs.append(z)

    # check of consistency
    if len(set(zs)) != 1:
        raise IndexError("Redshift do not have same value.")
        
    APP = [x.app for x in stack_infos]
    if len(set(APP)) != 1:
        raise IndexError("Different simulations.")

    # get data
    data_list_new, extrap_pk = get_pk_broad_k(data_list, stack_infos, get_extrap_pk=get_extrap_pk, cutoff_high=cutoff_high, lim_kmax=lim_kmax)
    return data_list_new, extrap_pk

def get_check_pk_diff_broad(stack_infos, cutoff_high=8, lim_kmax=2):
    data_lists = []
    zs = stack_infos[0].data["pk_supp_input"]["zs"]

    for idx in range(len(zs)):
        data_list_new, _ = get_check_pk_broad(stack_infos, idx, data_key='supp_z_eff', get_extrap_pk=False, cutoff_high=cutoff_high, lim_kmax=lim_kmax)
        data_lists.append(data_list_new)
    return np.array(data_lists), zs

def get_init_group(db, query, collection='data'):
    # group by app
    group = {'_id' : {'app' : '$app'}, 'ids' : {'$addToSet' : '$_id'}}
    pipeline = [
        {"$match" : query },
        {"$group" : group }
    ]

    # init all, make list
    groups = {}
    for doc in db[collection].aggregate(pipeline):
        app = doc['_id']['app']
        doc_ids = doc['ids']
        groups[app] = []
        for doc_id in doc_ids:
            stack_info = struct.StackInfo(db, doc_id, collection=collection)
            init_data(stack_info, get_pk=True)
            groups[app].append(stack_info)

    return groups

def get_plot_mlt_pk_broad(db, query=None, collection='data', out_dir='auto', z=0):
    # default to non-CHI StackInfo with NM = 1024
    if query is None:
        query = {'app' : {"$ne" : 'CHI'}, 'type' : 'stack_info', 'box_opt.mesh_num_pwr' : 1024}

    # get initializied groups
    groups = get_init_group(db, query, collection)

    # get all data
    zs, Pk_list_extrap, data_all = [], [], []
    for group in groups.values():
        # get idx of nearest redshift in zs
        zs_ = group[0].data["zs"]
        idx = find_nearest_idx(zs_, z)
        data_list_new, extrap_pk = get_check_pk_broad(group, idx, lim_kmax=2)
        zs.append(zs_[idx])
        Pk_list_extrap.append(extrap_pk["Pk_par"])
        data_all.append(data_list_new)
    
    # plot
    cosmo = group[0].sim.cosmo
    plot.plot_pwr_spec_comparison(Pk_list_extrap, data_all, zs, groups.keys(), cosmo, out_dir=out_dir, save=True, show=True)
    plot.plot_pwr_spec_comparison_ratio_nl(Pk_list_extrap, data_all, zs, groups.keys(), cosmo, out_dir=out_dir, save=True, show=True)

def get_plot_mlt_pk_diff_broad(stack_infos, plot_diff=True, out_dir='auto'):
    # sort according to app
    groups = get_init_group(stack_infos)

    # get all data
    for group in groups.values():
        for a_sim_info in group:
            get_supp_map(a_sim_info)
        data_array, zs = get_check_pk_diff_broad(group, cutoff_high=16, lim_kmax=2)
        correct_tza(group[0], data_array)
        
        # plot
        ut.print_info(group[0].app)
        plot.plot_pwr_spec_diff_map_from_data(
            data_array, zs, a_sim_info, out_dir=out_dir, add_app=True, 
            show_nyquist=False, save=True, show=True, shading='gouraud')
        
        if plot_diff:
            plot.plot_pwr_spec_diff_from_data(
                data_array, zs, a_sim_info, out_dir=out_dir, add_app=True, 
                show_nyquist=False, show_scales=False, save=True, show=True)