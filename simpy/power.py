
"""
'power.py' module serves mainly for interacting with C++ library fastsim.py

 - translate all power spectra, growth functions, correlations functions, etc.
   into C++ functions for speed
 - handles numpy arrays
 - cosmo == C++ class Cosmo_Param, accessible through SimInfo.sim.cosmo
 - FTYPE_t=[float, double, long double]
"""

import numpy as np
from scipy.optimize import brentq, curve_fit, minimize_scalar
from . import fastsim as fs

class Non_Linear_Cosmo(object):
    # private variables and methods
    _cosmo_emu = None
    _cosmo_halofit = None
    _sim_emu = None
    _sim_halofit = None

    @staticmethod
    def _init_cosmo(cosmo):
        if Non_Linear_Cosmo._cosmo_emu is None:
            Non_Linear_Cosmo._cosmo_emu = Non_Linear_Cosmo._copy_cosmo(cosmo, fs.ccl_emu, transfer_function=fs.ccl_emulator)

        if Non_Linear_Cosmo._cosmo_halofit is None:
            Non_Linear_Cosmo._cosmo_halofit = Non_Linear_Cosmo._copy_cosmo(cosmo, fs.ccl_halofit)

    @staticmethod
    def _init_sim(sim):
        Non_Linear_Cosmo._init_cosmo(sim.cosmo)

        if Non_Linear_Cosmo._sim_emu is None:
            Non_Linear_Cosmo._sim_emu = Non_Linear_Cosmo._copy_sim(sim, Non_Linear_Cosmo._cosmo_emu)

        if Non_Linear_Cosmo._sim_halofit is None:
            Non_Linear_Cosmo._sim_halofit = Non_Linear_Cosmo._copy_sim(sim, Non_Linear_Cosmo._cosmo_halofit)

    @staticmethod
    def _copy_cosmo(cosmo_from, matter_power_spectrum_method, transfer_function=None):
        # create empty Cosmo_Param
        cosmo_to = fs.Cosmo_Param()

        # copy basic parameterz
        cosmo_to.sigma8 = cosmo_from.sigma8
        cosmo_to.ns = cosmo_from.ns
        cosmo_to.k2_G = cosmo_from.k2_G
        cosmo_to.Omega_m = cosmo_from.Omega_m
        cosmo_to.Omega_b = cosmo_from.Omega_b
        cosmo_to.H0 = cosmo_from.H0
        cosmo_to.truncated_pk = cosmo_from.truncated_pk

        # copy ccl methods
        cosmo_to.config.baryons_power_spectrum_method = cosmo_from.config.baryons_power_spectrum_method
        cosmo_to.config.mass_function_method = cosmo_from.config.mass_function_method
        # -> transfer function and matter power spectrum is different
        if transfer_function is None:
            cosmo_to.config.transfer_function_method = cosmo_from.config.transfer_function_method
        else:
            cosmo_to.config.transfer_function_method = transfer_function

        cosmo_to.config.matter_power_spectrum_method = matter_power_spectrum_method

        # initialize Cosmo_Param
        cosmo_to.init()

        return cosmo_to

    @staticmethod
    def _copy_sim(sim_from, cosmo):
        sim_to = fs.Sim_Param()
        sim_to.cosmo = cosmo
        sim_to.box_opt = sim_from.box_opt
        sim_to.integ_opt = sim_from.integ_opt
        sim_to.out_opt = sim_from.out_opt
        sim_to.comp_app = sim_from.comp_app
        sim_to.app_opt = sim_from.app_opt
        sim_to.run_opt = sim_from.run_opt
        sim_to.other_par = sim_from.other_par
        sim_to.chi_opt = sim_from.chi_opt
        sim_to.test_opt = sim_from.test_opt
        return sim_to

    @staticmethod
    def non_lin_pow_spec(a, k, cosmo):
        """ return ndarray of nonlinear power spectrum """
        # initialize emulator cosmology (if not done already)
        Non_Linear_Cosmo._init_cosmo(cosmo)

        # for z < 2 use emulator, halofit otherwise
        if a < 0.3:
            cosmo_ = Non_Linear_Cosmo._cosmo_halofit
        else:
            cosmo_ = Non_Linear_Cosmo._cosmo_emu

        # call non-linear power spectrum
        k = np.array(k)
        if k.shape:
            return np.array([fs.non_lin_pow_spec(a, k_, cosmo_)  for k_ in k])
        else:
            return fs.non_lin_pow_spec(a, np.asscalar(k), cosmo_)

    @staticmethod
    def gen_func(sim, a, data, gen_func):
        # initialize emulator cosmology (if not done already)
        Non_Linear_Cosmo._init_sim(sim)

        # for z < 2 use emulator, halofit otherwise
        if a < 1./3.:
            sim_ = Non_Linear_Cosmo._sim_halofit
        else:
            sim_ = Non_Linear_Cosmo._sim_emu

        # call non-linear gen_func
        gen_func(sim_, a, data)

    @staticmethod
    def non_lin_corr_func(sim, a, data):
        """ return non-linear correlation function """
        Non_Linear_Cosmo.gen_func(sim, a, data, fs.gen_corr_func_binned_gsl_qawf_nl)

    @staticmethod
    def non_lin_sigma_func(sim, a, data):
        """ return non-linear amplitude of density fluctuations """
        Non_Linear_Cosmo.gen_func(sim, a, data, fs.gen_sigma_func_binned_gsl_qawf_nl)


def get_a_init_from_zs(zs):
    """ from list of redshifts returns initial scale factor, i.e. value after 'init' """
    for z in zs:
        if z != 'init':
            return 1/(1.+z)

def get_a_fom_zs(zs):
    try:
        iter(zs)
    except TypeError:
        return 1./(1+zs)
    else:
        a = [1./(z + 1) for z in zs if z != 'init']
        return np.array(a)

def get_z_from_a(a):
    try:
        iter(a)
    except TypeError:
        return 1./a - 1
    else:
        zs = [1./a_ - 1 for a_ in a]
        return np.array(zs)

def get_a_from_A(cosmo, A):
    """ return scale factor a at which amplitude of linear power spectrum is A (normalize as A=1 at a=1) """
    # 'f = 0' <=> A = D^2 (linear power grows as D^2)
    A = np.array(A)
    if A.shape:
        a_eff = []
        for A_ in A:
            f = lambda a : A_ - fs.growth_factor(a, cosmo)**2
            a_eff.append(brentq(f, 0, 2))
        return np.array(a_eff)
    else:
        f = lambda a : A - fs.growth_factor(a, cosmo)**2
        return brentq(f, 0, 2)

def get_ndarray(Data_Vec):
    """ copy C++ class Data_Vec<FTYPE_t, N> into numpy array """
    dim = Data_Vec.dim()
    data = [[x for x in Data_Vec[i]] for i in range(dim)]
    return np.array(data)

def get_Data_vec(data):
    """ copy 2D data 'dim x size' into C++ class Data_Vec<FTYPE_t, dim> """
    dim = len(data)
    size = len(data[0])
    if dim == 2:
        Data_Vec = fs.Data_Vec_2(size)
    elif dim == 3:
        Data_Vec = fs.Data_Vec_3(size)
    else:
        raise IndexError("only Data_Vec<FTYPE_t, dim> of 'dim' 2 or 3 supported")
    for j in range(dim):
        for i in range(size):
            Data_Vec[j][i] = data[j][i]
    return Data_Vec

def non_lin_pow_spec(a, k, cosmo):
    """ return ndarray of nonlinear power spectrum """
    # special wrapper -- use emulator power spectrum regardless of the one in cosmo
    return Non_Linear_Cosmo.non_lin_pow_spec(a, k, cosmo)

def lin_pow_spec(a, k, cosmo):
    """ return ndarray of linear power spectrum """
    k = np.array(k)
    if k.shape:
        return np.array([fs.lin_pow_spec(a, k_, cosmo)  for k_ in k])
    else:
        return fs.lin_pow_spec(a, np.asscalar(k), cosmo)

def chi_bulk_a(a, chi_opt, MPL=1, CHI_A_UNITS=True):
    """ return bulk value of chameleon field at background level """
    if CHI_A_UNITS: return 1
    chi_0 = 2*chi_opt["beta"]*MPL*chi_opt["phi"]
    n = chi_opt["n"]
    return chi_0*pow(a, 3/(1-n))

def chi_bulk_a_n(a, chi_opt, MPL=1, CHI_A_UNITS=True):
    """ return bulk value of chameleon field at background level divided by (1-n), i.e. common factor """
    n = chi_opt["n"]
    return chi_bulk_a(a, chi_opt, MPL=MPL, CHI_A_UNITS=CHI_A_UNITS)/(1-n)

def chi_psi_a(a, chi_opt):
    """ return value of screening potential at given time """
    phi = chi_opt["phi"]
    n = chi_opt["n"]
    phi *= pow(a, (5.-2*n)/(1.-n))
    return phi

def phi_G_prefactor(cosmo, c_kms=299792.458):
    mu_ = 3./2*cosmo.Omega_m * pow(cosmo.H0 * cosmo.h / c_kms, 2)
    return 1/mu_

def phi_G_k(a, k, cosmo, c_kms=299792.458):
    Pk = lin_pow_spec(a, k, cosmo)
    Pk_til = Pk*pow(k/(2*np.pi), 3)
    drho = np.sqrt(Pk_til)
    mu = phi_G_prefactor(cosmo, c_kms=c_kms)
    phi = drho/(mu*a*k*k)
    return phi

def chi_psi_k_a_single(a, cosmo, chi_opt, k_min=1e-5, k_max=1e3, rel_tol=1e-1):
    """ return scale at which hravitational potential is equal to screening potential """
    psi_scr_a = chi_psi_a(a, chi_opt)
    f = lambda k : np.abs(psi_scr_a - phi_G_k(a, k, cosmo))
    k_scr = minimize_scalar(f, bracket=(k_min, k_max), bounds=(k_min, np.inf), tol=1e-12).x

    if f(k_scr)/psi_scr_a < rel_tol:
        return k_scr
    else:
        return 0

def chi_psi_k_a(a, cosmo, chi_opt, k_min=1e-5, k_max=1e3, rel_tol=1e-1):
    a = np.array(a)
    fce = lambda a_ : chi_psi_k_a_single(a_, cosmo, chi_opt, k_min=k_min, k_max=k_max, rel_tol=rel_tol)
    if a.shape:
        return np.array([fce(a_)  for a_ in a])
    else:
        return fce(np.asscalar(a))

def chi_mass_sq(a, cosmo, chi_opt, MPL=1, c_kms=299792.458):
    """ return mass squared of chameleon field sitting at chi_bulk(a, 0) """
    prefactor = (3*MPL*chi_opt["beta"]*cosmo.Omega_m *pow(cosmo.H0 # beta*rho_m,0 / Mpl
               * cosmo.h / c_kms # units factor for 'c = 1' and [L] = Mpc / h
               ,2))
    # evolve rho_m,0 -> rho_m
    a = np.array(a)
    prefactor /= pow(a, 3)
    return prefactor/chi_bulk_a_n(a, chi_opt, MPL=MPL, CHI_A_UNITS=False)

def chi_compton_wavelength(a, cosmo, chi_opt, MPL=1, c_kms=299792.458):
    m_sq = chi_mass_sq(a, cosmo, chi_opt, MPL=MPL, c_kms=c_kms)
    return 1/np.sqrt(m_sq)

def chi_lin_pow_spec(a, k, cosmo, chi_opt, MPL=1, c_kms=299792.458):
    """ return ndarray of linear power spectrum for chameleon in units of chi_prefactor """
    mass_sq = chi_mass_sq(a, cosmo, chi_opt, MPL=MPL, c_kms=c_kms)
    k = np.array(k)
    chi_mod = pow(mass_sq/(mass_sq+k*k), 2)

    if k.shape:
        return chi_mod*np.array([fs.lin_pow_spec(a, k_, cosmo)  for k_ in k])
    else:
        return chi_mod*fs.lin_pow_spec(a, np.asscalar(k), cosmo)

def chi_trans_to_supp(a, k, Pk, cosmo, chi_opt, CHI_A_UNITS=True):
    """ transform input chameleon power spectrum to suppression according to linear prediction """
    Pk_lin = chi_lin_pow_spec(a, k, cosmo, chi_opt)
    return Pk/ (Pk_lin * pow(chi_bulk_a_n(a, chi_opt), 2)) # chi_bulk for normalization of Pk_lin

def chi_trans_to_init(data_list, zeropoint=0):
    """ transform supp (ref: lin) to supp (ref: init) """
    reversed_data_list = data_list[::-1]
    for data in reversed_data_list:
        data[1] += zeropoint - reversed_data_list[-1][1]


def hybrid_pow_spec(a, k, A, cosmo):
    """ return 'hybrid' power spectrum: (1-A)*P_lin(k, a) + A*P_nl """
    return (1-A)*lin_pow_spec(a, k, cosmo) + A*non_lin_pow_spec(a, k, cosmo)

def gen_func(sim, fc_par, fce_lin, fc_nl, Pk=None, z=None, non_lin=False):
    data = fs.Data_Vec_2()
    if Pk: # compute function from given continuous power spectrum
        try:
            fc_par(sim, Pk, data)
        # GSL integration error
        except RuntimeError:
            return None
    elif z is not None: # compute (non-)linear function
        a = 1./(1.+z) if z != 'init' else 1.0
        if non_lin:
            fc_nl(sim, a, data)
        else:
            fce_lin(sim, a, data)
    else:
        raise KeyError("Function 'gen_func' called without arguments.")
    return get_ndarray(data)

def corr_func(sim, Pk=None, z=None, non_lin=False):
    """ return correlation function
    if given Pk -- C++ class Extrap_Pk or Extrap_Pk_Nl -- computes its corr. func.
    if given redshift, computes linear or non-linear (emulator) correlation function  """
    fc_par = fs.gen_corr_func_binned_gsl_qawf
    fce_lin = fs.gen_corr_func_binned_gsl_qawf_lin
    fc_nl = Non_Linear_Cosmo.non_lin_corr_func
    return gen_func(sim, fc_par, fce_lin, fc_nl, Pk=Pk, z=z, non_lin=non_lin)

def sigma_R(sim, Pk=None, z=None, non_lin=False):
    """ return amplitude of density fluctuations
    if given Pk -- C++ class Extrap_Pk or Extrap_Pk_Nl -- computes its sigma_R.
    if given redshift, computes linear or non-linear (emulator) amplitude of density fluctuations  """
    fc_par = fs.gen_sigma_binned_gsl_qawf
    fce_lin = fs.gen_sigma_func_binned_gsl_qawf_lin
    fc_nl = Non_Linear_Cosmo.non_lin_sigma_func
    return gen_func(sim, fc_par, fce_lin, fc_nl, Pk=Pk, z=z, non_lin=non_lin)

def get_hybrid_pow_spec_amp(sim, data, k_nyquist_par, a=None, fit_lin=False):
    """ fit data [k, Pk, std] to hybrid power spectrum (1-A)*P_lin(k) + A*P_nl(k)
    return dictionary with C++ class Extrap_Pk_Nl, fit values and covariance.
    If 'fit_lin' is True, fit linear power spectrum and use C++ class Extrap_Pk instead.
     """
    # extract data
    kk, Pk = np.array(data[0]), np.array(data[1])
    
    # get proper slice of data -- last decade before half of particle nyquist
    idx = (np.abs(kk-0.5*k_nyquist_par)).argmin()
    idx = slice(idx - sim.out_opt.bins_per_decade, idx)

    # do we have errors?
    sigma = data[2][idx] if len(data) > 2 else None

    # get data vector
    data_vec = get_Data_vec(data)

    # are fitting only linear power spectrum?
    if fit_lin:
        Pk_par = fs.Extrap_Pk_2 if sigma is None else fs.Extrap_Pk_3
        Pk_par = Pk_par(data_vec, sim)
        popt = pcov = perr = pcor = None

        # check proper slope of power spectrum : n_s < -1.5
        if Pk_par.n_s > -1.5:
            P_0 = Pk_par(Pk_par.k_max)
            Pk_par.n_s = -1
            Pk_par.A_up *= P_0 / Pk_par(Pk_par.k_max)

    # fit hybrid power spectrum
    else:
        # define functions which will be used in fitting (whether 'a' is free parameter or not)
        if a is None:
            pk_hyb_func = lambda k, A, a: hybrid_pow_spec(a, k, A, sim.cosmo)
            p0 = (0.05, 0.5)
        else:
            pk_hyb_func = lambda k, A: hybrid_pow_spec(a, k, A, sim.cosmo)
            p0 = 0.05

        # fit data, a = <0, 1>, A = <0, 1>        
        bounds = (0, 1)
        
        popt, pcov = curve_fit(pk_hyb_func, kk[idx], Pk[idx], p0=p0, bounds=bounds, sigma=sigma)
        perr, pcor = get_perr_pcor(pcov)

        # get hybrid Extrap
        Pk_par = fs.Extrap_Pk_Nl_2 if sigma is None else fs.Extrap_Pk_Nl_3
        A = popt[0]
        a = popt[1] if a is None else a
        Pk_par = Pk_par(data_vec, sim, A, a)

    # return all info in dict
    return {"Pk_par" : Pk_par, "popt" : popt, "pcov" : pcov, 'perr' : perr, 'pcor' : pcor}

def get_perr_pcor(pcov):
    """convert covariance matrix into standard errors and correlation matrix"""
    perr = np.sqrt(np.diag(pcov))
    inv = np.diag([1/x if x else 0 for x in perr]) # cut out values with 0 variance
    pcor = np.dot(np.dot(inv, pcov), inv)
    return perr, pcor

def gaussian(x, A, mu, sigma):
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def get_bao_peak(corr, r_low=80, r_high=120):
    # get data between r_low and r_high
    r, xi = corr
    idx = np.where(r > r_low)
    x = r[idx]
    y = xi[idx]
    idx = np.where(x < r_high)
    x = x[idx]
    y = x*x*y[idx]
    
    # fit gaussian to the data
    p0 = (y[0], 100, 10)
    popt, pcov = curve_fit(gaussian, x, y, p0=p0)
    perr, pcor = get_perr_pcor(pcov)
    return {"popt" : popt, "pcov" : pcov, 'perr' : perr, 'pcor' : pcor}

def get_truncation(k, cosmo_par):
    k2_G = cosmo_par["smoothing_k"]
    return np.exp(-k*k/k2_G)

def growth_factor(a, cosmo):
    """ return growth factor D at scale factor a, accepts ndarray """
    a = np.array(a)
    if a.shape:
        return np.array([fs.growth_factor(a_, cosmo) for a_ in a])
    else:
        return fs.growth_factor(np.asscalar(a), cosmo)

def growth_rate(a, cosmo):
    """ return growth rate 'f= d ln D/ d ln a' at scale factor a, accepts ndarray """
    a = np.array(a)
    if a.shape:
        return np.array([fs.growth_rate(a_, cosmo) for a_ in a])
    else:
        return fs.growth_rate(np.asscalar(a), cosmo)

def growth_change(a, cosmo):
    """ return growth change 'd D/d a' at scale factor a, accepts ndarray """
    a = np.array(a)
    if a.shape:
        return np.array([fs.growth_change(a_, cosmo) for a_ in a])
    else:
        return fs.growth_change(np.asscalar(a), cosmo)

def get_a_from_growth(D, cosmo):
    """ get scale factor at which growth factor is D """
    D = np.array(D)
    if D.shape:
        a_eff = []
        for D_ in D:
            f = lambda a : D_ - growth_factor(a, cosmo)
            a_eff.append(brentq(f, 0, 2))
        return np.array(a_eff)
    else:
        f = lambda a : D - growth_factor(a, cosmo)
        return brentq(f, 0, 2)
