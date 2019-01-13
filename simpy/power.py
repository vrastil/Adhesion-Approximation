
"""
'power.py' module serves mainly for interacting with C++ library fastsim.py

 - translate all power spectra, growth functions, correlations functions, etc.
   into C++ functions for speed
 - handles numpy arrays
 - cosmo == C++ class Cosmo_Param, accessible through SimInfo.sim.cosmo
 - FTYPE_t=[float, double, long double]
"""

import numpy as np
from scipy.optimize import brentq
from scipy.signal import argrelextrema
from . import fastsim as fs

def get_a_from_A(cosmo, A):
    """ return scale factor a at which amplitude of linear power spectrum is A (normalize as A=1 at a=1) """
    # 'f = 0' <=> A = D^2 (linear power grows as D^2)
    f = lambda a : A - fs.growth_factor(a, cosmo)**2
    return brentq(f, 0, 1)

def get_ndarray(Data_Vec):
    """ copy C++ class Data_Vec<FTYPE_t, N> into numpy array """
    dim = Data_Vec.dim()
    data = [[x for x in Data_Vec[i]] for i in xrange(dim)]
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
    for j in xrange(dim):
        for i in xrange(size):
            Data_Vec[j][i] = data[j][i]
    return Data_Vec

def non_lin_pow_spec(a, k, cosmo):
    """ return ndarray of nonlinear power spectrum """
    k = np.array(k)
    if k.shape:
        return np.array([fs.non_lin_pow_spec(a, k_, cosmo)  for k_ in k])
    else:
        return fs.non_lin_pow_spec(a, np.asscalar(k), cosmo)

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

def chi_mass_sq(a, cosmo, chi_opt, MPL=1, c_kms=299792.458):
    """ return mass squared of chameleon field sitting at chi_bulk(a, 0) """
    prefactor = (3*MPL*chi_opt["beta"]*cosmo.Omega_m *pow(cosmo.H0 # beta*rho_m,0 / Mpl
               * cosmo.h / c_kms # units factor for 'c = 1' and [L] = Mpc / h
               ,2))
    # evolve rho_m,0 -> rho_m
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
        fc_par(sim, Pk, data)
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
    fc_nl = fs.gen_corr_func_binned_gsl_qawf_nl
    return gen_func(sim, fc_par, fce_lin, fc_nl, Pk=Pk, z=z, non_lin=non_lin)

def sigma_R(sim, Pk=None, z=None, non_lin=False):
    """ return amplitude of density fluctuations
    if given Pk -- C++ class Extrap_Pk or Extrap_Pk_Nl -- computes its sigma_R.
    if given redshift, computes linear or non-linear (emulator) amplitude of density fluctuations  """
    fc_par = fs.gen_sigma_binned_gsl_qawf
    fce_lin = fs.gen_sigma_func_binned_gsl_qawf_lin
    fc_nl = fs.gen_sigma_func_binned_gsl_qawf_nl
    return gen_func(sim, fc_par, fce_lin, fc_nl, Pk=Pk, z=z, non_lin=non_lin)

def get_bao_peak(corr, cutof=80):
    r, xi = corr
    
    # get id of local maxima
    idx = argrelextrema(xi, np.greater)
    
    # get first maximum after 25 Mpc/h (default cutof)
    peak_id = [x for x in idx[0] if r[x] > cutof][0]
    
    # return peak location and amplitude
    return r[peak_id], xi[peak_id]

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
            a_eff.append(brentq(f, 0, 1))
        return np.array(a_eff)
    else:
        f = lambda a : D - growth_factor(a, cosmo)
        return brentq(f, 0, 1)
