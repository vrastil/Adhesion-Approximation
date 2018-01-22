
"""
'power.py' module serves mainly for interacting with C++ library fastsim.py

 - translate all power spectra, growth functions, correlations functions, etc.
   into C++ functions for speed
 - handles numpy arrays
 - cosmo == C++ class Cosmo_Param, accessible through SimInfo.sim.cosmo
 - FTYPE=[float, double, long double]
"""

import numpy as np
from scipy.optimize import brentq
from . import fastsim as fs

def get_a_from_A(cosmo, A):
    """ return scale factor a at which amplitude of linear power spectrum is A (normalize as A=1 at a=1) """
    # 'f = 0' <=> A = D^2 (linear power grows as D^2)
    f = lambda a : A - fs.growth_factor(a, cosmo)**2
    return brentq(f, 0, 1)

def get_ndarray(Data_Vec):
    """ copy C++ class Data_Vec<FTYPE, N> into numpy array """
    dim = Data_Vec.dim()
    data = [[x for x in Data_Vec[i]] for i in xrange(dim)]
    return np.array(data)

def get_Data_vec(data):
    """ copy 2D data 'dim x size' into C++ class Data_Vec<FTYPE, dim> """
    dim = len(data)
    size = len(data[0])
    if dim == 2:
        Data_Vec = fs.Data_Vec_2(size)
    elif dim == 3:
        Data_Vec = fs.Data_Vec_3(size)
    else:
        raise IndexError("only Data_Vec<FTYPE, dim> of 'dim' 2 or 3 supported")
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

def hybrid_pow_spec(a, k, A, cosmo):
    """ return 'hybrid' power spectrum: (1-A)*P_lin(k, a) + A*P_nl """
    return (1-A)*lin_pow_spec(a, k, cosmo) + A*non_lin_pow_spec(a, k, cosmo)

def corr_func(sim, Pk=None, z=None, non_lin=False):
    """ return correlation function
    if given Pk -- C++ class Extrap_Pk or Extrap_Pk_Nl -- computes its corr. func.
    if given redshift, computes linear or non-linear (emulator) correlation function  """
    corr = fs.Data_Vec_2()
    if Pk: # compute correlation function from given continuous power spectrum
        fs.gen_corr_func_binned_gsl_qawf(sim, Pk, corr)
    elif z is not None: # compute (non-)linear correlation function
        a = 1./(1.+z) if z != 'init' else 1.0
        if non_lin:
            fs.gen_corr_func_binned_gsl_qawf_nl(sim, a, corr)
        else:
            fs.gen_corr_func_binned_gsl_qawf_lin(sim, a, corr)
    else:
        raise KeyError("Function 'corr_func' called without arguments.")
    return get_ndarray(corr)

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
