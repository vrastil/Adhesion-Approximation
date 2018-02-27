/**
 * @file:	core_power.h
 * @brief:	functions handling operations with power spectra
 */
#pragma once
#include "stdafx.h"
#include "params.hpp"

void norm_pwr(Cosmo_Param& cosmo);
FTYPE norm_growth_factor(const Cosmo_Param& cosmo);
FTYPE growth_factor(FTYPE a, const Cosmo_Param& cosmo);
FTYPE growth_rate(FTYPE a, const Cosmo_Param& cosmo);
FTYPE growth_change(FTYPE a, const Cosmo_Param& cosmo);
FTYPE Omega_lambda(FTYPE a, const Cosmo_Param& cosmo);
FTYPE lin_pow_spec(FTYPE a, FTYPE k, const Cosmo_Param& cosmo);
FTYPE non_lin_pow_spec(FTYPE a, FTYPE k, const Cosmo_Param& cosmo);

/**
 * @class:	ODE_Solver
 * @brief:	Explicit embedded Runge-Kutta Prince-Dormand (8, 9) method.
 */

class ODE_Solver
{
public:
    ODE_Solver(int (* function) (double t, const double y[], double dydt[], void * params), size_t dim, void* params,
               const double hstart = 1e-6, const double epsabs = 1e-6, const double epsrel = 0);
    ~ODE_Solver();
    void update(double &t, const double t1, double y[]);
    int status;
    gsl_odeiv2_system sys;
    gsl_odeiv2_driver* d;
};

/**
 * @class:	Interp_obj
 * @brief:	linear interpolation of data [x, y]
 */

class Interp_obj
{// Steffen interpolation of data [x, y]
public:
    Interp_obj(): is_init(false) {}
    ~Interp_obj();
    double operator()(double x) const;
    template <typename T, unsigned N>
    void init(const Data_Vec<T, N>& data);
    double x_min, x_max;

private:
    bool is_init;
    gsl_spline* spline;
    gsl_interp_accel* acc;
};

/**
 * @class:	Extrap_Pk
 * @brief:	linear interpolation of data [k, P(k)] within 'useful' range
            fit to primordial P_i(k) below the 'useful' range
            fit to Padé approximant R [0/3] above the 'useful' range
 */

template <typename T, unsigned N>
class Extrap_Pk : public Interp_obj
{ /*
    Steffen interpolation of data [k, P(k)] within range k_min, k_max
    fit to primordial P_i(k) below this range, fit A*k^ns above
*/
public:
    Extrap_Pk(const Data_Vec<T, N>& data, const Sim_Param& sim);
    Extrap_Pk(const Data_Vec<T, N>& data, const Sim_Param& sim, const unsigned m_l, const unsigned n_u);
    Extrap_Pk(const Data_Vec<T, N>& data, const Sim_Param& sim, const unsigned m_l, const unsigned n_l,
              const unsigned m_u, const unsigned n_u);
    double operator()(double k) const;

    void fit_lin(const Data_Vec<T, N>& data, const unsigned m, const unsigned n, double& A);
    void fit_power_law(const Data_Vec<T, N>& data, const unsigned m, const unsigned n, double& A, double& n_s);

    double A_low; // amplitude of linear power in lower range
    const Cosmo_Param& cosmo;
    double A_up, n_s; // scale-free power spectrum in upper range
    T k_min, k_max; // interpolation range
};

/**
 * @class:	Extrap_Pk_Nl
 * @brief:	creates Extrapolate object (linear power spectrum) from data and store non-linear parameters
            call 'operator()(k)' based on k_split (upper range of the linear)
 */

template <typename T, unsigned N>
class Extrap_Pk_Nl : public Extrap_Pk<T, N>
{
public:
    Extrap_Pk_Nl(const Data_Vec<T, N>& data, const Sim_Param &sim, T A_nl, T a_eff);
    const T A_nl, a_eff, k_split;
    double operator()(double k) const;
};

template<class P> // everything callable P_k(k)
void gen_corr_func_binned_gsl_qawf(const Sim_Param &sim, const P& P_k, Data_Vec<FTYPE, 2>& corr_func_binned);
void gen_corr_func_binned_gsl_qawf_lin(const Sim_Param &sim, FTYPE a, Data_Vec<FTYPE, 2>& corr_func_binned);
void gen_corr_func_binned_gsl_qawf_nl(const Sim_Param &sim, FTYPE a, Data_Vec<FTYPE, 2>& corr_func_binned);