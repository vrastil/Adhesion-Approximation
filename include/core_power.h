/**
 * @file:	core_power.h
 * @brief:	functions handling operations with power spectra
 */
#pragma once
 
#include "stdafx.h"
#include "core.h"

void norm_pwr(Cosmo_Param& cosmo);
double norm_growth_factor(const Cosmo_Param& cosmo);
double growth_factor(double a, const Cosmo_Param& cosmo);
double growth_rate(double a, const Cosmo_Param& cosmo);
double growth_change(double a, const Cosmo_Param& cosmo);
double Omega_lambda(double a, const Cosmo_Param& cosmo);
double lin_pow_spec(double a, double k, const Cosmo_Param& cosmo);
double non_lin_pow_spec(double a, double k, const Cosmo_Param& cosmo);

/**
 * @class:	Interp_obj
 * @brief:	linear interpolation of data [x, y]
 */

class Interp_obj
{// Steffen interpolation of data [x, y]
public:
    Interp_obj(): is_init(false) {}
    template<unsigned N>
    Interp_obj(const Data_Vec<double, N>& data);
    ~Interp_obj();
    double operator()(double x) const;
    template<unsigned N>
    void init(const Data_Vec<double, N>& data);

    double x_min, x_max;

private:
    bool is_init;
    gsl_spline* spline;
    gsl_interp_accel* acc;
};

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
 * @class:	Extrap_Pk
 * @brief:	linear interpolation of data [k, P(k)] within 'useful' range
            fit to primordial P_i(k) below the 'useful' range
            fit to Padé approximant R [0/3] above the 'useful' range
 */

class Extrap_Pk : public Interp_obj
{ /*
    Steffen interpolation of data [k, P(k)] within range k_min, k_max
    fit to primordial P_i(k) below this range, fit A*k^ns above
*/
public:
    template<unsigned N>
    Extrap_Pk(const Data_Vec<double, N>& data, const Sim_Param& sim);
    template<unsigned N>
    Extrap_Pk(const Data_Vec<double, N>& data, const Sim_Param& sim, const unsigned m_l, const unsigned n_u);
    template<unsigned N>
    Extrap_Pk(const Data_Vec<double, N>& data, const Sim_Param& sim, const unsigned m_l, const unsigned n_l,
              const unsigned m_u, const unsigned n_u);
    double operator()(double k) const;

    template<unsigned N>
    void fit_lin(const Data_Vec<double, N>& data, const unsigned m, const unsigned n, double& A);
    template<unsigned N>
    void fit_power_law(const Data_Vec<double, N>& data, const unsigned m, const unsigned n, double& A, double& n_s);

    double A_low; // amplitude of linear power in lower range
    const Cosmo_Param& cosmo;
    double A_up, n_s; // scale-free power spectrum in upper range
    double k_min, k_max; // interpolation range
};

/**
 * @class:	Extrap_Pk_Nl
 * @brief:	creates Extrapolate object (linear power spectrum) from data and store non-linear parameters
            call 'operator()(k)' based on k_split (upper range of the linear)
 */

class Extrap_Pk_Nl
{
public:
    Extrap_Pk_Nl(const Data_Vec<double, 3>& data, const Sim_Param &sim, double A_nl, double a_eff);
    const Extrap_Pk Pk_lin;
    const double A_nl, a_eff, k_split;
    double operator()(double k) const;
};

template<class P> // everything callable P_k(k)
void gen_corr_func_binned_gsl_qawf(const Sim_Param &sim, const P& P_k, Data_Vec<double, 2>* corr_func_binned);
void gen_corr_func_binned_gsl_qawf_lin(const Sim_Param &sim, double a, Data_Vec<double, 2>* corr_func_binned);