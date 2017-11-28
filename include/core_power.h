/**
 * @file:	core_power.h
 * @brief:	functions handling operations with power spectra
 */
#pragma once
 
#include "stdafx.h"
#include "core.h"


void norm_pwr(Cosmo_Param* pwr_par);
double lin_pow_spec(double k, const Cosmo_Param& pwr_par, double a);
double lin_pow_spec(double k, const Cosmo_Param& parameters);
double growth_factor(double a, const Cosmo_Param& pwr_par);
double growth_rate(double a, const Cosmo_Param& pwr_par);
double growth_change(double a, const Cosmo_Param& pwr_par);
double  get_max_Pk(Sim_Param* sim);

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
    double eval(double x) const;
    template<unsigned N>
    void init(const Data_Vec<double, N>& data);

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

class Extrap_Pk : public Interp_obj
{ /*
    linear interpolation of data [k, P(k)] within 'useful' range
    fit to primordial P_i(k) below the 'useful' range
    fit to Padé approximant R [0/3] above the 'useful' range
*/
public:
    template<unsigned N>
    Extrap_Pk(const Data_Vec<double, N>& data, const Sim_Param& sim);
    template<unsigned N>
    Extrap_Pk(const Data_Vec<double, N>& data, const Sim_Param& sim, const unsigned m_l, const unsigned n_u);
    template<unsigned N>
    Extrap_Pk(const Data_Vec<double, N>& data, const Sim_Param& sim, const unsigned m_l, const unsigned n_l,
              const unsigned m_u, const unsigned n_u);
    double eval(double k) const;

private:
    template<unsigned N>
    void fit_lin(const Data_Vec<double, N>& data, const unsigned m, const unsigned n, double& A);
    template<unsigned N>
    void fit_power_law(const Data_Vec<double, N>& data, const unsigned m, const unsigned n, double& A, double& n_s);

    double A_low; // amplitude of linear power in lower range
    const Cosmo_Param& cosmo;
    double A_up, n_s; // scale-free power spectrum in upper range
    double k_min, k_max; // interpolation range
};

void gen_corr_func_binned_gsl_qawf(const Sim_Param &sim, const Extrap_Pk& P_k, Data_Vec<double, 2>* corr_func_binned);
void gen_corr_func_binned_gsl_qawf_lin(const Sim_Param &sim, double a, Data_Vec<double, 2>* corr_func_binned);