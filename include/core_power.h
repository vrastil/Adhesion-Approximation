/**
 * @file:	core_power.h
 * @brief:	functions handling operations with power spectra
 */
#pragma once
 
#include "stdafx.h"
#include "core.h"
#include "emu.h"


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

class Extrap_Pk_Nl
{/*
    creates Extrapolate object (linear power spectrum) from data and store non-linear parameters
    call 'operator()(k)' based on k_split (upper range of the linear)
*/
public:
    Extrap_Pk_Nl(const Data_Vec<double, 3>& data, const Sim_Param &sim, double A_nl, double z_eff):
        Pk_lin(data, sim), Pk_nl(emu::init_emu(sim, z_eff > 2.02 ? 2.02 : z_eff < 0 ? 0 : z_eff), sim, 0, 10, 341, 351),
        A_nl(A_nl), k_split(Pk_lin.k_max), D(growth_factor(1./(1.+z_eff), sim.cosmo)) {}
    const Extrap_Pk Pk_lin, Pk_nl;
    const double A_nl, k_split, D;
    double operator()(double k) const { 
        return k < k_split ? Pk_lin(k) : A_nl*Pk_nl(k) + D*D*(1-A_nl)*lin_pow_spec(k, Pk_nl.cosmo);
        }
};

template<class P> // everything callable P_k(k)
void gen_corr_func_binned_gsl_qawf(const Sim_Param &sim, const P& P_k, Data_Vec<double, 2>* corr_func_binned);
void gen_corr_func_binned_gsl_qawf_lin(const Sim_Param &sim, double a, Data_Vec<double, 2>* corr_func_binned);