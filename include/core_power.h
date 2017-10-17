/**
 * @file:	core_power.h
 * @brief:	functions handling operations with power spectra
 */
#pragma once
 
#include "stdafx.h"
#include "core.h"


void norm_pwr(Pow_Spec_Param* pwr_par);
double lin_pow_spec(double k, const Pow_Spec_Param& parameters);
double lin_pow_spec(double k, void* parameters);
double  get_max_Pk(Sim_Param* sim);

/**
 * @class:	Interp_obj
 * @brief:	linear interpolation of data [x, y]
 */

class Interp_obj
{// linear interpolation of data [x, y]
public:
    Interp_obj(const Data_x_y<double>& data);
    ~Interp_obj();
    double eval(double x) const;

private:
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
    Extrap_Pk(const Data_x_y<double>& data, const Sim_Param& sim);
    double eval(double k) const;

private:
    double n_s, A; // lower range, priomordial
    std::vector<double> a_m, b_n; // upper range, Pade approximant
    double k_min, k_max; // interpolation range
};

void gen_corr_func_binned_gsl_qagi(const Sim_Param &sim, const Extrap_Pk& P_k, Data_x_y<double>* corr_func_binned);
void gen_corr_func_binned_gsl_qawo(const Sim_Param &sim, const Extrap_Pk& P_k, Data_x_y<double>* corr_func_binned);
void gen_corr_func_binned_gsl_qawf(const Sim_Param &sim, const Extrap_Pk& P_k, Data_x_y<double>* corr_func_binned);