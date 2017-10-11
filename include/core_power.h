/**
 * @file:	core_power.h
 * @brief:	functions handling operations with power spectra
 */
 
#include "stdafx.h"
#include "core.h"

void norm_pwr(Pow_Spec_Param* pwr_par);
double lin_pow_spec(const Pow_Spec_Param* pwr_par, double k);
void gen_corr_func_binned_gsl(const Sim_Param &sim, const double x_min, const double x_max, const Data_x_y<double>& pwr_spec_binned, Data_x_y<double>* corr_func_binned);