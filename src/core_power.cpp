
#include "stdafx.h"
#include "core.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>

using namespace std;
const double PI = acos(-1.);

static const double Omega_0 = 1;
static const double h = 0.67;

double transfer_function_2(double k, void* parameters)
{
	if (k == 0) return 1.;
	double q = k / (Omega_0*h);
	double T_k =	log(1+2.34*q)/(2.34*q)*
					pow(1 + 3.89*q + pow(16.1*q, 2.) + pow(5.4*q, 3.) + pow(6.71*q, 4.)
					, -1./4.);
	return pow(T_k, 2.);
}

double power_spectrum_T(double k, double* parameters)
{
	if (k == 0) return 0;
	double A = parameters[0];
	double ns = parameters[1];
	return A*pow(k, ns)*transfer_function_2(k, parameters);
}

double power_spectrum_scale_free(double k, double* parameters)
{
	if (k == 0) return 0;
	double A = parameters[0];
	double ns = parameters[1];
	return A*pow(k, ns);
}

double flat_power_spectrum(double k, double* parameters)
{
	double A = parameters[0];
	return A;
}

double single_power_spectrum_T(double k, double* parameters)
{
	if ((k > 0.01) and (k < 0.04)) return power_spectrum_T(k, parameters);
	else return 0.;
}

double power_spectrum(double k, double* parameters)
{
	e_power_spec pwr_type = static_cast<e_power_spec>(parameters[2]);
	switch (pwr_type)
	{
		case power_law_T: return power_spectrum_T(k, parameters);
		case power_law: return power_spectrum_scale_free(k, parameters);
		case flat: return flat_power_spectrum(k, parameters);
		case single: return single_power_spectrum_T(k, parameters);
		default: return power_spectrum_T(k, parameters);
	}
}


double power_spectrum_s8(double k, void* parameters)
{
	return	k*k/(2.*PI*PI)* // spherical factor
			power_spectrum(k, (double*)parameters)* // P(k)
			pow(3.*gsl_sf_bessel_j1(k*8)/(k*8), 2.); // window function R = 8 Mpc/h
}
	
void norm_pwr(Pow_Spec_Param* pwr_par)
{
	/* Normalize the power spectrum */
	printf("Computing normalization of the given power spectrum...\n");
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	double result, error;
	
	double parameters[3] = {pwr_par->A, pwr_par->ns, static_cast<double>(pwr_par->pwr_type)};
	gsl_function F;
	F.function = &power_spectrum_s8;
	F.params = parameters;
	gsl_integration_qagiu (&F, 0, 0, 1e-7, 1000, w, &result, &error); 
	
	gsl_integration_workspace_free (w);
	pwr_par->A = pwr_par->s8*pwr_par->s8/result;
}

double lin_pow_spec(Pow_Spec_Param pwr_par, double k)
{
	double parameters[3] = {pwr_par.A, pwr_par.ns, static_cast<double>(pwr_par.pwr_type)};
	return power_spectrum(k, parameters);
}