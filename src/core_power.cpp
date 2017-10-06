
#include "stdafx.h"
#include "core.h"
#include "core_mesh.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_spline.h>

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
	const double ns = parameters[1];
	return pow(k, ns)*transfer_function_2(k, parameters);
}

double power_spectrum_scale_free(double k, double* parameters)
{
	if (k == 0) return 0;
	const double ns = parameters[1];
	return pow(k, ns);
}

double flat_power_spectrum(double k, double* parameters)
{
	return 1;
}

double single_power_spectrum_T(double k, double* parameters)
{
	if ((k > 0.01) and (k < 0.04)) return power_spectrum_T(k, parameters);
	else return 0.;
}

double power_spectrum(double k, double* parameters)
{
    const double A = parameters[0];
    const double supp = parameters[3] ? exp(-k*k/parameters[3]) : 1;
	e_power_spec pwr_type = static_cast<e_power_spec>(parameters[2]);
	switch (pwr_type)
	{
		case power_law_T: return supp*A*power_spectrum_T(k, parameters);
		case power_law: return supp*A*power_spectrum_scale_free(k, parameters);
		case flat: return supp*A*flat_power_spectrum(k, parameters);
		case single: return supp*A*single_power_spectrum_T(k, parameters);
		default: return supp*A*power_spectrum_T(k, parameters);
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
	
	double parameters[4] = {1., pwr_par->ns, static_cast<double>(pwr_par->pwr_type), pwr_par->k2_G};
	gsl_function F;
	F.function = &power_spectrum_s8;
	F.params = parameters;
	gsl_integration_qagiu (&F, 0, 0, 1e-7, 1000, w, &result, &error); 
	
	gsl_integration_workspace_free (w);
	pwr_par->A = pwr_par->s8*pwr_par->s8/result;
}

double lin_pow_spec(Pow_Spec_Param pwr_par, double k)
{
	double parameters[4] = {pwr_par.A, pwr_par.ns, static_cast<double>(pwr_par.pwr_type), pwr_par.k2_G};
	return power_spectrum(k, parameters);
}

class Interp_obj
{
public:
    // CONSTRUCTOR
    Interp_obj(const Data_x_y<double>& data)
    {
        acc = gsl_interp_accel_alloc ();
        spline = gsl_spline_alloc (gsl_interp_linear, data.size());
        gsl_spline_init (spline, data.x.data(), data.y.data(), data.size());        
    }
    // DESTRUCTOR
    ~Interp_obj()
    {
        gsl_spline_free(spline);
        gsl_interp_accel_free (acc);
    }
    // METHODS
    double eval(double x) const{ return gsl_spline_eval(spline, x, acc); }

private:
    // VARIABLES
    gsl_spline* spline;
    gsl_interp_accel* acc;
};

template <class spec_params>
// class [spec_params] has to have parameter r;
class Integr_obj_qawo
{
public:
    // CONSTRUCTOR
    Integr_obj_qawo(double(*f) (double, void*), spec_params* params,  double a, double b, size_t limit, size_t n):
    a(a), L(b-a), limit(limit), params(params)
    {
        w = gsl_integration_workspace_alloc (limit);
        t = gsl_integration_qawo_table_alloc(1, 1, GSL_INTEG_SINE, n);
        F.function = f;
    }
    // DESTRUCTOR
    ~Integr_obj_qawo()
    {
        gsl_integration_workspace_free (w);
        gsl_integration_qawo_table_free(t);
    }
    // METHODS
    double operator()(double r)
    {
        gsl_integration_qawo_table_set(t, r, L, GSL_INTEG_SINE);
        params->r = r;
        F.params = params;
        gsl_integration_qawo(&F, a, 1e-3, 1e-3, limit, w, t, &result, &error);
        return result;
    }

    // VARIABLES
    double result, error;
private:
    const double a, L;
    const size_t limit;
    spec_params* params;
    gsl_function F;
    gsl_integration_qawo_table* t;
    gsl_integration_workspace* w;
};

struct xi_integrand_param
{
    double r;
    Interp_obj* P_k;
};

double xi_integrand(double k, void* params){
    xi_integrand_param* my_par = (xi_integrand_param*) params;
    const double r = my_par->r;
    Interp_obj* P_k = my_par->P_k;
    return 1/(2*PI*PI)*k/r*P_k->eval(k);
};

void gen_corr_func_binned_gsl(const double x_min, const double x_max, Data_x_y<double>* corr_func_binned)
{
    Data_x_y<double> pwr_spec_binned = *corr_func_binned; // prevent overwriting

    printf("Allocationg space for interpolation function P(k)\n");
    Interp_obj P_k(pwr_spec_binned);

    printf("Allocationg space for integration via [QAWO adaptive integration for oscillatory functions]\n");
    const double k_min = min(pwr_spec_binned.x);
    const double k_max = max(pwr_spec_binned.x);
    xi_integrand_param my_param;
    my_param.P_k = &P_k;
    Integr_obj_qawo<xi_integrand_param> xi_r(&xi_integrand, &my_param, k_min, k_max, 1000, 25);

    printf("Computing correlation function via [QAWO adaptive integration for oscillatory functions]...\n");
    const unsigned int N = corr_func_binned->size();
    const double lin_bin = (x_max - x_min)/N;
	double r;
	for(unsigned i = 0; i < N; i++){
        r = x_min + i*lin_bin;
        corr_func_binned->x[i] = r;
        corr_func_binned->y[i] = xi_r(r);
    }
}