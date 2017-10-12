
#include "stdafx.h"
#include "core.h"
#include "core_mesh.h"
#include "core_power.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_linalg.h>
#include <ccl.h>

using namespace std;

double transfer_function_2(double k, const Pow_Spec_Param* parameters)
{
    if (k == 0) return 1.;

	const double q = k / (parameters->Omega_m()*parameters->h);
	double T_k =	log(1+2.34*q)/(2.34*q)*
					pow(1 + 3.89*q + pow(16.2*q, 2.) + pow(5.47*q, 3.) + pow(6.71*q, 4.)
					, -1./4.);
	return pow(T_k, 2.);
}

double power_spectrum_T(double k, const Pow_Spec_Param* parameters)
{
	if (k == 0) return 0;
	const double ns = parameters->ns;
	return pow(k, ns)*transfer_function_2(k, parameters);
}

double power_spectrum_scale_free(double k, const Pow_Spec_Param* parameters)
{
	if (k == 0) return 0;
	const double ns = parameters->ns;
	return pow(k, ns);
}

double flat_power_spectrum(double k, const Pow_Spec_Param* parameters)
{
	return 1;
}

double single_power_spectrum_T(double k, const Pow_Spec_Param* parameters)
{
	if ((k > 0.01) and (k < 0.04)) return power_spectrum_T(k, parameters);
	else return 0.;
}

double power_spectrum(double k, const Pow_Spec_Param* parameters)
{
    const double A = parameters->A;
    const double supp = parameters->k2_G ? exp(-k*k/parameters->k2_G) : 1;
	switch (parameters->pwr_type)
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
	return	k*k/(2.*PI*PI) // spherical factor
			*power_spectrum(k, static_cast<Pow_Spec_Param*>(parameters)) // P(k)
			*pow(3.*gsl_sf_bessel_j1(k*8)/(k*8), 2.); // window function R = 8 Mpc/h
}

void norm_pwr_gsl(Pow_Spec_Param* pwr_par)
{
	/* Normalize the power spectrum */
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	double result, error;
	
	gsl_function F;
	F.function = &power_spectrum_s8;
	F.params = pwr_par;
	gsl_integration_qagiu (&F, 0, 0, 1e-7, 1000, w, &result, &error); 
	
	gsl_integration_workspace_free (w);
	pwr_par->A = pwr_par->sigma8*pwr_par->sigma8/result;
}

void norm_pwr_ccl(Pow_Spec_Param* pwr_par)
{
    /* Normalize the power spectrum */
    int status = 0;
    ccl_sigma8(pwr_par->cosmo, &status);
}

void norm_pwr(Pow_Spec_Param* pwr_par)
{
    printf("Initializing given power spectrum...\n");
    if (pwr_par->pwr_type < 4) norm_pwr_gsl(pwr_par);
    else norm_pwr_ccl(pwr_par);
}

double lin_pow_spec(const Pow_Spec_Param* pwr_par, double k)
{
    if (pwr_par->pwr_type < 4){
        return power_spectrum(k, pwr_par);
    } else {
        int status = 0;
        return ccl_linear_matter_power(pwr_par->cosmo, k*pwr_par->h, 1, &status)*pow(pwr_par->h, 3);
    }
}

double pade_approx(double x, unsigned m, const double* a, unsigned n, const double* b)
{
    double numerator = 0;
    double denominator = 1;
    for(unsigned i = 0; i <= m; i++){
        numerator += a[i]*pow(x, i);
    }
    for(unsigned i = 0; i <= n; i++){
        denominator += b[i]*pow(x, i+1); // note that b[0] = b1
    }
    return numerator/denominator;
}

double pade_approx(double x, const vector<double>& a, const vector<double>& b)
{
    return pade_approx(x, a.size(), a.data(), b.size(), b.data());
}

int get_nearest(const double val, const vector<double>& vec)
{
    int pos = 0;
    double dv = abs(vec[0] - val);
    double dv_;

    for(unsigned i = 1; i <vec.size(); i++)
    {
        dv_ = abs(vec[i] - val);
        if (dv_ < dv)
        {
            dv = dv_;
            pos = i;
        }
    }
    return pos;
}

/**
 * @class:	Interp_obj
 * @brief:	linear interpolation of data [x, y]
 */

Interp_obj::Interp_obj(const Data_x_y<double>& data)
{
    acc = gsl_interp_accel_alloc ();
    spline = gsl_spline_alloc (gsl_interp_linear, data.size());
    gsl_spline_init (spline, data.x.data(), data.y.data(), data.size());
}

Interp_obj::~Interp_obj()
{
    gsl_spline_free(spline);
    gsl_interp_accel_free (acc);
}

double Interp_obj::eval(double x) const{ return gsl_spline_eval(spline, x, acc); }

/**
 * @class:	Extrap_Pk
 * @brief:	linear interpolation of data [k, P(k)] within 'useful' range
            fit to primordial P_i(k) below the 'useful' range
            fit to Padé approximant R [0/3] above the 'useful' range
 */

Extrap_Pk::Extrap_Pk(const Data_x_y<double>& data, const Sim_Param& sim):
Interp_obj(data), n_s(sim.power.ns)
{
    {   // LOWER RANGE -- fit Ak^ns to data[0:n]
        printf("Fitting amplitude of P(k) in lower range.\n");
        constexpr int n = 10;
        k_min = data.x[n];
        vector<double> k, Pk;

        for(unsigned i = 0; i < n; i++){
            k.push_back(pow(data.x[i], n_s));
            Pk.push_back(data.y[i]);
        }
        double A_sigma2, sumsq;
        gsl_fit_mul(k.data(), 1, Pk.data(), 1, n, &A, &A_sigma2, &sumsq);
        printf("\t[fit A = %.12f, sigma = %.12f, sumsq = %.12f]\n", A, sqrt(A_sigma2), sumsq);
    }
    {   // UPPER RANGE -- solve Ax=b to get Pade approximant
        printf("Computing Padé approximant of P(k) in upper range.\n");
        unsigned order = sim.k_par.pade_order; // Pade approximant R [0/order-1]

        vector<double> R_0m(order);
        k_max = sim.k_par.k_interp.upper;
        const int n = get_nearest(k_max, data.x);
        printf("\t[k_max = %.5e, n = %i, k[n] = %.5e]\n", k_max, n, data.x[n]);

        vector<double> A, b;
        double k, Pk;

        // TEMPORARY!!!
        sim.k_par.k_pade.clear();

        for(unsigned i = 0; i < order; i++)
        {
            k = data.x[n-i*7]; // TEMPORARY!!! these two lines are to be other way around
            sim.k_par.k_pade.push_back(k); // i.e. load k from simulation parameters

            Pk = data.y[n-i*7];
            A.push_back(1.); // a0
            for(unsigned j = 1; j < order; j++)
            {
                A.push_back(-pow(k, j)*Pk); // b_i
            }
            b.push_back(Pk); // P(k)
        }

        gsl_matrix_view A_ = gsl_matrix_view_array (A.data(), order, order);
        gsl_vector_view b_ = gsl_vector_view_array (b.data(), order);
        gsl_vector_view x_ = gsl_vector_view_array (R_0m.data(), order);
        gsl_permutation * p = gsl_permutation_alloc (order);
        int s;
    
        gsl_linalg_LU_decomp (&A_.matrix, p, &s);
        gsl_linalg_LU_solve (&A_.matrix, p, &b_.vector, &x_.vector);

        a_m.push_back(R_0m[0]);
        b_n = vector<double>(R_0m.begin() + 1, R_0m.end());
        gsl_permutation_free (p);
    }
    printf("\t[a0 = %f", a_m[0]);
    for(unsigned j = 0; j < b_n.size(); j++) printf(", b%i = %f", j+1, b_n[j]);
    printf("]\n"); 
}


double Extrap_Pk::eval(double k) const
{
    double val;
    if (k < k_min)
    {
        return A*pow(k, n_s);
    }
    else if (k <= k_max) return Interp_obj::eval(k);
    else
    {
        return pade_approx(k, a_m, b_n);
    }
}

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
    const Extrap_Pk* P_k;
};

double xi_integrand(double k, void* params){
    xi_integrand_param* my_par = (xi_integrand_param*) params;
    const double r = my_par->r;
    const Extrap_Pk* P_k = my_par->P_k;
    return 1/(2*PI*PI)*k/r*P_k->eval(k);
};

void gen_corr_func_binned_gsl(const Sim_Param &sim, const double x_min, const double x_max, const Extrap_Pk& P_k, Data_x_y<double>* corr_func_binned)
{
    const double k_min = 0;
    const double k_max = 1e6;

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