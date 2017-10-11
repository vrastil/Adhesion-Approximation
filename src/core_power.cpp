
#include "stdafx.h"
#include "core.h"
#include "core_mesh.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_poly.h>
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

void polynomialfit(int obs, int degree, double *dx, double *dy, double *store) /* n, p */
{ // SOURCE: http://rosettacode.org/wiki/Polynomial_regression#C
    gsl_multifit_linear_workspace *ws;
    gsl_matrix *cov, *X;
    gsl_vector *y, *c;
    double chisq;

    int i, j;

    X = gsl_matrix_alloc(obs, degree);
    y = gsl_vector_alloc(obs);
    c = gsl_vector_alloc(degree);
    cov = gsl_matrix_alloc(degree, degree);

    for(i=0; i < obs; i++){
        for(j=0; j < degree; j++)
        {
            gsl_matrix_set(X, i, j, pow(dx[i], j));
        }
        gsl_vector_set(y, i, dy[i]);
    }

    ws = gsl_multifit_linear_alloc(obs, degree);
    gsl_multifit_linear(X, y, c, cov, &chisq, ws);

    /* store result ... */
    for(i=0; i < degree; i++)
    {
        store[i] = gsl_vector_get(c, i);
    }

    gsl_multifit_linear_free(ws);
    gsl_matrix_free(X);
    gsl_matrix_free(cov);
    gsl_vector_free(y);
    gsl_vector_free(c);
}

class Interp_obj
{// linear interpolation of data [x, y]
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

class Extrap_obj : public Interp_obj
{ /*
    linear interpolation of data [x, y] within range <x_inter0, x_inter1>
    polynomial extrapolation of degree n_0 / n_1 outside this range
    if n_0 or n_1 is negative, then [p(x, |n|)]^-1 is fitted
    fitting is performed in data range <x_fit00, x_fit01> / <x_fot10, x_fit11>
        using n_fit equidistant values in logspace 
*/
public:
    // CONSTRUCTOR
    Extrap_obj(const Data_x_y<double>& data, double x_inter0, double x_inter1,
               int n_0, double x_fit00, double x_fit01,
               int n_1, double x_fit10, double x_fit11,
               unsigned n_fit, bool const_term_0, bool const_term_1):
    Interp_obj(data), x_inter0(x_inter0), x_inter1(x_inter1),
    n_0(n_0), x_fit00(x_fit00), x_fit01(x_fit01), 
    n_1(n_1), x_fit10(x_fit10), x_fit11(x_fit11),
    coeff_0(abs(n_0)), coeff_1(abs(n_1)),
    const_term_0(const_term_0), const_term_1(const_term_1)
    {
        unsigned i;
        double x_, y_, logstep;
        vector<double> x, y;
    
        // LOWER RANGE
        for(i = 0, x_ = x_fit00, logstep = pow(x_fit01/x_fit00, 1./n_fit);
            i < n_fit; i++, x_*= logstep)
        {
            x.push_back(x_);
            y_ = (n_0 > 0) ? Interp_obj::eval(x_) : 1./Interp_obj::eval(x_);
            if (!const_term_0) y_ /= x_;
            y.push_back(y_);
        }
        if (const_term_0) polynomialfit(n_fit, abs(n_0), x.data(), y.data(), coeff_0.data());
        else polynomialfit(n_fit, abs(n_0)-1, x.data(), y.data(), coeff_0.data());

        // UPPER RANGE
        x.clear();
        y.clear();
        for(i = 0, x_ = x_fit10, logstep = pow(x_fit11/x_fit10, 1./n_fit);
        i < n_fit; i++, x_*= logstep)
        {
            x.push_back(x_);
            y_ = (n_1 > 0) ? Interp_obj::eval(x_) : 1./Interp_obj::eval(x_);
            if (!const_term_1) y_ /= x_;
            y.push_back(y_);
        }
        if (const_term_1) polynomialfit(n_fit, abs(n_1), x.data(), y.data(), coeff_1.data());
        else polynomialfit(n_fit, abs(n_1)-1, x.data(), y.data(), coeff_1.data());
    }

    // METHODS
    double eval(double x) const
    {
        double val;
        if (x < x_inter0)
        {
            if (const_term_0) val = gsl_poly_eval(coeff_0.data(), abs(n_0), x);
            else val = x*gsl_poly_eval(coeff_0.data(), abs(n_0)-1, x);
            return (n_0 > 0) ? val : 1./val;
        }
        else if (x <= x_inter1) return Interp_obj::eval(x);
        else
        {
            if (const_term_1) val = gsl_poly_eval(coeff_1.data(), abs(n_1), x);
            else val = x*gsl_poly_eval(coeff_1.data(), abs(n_1)-1, x);
            return (n_1 > 0) ? val : 1./val;
        }
    }

private:
    // VARIABLES
    const int n_0, n_1;
    const double x_inter0, x_inter1, x_fit00, x_fit01, x_fit10, x_fit11;
    vector<double> coeff_0, coeff_1;
    const bool const_term_0, const_term_1;
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

void gen_corr_func_binned_gsl(const double x_min, const double x_max, const Data_x_y<double>& pwr_spec_binned, Data_x_y<double>* corr_func_binned)
{
    Data_x_y<double> pwr_spec_binned_cp = pwr_spec_binned; // prevent overwriting when 

    // printf("Allocationg space for interpolation function P(k)\n");
    // Interp_obj P_k(pwr_spec_binned_cp);

    // printf("Allocationg space for integration via [QAWO adaptive integration for oscillatory functions]\n");
    // const double k_min = min(pwr_spec_binned_cp.x);
    // const double k_max = max(pwr_spec_binned_cp.x) / 4;
    // xi_integrand_param my_param;
    // my_param.P_k = &P_k;
    // Integr_obj_qawo<xi_integrand_param> xi_r(&xi_integrand, &my_param, k_min, k_max, 1000, 25);

    // printf("Computing correlation function via [QAWO adaptive integration for oscillatory functions]...\n");
    // const unsigned int N = corr_func_binned->size();
    // const double lin_bin = (x_max - x_min)/N;
	// double r;
	// for(unsigned i = 0; i < N; i++){
    //     r = x_min + i*lin_bin;
    //     corr_func_binned->x[i] = r;
    //     corr_func_binned->y[i] = xi_r(r);
    // }


    /* TEST OF INTERPOLATION */
    double k_min = min(pwr_spec_binned_cp.x);
    const double k_max = max(pwr_spec_binned_cp.x);
    Extrap_obj P_k(pwr_spec_binned_cp, k_min*5, k_max/5, 2, k_min, k_min*5, -3, k_max/20, k_max/10, 10, false, true);


    corr_func_binned->resize(200);
    const double log_bin = pow(1000*k_max/k_min, 1./corr_func_binned->size());
    double k;
    k_min *=sqrt(log_bin)/10;

    #pragma omp parallel for private(k)
	for (unsigned j = 0; j < corr_func_binned->size(); j++){
		k = k_min*pow(log_bin, j);
		corr_func_binned->x[j] = k;
        corr_func_binned->y[j] = P_k.eval(k);
    }
    /* END OF TEST*/
}