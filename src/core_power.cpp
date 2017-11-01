
#include "stdafx.h"
#include "core.h"
#include "core_mesh.h"
#include "core_power.h"

using namespace std;

class Integr_obj
{
public:
    // CONSTRUCTOR
    Integr_obj(double(*f) (double, void*),  const double a, const double b,
               const double epsabs, const double epsrel, size_t limit):
    a(a), b(b), L(b-a), epsabs(epsabs), epsrel(epsrel),  limit(limit)
    {
        w = gsl_integration_workspace_alloc (limit);
        F.function = f;
    }
    // DESTRUCTOR
    ~Integr_obj()
    {
        gsl_integration_workspace_free (w);
    }
    // METHODS
    void set_a(double a_new) { a = a_new; L = a -b; }
    void set_b(double b_new) { b = b_new; L = a -b; }

protected:
    // VARIABLES
    double result, error;
    double a, b, L;
    double epsabs, epsrel;
    size_t limit;
    gsl_function F;
    gsl_integration_workspace* w;
    int gsl_errno;
};

class Integr_obj_qag : public Integr_obj
{
public:
    // CONSTRUCTOR
    Integr_obj_qag(double(*f) (double, void*),  double a, double b,
        const double epsabs, const double epsrel, size_t limit, int key):
    Integr_obj(f, a, b, epsabs, epsrel, limit), key(key) {}

    // METHODS
    double eval(void* params)
    {
        F.params = params;
        gsl_errno = gsl_integration_qag(&F, a, b, epsabs, epsrel, limit, key, w, &result, &error);
        if (gsl_errno) printf ("GSL integration error: %s\n", gsl_strerror (gsl_errno));
        return result;
    }
protected:
    // VARIABLES
    int key;
};

class Integr_obj_qagiu : public Integr_obj
{
public:
    // CONSTRUCTOR
    Integr_obj_qagiu(double(*f) (double, void*),  double a,
        const double epsabs, const double epsrel, size_t limit):
    Integr_obj(f, a, 0, epsabs, epsrel, limit) {}

    // METHODS
    double eval(void* params)
    {
        F.params = params;
        gsl_errno = gsl_integration_qagiu(&F, a, epsabs, epsrel, limit, w, &result, &error);
        return result;
    }
};

class Integr_obj_qawo : public Integr_obj
{
public:
    // CONSTRUCTOR
    Integr_obj_qawo(double(*f) (double, void*),  double a, double b,
                    const double epsabs, const double epsrel, size_t limit, size_t n):
    Integr_obj(f, a, b, epsabs, epsrel, limit)
    {
        wf = gsl_integration_qawo_table_alloc(1, 1, GSL_INTEG_SINE, n);
    }
    // DESTRUCTOR
    ~Integr_obj_qawo()
    {
        gsl_integration_qawo_table_free(wf);
    }
    // METHODS
    double eval(double r, void* params)
    {
        gsl_integration_qawo_table_set(wf, r, L, GSL_INTEG_SINE);
        F.params = params;
        gsl_integration_qawo(&F, a, epsabs, epsrel, limit, w, wf, &result, &error);
        return result;
    }
protected:
    gsl_integration_qawo_table* wf;
};

class Integr_obj_qawf : public Integr_obj_qawo
{
public:
    // CONSTRUCTOR
    Integr_obj_qawf(double(*f) (double, void*),  double a,
                    const double epsabs, size_t limit, size_t n):
    Integr_obj_qawo(f, a, 0, epsabs, 0, limit, n)
    {
        wc = gsl_integration_workspace_alloc (limit);
    }
    // DESTRUCTOR
    ~Integr_obj_qawf()
    {
        gsl_integration_workspace_free (wc);
    }
    // METHODS
    double eval(double r, void* params)
    {
        gsl_integration_qawo_table_set(wf, r, L, GSL_INTEG_SINE);
        F.params = params;
        gsl_integration_qawf(&F, a, epsabs, limit, w, wc, wf, &result, &error);
        return result;
    }
protected:
    gsl_integration_workspace* wc;
};

double transfer_function_2(double k, const Pow_Spec_Param& parameters)
{
    if (k == 0) return 1.;

	const double q = k / (parameters.Omega_m*parameters.h);
	double T_k =	log(1+2.34*q)/(2.34*q)*
					pow(1 + 3.89*q + pow(16.2*q, 2.) + pow(5.47*q, 3.) + pow(6.71*q, 4.)
					, -1./4.);
	return pow(T_k, 2.);
}

double power_spectrum_T(double k, const Pow_Spec_Param& parameters)
{
	if (k == 0) return 0;
	const double ns = parameters.ns;
	return pow(k, ns)*transfer_function_2(k, parameters);
}

double power_spectrum_scale_free(double k, const Pow_Spec_Param& parameters)
{
	if (k == 0) return 0;
	const double ns = parameters.ns;
	return pow(k, ns);
}

double flat_power_spectrum(double k, const Pow_Spec_Param& parameters)
{
	return 1;
}

double single_power_spectrum_T(double k, const Pow_Spec_Param& parameters)
{
	if ((k > 0.01) and (k < 0.04)) return power_spectrum_T(k, parameters);
	else return 0.;
}

double power_spectrum(double k, const Pow_Spec_Param& parameters)
{
    const double A = parameters.A;
    const double supp = parameters.k2_G ? exp(-k*k/parameters.k2_G) : 1;
	switch (parameters.pwr_type)
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
			*power_spectrum(k, *static_cast<const Pow_Spec_Param*>(parameters)) // P(k)
			*pow(3.*gsl_sf_bessel_j1(k*8)/(k*8), 2.); // window function R = 8 Mpc/h
}

void norm_pwr_gsl(Pow_Spec_Param* pwr_par)
{
    /* Normalize the power spectrum */
    Integr_obj_qagiu s8(&power_spectrum_s8, 0, 0, 1e-7, 1000);
	pwr_par->A = pow(pwr_par->sigma8, 2)/s8.eval(pwr_par);
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

double hubble_param(double a, const Pow_Spec_Param& pwr_par)
{
    const double Om = pwr_par.Omega_m;
    const double OL = pwr_par.Omega_L();
    return sqrt(Om*pow(a, -3) + OL);
}

struct growth_factor_integrand_param { const Pow_Spec_Param& pwr_par; };

double growth_factor_integrand(double a, void* parameters)
{    
    return pow(a*hubble_param(a, static_cast<growth_factor_integrand_param*>(parameters)->pwr_par), -3);
}

double growth_factor(double a, const Pow_Spec_Param& pwr_par)
{
    #ifdef CCL_GROWTH
    int status = 0;
    return ccl_growth_factor(pwr_par.cosmo, a, &status);
    #else
    if (a == 0) return 0;
    Integr_obj_qag D(&growth_factor_integrand, 0, a, 0, 1e-7, 1000, GSL_INTEG_GAUSS61);
    growth_factor_integrand_param param = {pwr_par};
    return hubble_param(a, pwr_par)*D.eval(&param)/pwr_par.D_norm;
    #endif
}

double growth_factor(double a, void* parameters)
{
    return growth_factor(a, static_cast<growth_factor_integrand_param*>(parameters)->pwr_par);
}

double ln_growth_factor(double log_a, void* parameters)
{
    return log(growth_factor(exp(log_a), static_cast<growth_factor_integrand_param*>(parameters)->pwr_par));
}

double growth_rate(double a, const Pow_Spec_Param& pwr_par)
{
    #ifdef CCL_GROWTH
    int status = 0;
    return ccl_growth_rate(pwr_par.cosmo, a, &status);
    #else
    if (a == 0) return 1;
    growth_factor_integrand_param param = {pwr_par};
    gsl_function F = {&ln_growth_factor, &param};
    double dDda, error;
    gsl_deriv_central(&F, log(a), 1e-6, &dDda, &error);
    return dDda;
    #endif
}

double growth_change(double a, const Pow_Spec_Param& pwr_par)
{ // normal derivative dD/da, not logarithmic one
    growth_factor_integrand_param param = {pwr_par};
    gsl_function F = {&growth_factor, &param};
    double dDda, error;
    if (a == 0) gsl_deriv_forward(&F, a, 1e-6, &dDda, &error);
    else gsl_deriv_central(&F, a, 1e-6, &dDda, &error);
    return dDda;
}

double lin_pow_spec(double k, const Pow_Spec_Param& pwr_par, double a)
{
    if (pwr_par.pwr_type < 4){
        return power_spectrum(k, pwr_par);
    } else {
        int status = 0;
        return ccl_linear_matter_power(pwr_par.cosmo, k*pwr_par.h, a, &status)*pow(pwr_par.h, 3);
    }
}

double lin_pow_spec(double k, const Pow_Spec_Param& pwr_par)
{
    return lin_pow_spec(k, pwr_par, 1.);
}

double find_pk_max(double k, void* parameters)
{
    return -lin_pow_spec(k, *static_cast<Pow_Spec_Param*>(parameters));
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
    Interp_obj(data), power(sim.power), n_s(0)
{
    unsigned m, n;
    // LOWER RANGE -- fit linear power spectrum to data[m:n)
    printf("Fitting amplitude of P(k) in lower range.\n");
    m = 10 + get_nearest(2*PI/sim.box_size, data.x); // discard first 10 values due to undersampling
    n = m + 20; // fit over 20 points
    k_min = data.x[n];
    fit_lin(data, m, n, A_low);
    // UPPER RANGE -- ffit linear power spectrum to data[m:n)
    printf("Fitting amplitude of P(k) in upper range.\n");
    k_max = sim.k_par.k_interp.upper;
    n = get_nearest(k_max, data.x) + 1;
    m = n - 5; // fit over last 5 values
    fit_lin(data, m, n, A_up);
}

Extrap_Pk::Extrap_Pk(const Data_x_y<double>& data, const Sim_Param& sim, const unsigned m_l, const unsigned n_l,
    const unsigned m_u, const unsigned n_u, double n_s):
    Interp_obj(data), power(sim.power), n_s(n_s)
{
    // LOWER RANGE -- fit linear power spectrum to data[m:n)
    printf("Fitting amplitude of P(k) in lower range.\n");
    k_min = data.x[n_l-1];
    fit_lin(data, m_l, n_l, A_low);

    // UPPER RANGE -- fit Ak^ns to data[m,n)
    printf("Fitting amplitude of P(k) in upper range.\n");
    k_max =  data.x[n_u-1];
    fit_power_law(data, m_u, n_u, A_up);
}

void Extrap_Pk::fit_lin(const Data_x_y<double>& data, const unsigned m, const unsigned n, double &A)
{
    // FIT linear power spectrum to data[m:n)
    vector<double> k, Pk;
    for(unsigned i = m; i < n; i++){
        k.push_back(lin_pow_spec(data.x[i], power));
        Pk.push_back(data.y[i]);
    }
    double A_sigma2, sumsq;
    gsl_fit_mul(k.data(), 1, Pk.data(), 1, n-m, &A, &A_sigma2, &sumsq);
    printf("\t[fit A = %e, sigma = %f, sumsq = %f]\n", A, sqrt(A_sigma2), sumsq);
}

void Extrap_Pk::fit_power_law(const Data_x_y<double>& data, const unsigned m, const unsigned n, double &A)
{
    // FIT Ak^ns to data[m,n)
    vector<double> k, Pk;
    for(unsigned i = m; i < n; i++){
        k.push_back(pow(data.x[i], n_s));
        Pk.push_back(data.y[i]);
    }
    double A_sigma2, sumsq;
    gsl_fit_mul(k.data(), 1, Pk.data(), 1, n-m, &A, &A_sigma2, &sumsq);
    printf("\t[fit A = %e, n_s = %.3f, sigma = %f, sumsq = %f]\n", A, n_s, sqrt(A_sigma2), sumsq);
}

double Extrap_Pk::eval(double k) const
{
    if (k < k_min)
    {
        return A_low*lin_pow_spec(k, power);
    }
    else if (k <= k_max) return Interp_obj::eval(k);
    else
    {
        if (n_s == 0) return A_up*lin_pow_spec(k, power);
        else return A_up*pow(k, n_s);
    }
}

struct xi_integrand_param
{
    double r;
    const Extrap_Pk* P_k;
};

struct xi_integrand_param_lin
{
    double r;
    const Pow_Spec_Param& pwr_par;
};

/**
 * @brief integrand for correlation function when weight-function 'sin(kr)' is used in integration
 */
double xi_integrand_W(double k, void* params){
    xi_integrand_param* my_par = (xi_integrand_param*) params;
    const double r = my_par->r;
    const Extrap_Pk* P_k = my_par->P_k;
    return 1/(2*PI*PI)*k/r*P_k->eval(k);
};

double xi_integrand_W_lin(double k, void* params){
    xi_integrand_param_lin* my_par = (xi_integrand_param_lin*) params;
    const double r = my_par->r;
    const double P_k = lin_pow_spec(k, my_par->pwr_par);
    return 1/(2*PI*PI)*k/r*P_k;
};

/**
 * @brief integrand for correlation function when non-weighted integration is used
 */
double xi_integrand_G(double k, void* params){
    xi_integrand_param* my_par = (xi_integrand_param*) params;
    const double r = my_par->r;
    double j0;
    if (k*r < 1e-3) j0 = 1 - k*k*r*r/6.;
    else j0 = sin(k*r) / (k*r);
    const Extrap_Pk* P_k = my_par->P_k;
    return 1/(2*PI*PI)*k*k*j0*P_k->eval(k);
};

template <class T>
void gen_corr_func_binned_gsl(const Sim_Param &sim, const Extrap_Pk& P_k, Data_x_y<double>* corr_func_binned, T* xi_r)
{
    const double x_min = sim.x_corr.lower;
    const double x_max = sim.x_corr.upper;

    xi_integrand_param my_param = {0, &P_k};

    const unsigned int N = corr_func_binned->size();
    const double lin_bin = (x_max - x_min)/N;
	double r;
	for(unsigned i = 0; i < N; i++){
        r = x_min + i*lin_bin;
        my_param.r = r;
        corr_func_binned->x[i] = r;
        corr_func_binned->y[i] = xi_r->eval(r, &my_param);
    }
}

template <class T>
void gen_corr_func_binned_gsl_lin(const Sim_Param &sim, double a, Data_x_y<double>* corr_func_binned, T* xi_r)
{
    const double x_min = sim.x_corr.lower;
    const double x_max = sim.x_corr.upper;

    xi_integrand_param_lin my_param = {0, sim.power};

    const unsigned int N = corr_func_binned->size();
    const double lin_bin = (x_max - x_min)/N;
    const double D2 = pow(growth_factor(a, sim.power), 2);
	double r;
	for(unsigned i = 0; i < N; i++){
        r = x_min + i*lin_bin;
        my_param.r = r;
        corr_func_binned->x[i] = r;
        corr_func_binned->y[i] = D2*xi_r->eval(r, &my_param);
    }
}

#define EPSABS 1e-7
#define EPSREL 1e-4

void gen_corr_func_binned_gsl_qagi(const Sim_Param &sim, const Extrap_Pk& P_k, Data_x_y<double>* corr_func_binned)
{
    printf("Computing correlation function via GSL integration QAGI of extrapolated power spectrum...\n");
    
}

void gen_corr_func_binned_gsl_qawo(const Sim_Param &sim, const Extrap_Pk& P_k, Data_x_y<double>* corr_func_binned)
{
    printf("Computing correlation function via GSL integration QAWO of extrapolated power spectrum...\n");
    const double k_min = 3e-3;
    const double k_max = 1e-1;
    Integr_obj_qawo xi_r(&xi_integrand_W, k_min, k_max, EPSABS, EPSREL,  4000, 50);
    gen_corr_func_binned_gsl(sim, P_k, corr_func_binned, &xi_r);
}

void gen_corr_func_binned_gsl_qawf(const Sim_Param &sim, const Extrap_Pk& P_k, Data_x_y<double>* corr_func_binned)
{
    printf("Computing correlation function via GSL integration QAWF of extrapolated power spectrum...\n");

    Integr_obj_qawf xi_r(&xi_integrand_W, 0, EPSABS,  4000, 50);
    gen_corr_func_binned_gsl(sim, P_k, corr_func_binned, &xi_r);
}

void gen_corr_func_binned_gsl_qawf_lin(const Sim_Param &sim, double a, Data_x_y<double>* corr_func_binned)
{
    printf("Computing correlation function via GSL integration QAWF of linear power spectrum...\n");

    Integr_obj_qawf xi_r(&xi_integrand_W_lin, 0, EPSABS,  4000, 50);
    gen_corr_func_binned_gsl_lin(sim, a, corr_func_binned, &xi_r);
}

double  get_max_Pk(Sim_Param* sim)
{
    const gsl_min_fminimizer_type * T = gsl_min_fminimizer_brent;
    gsl_min_fminimizer * s = gsl_min_fminimizer_alloc (T);
    gsl_function F;
    F.function = &find_pk_max;
    F.params = &sim->power;
    double k_l = sim->k_par.k_print.lower;
    double k_u = sim->k_par.k_print.upper;
    double k = (k_l + k_u) / 2.;
    gsl_min_fminimizer_set(s, &F, k, k_l, k_u);

    int status;
    int iter = 0, max_iter = 100;
    do
    {
        iter++;
        status = gsl_min_fminimizer_iterate (s);

        k = gsl_min_fminimizer_x_minimum (s);
        k_l = gsl_min_fminimizer_x_lower (s);
        k_u = gsl_min_fminimizer_x_upper (s);

        status = gsl_min_test_interval (k_l, k_u, 0, 1e-3);

    } while ((status == GSL_CONTINUE) && (iter < max_iter));
    gsl_min_fminimizer_free(s);
    return k;
}