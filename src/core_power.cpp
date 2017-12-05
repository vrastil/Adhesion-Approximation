
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
    double operator()(void* params)
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
    double operator()(void* params)
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
    double operator()(double r, void* params)
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
    double operator()(double r, void* params)
    {
        gsl_integration_qawo_table_set(wf, r, L, GSL_INTEG_SINE);
        F.params = params;
        gsl_integration_qawf(&F, a, epsabs, limit, w, wc, wf, &result, &error);
        return result;
    }
protected:
    gsl_integration_workspace* wc;
};

/**
 * @class:	ODE_Solver
 * @brief:	Explicit embedded Runge-Kutta Prince-Dormand (8, 9) method.
 */


ODE_Solver::ODE_Solver(int (* function) (double, const double[], double[], void*), size_t dim, void* params,
            const double hstart, const double epsabs, const double epsrel):
    sys({function, NULL, dim, params}),
    d(gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, hstart, epsabs, epsrel)) {}

ODE_Solver::~ODE_Solver()
{
    gsl_odeiv2_driver_free (d);
}

void ODE_Solver::update(double &t, const double t1, double y[])
{
    status = gsl_odeiv2_driver_apply (d, &t, t1, y);
    if (status == GSL_SUCCESS) return;
    else if (status == GSL_EMAXITER) throw runtime_error("maximum number of steps reached");
    else if (status == GSL_ENOPROG) throw runtime_error("the step size droped below minimum value");
    else if (status == GSL_EBADFUNC) throw runtime_error("uknown error in 'gsl_odeiv2'");
}

double transfer_function_2(double k, const Cosmo_Param& parameters)
{
    if (k == 0) return 1.;

	const double q = k / (parameters.Omega_m*parameters.h);
	double T_k =	log(1+2.34*q)/(2.34*q)*
					pow(1 + 3.89*q + pow(16.2*q, 2.) + pow(5.47*q, 3.) + pow(6.71*q, 4.)
					, -1./4.);
	return pow(T_k, 2.);
}

double power_spectrum_T(double k, const Cosmo_Param& parameters)
{
	if (k == 0) return 0;
	const double ns = parameters.ns;
	return pow(k, ns)*transfer_function_2(k, parameters);
}

double power_spectrum_scale_free(double k, const Cosmo_Param& parameters)
{
	if (k == 0) return 0;
	const double ns = parameters.ns;
	return pow(k, ns);
}

double flat_power_spectrum(double k, const Cosmo_Param& parameters)
{
	return 1;
}

double single_power_spectrum_T(double k, const Cosmo_Param& parameters)
{
	if ((k > 0.01) and (k < 0.04)) return power_spectrum_T(k, parameters);
	else return 0.;
}

double power_spectrum(double k, const Cosmo_Param& parameters)
{
    const double A = parameters.A;
    const double supp = parameters.k2_G ? exp(-k*k/parameters.k2_G) : 1;
	switch (parameters.pwr_type)
	{
		case e_power_spec::power_law_T: return supp*A*power_spectrum_T(k, parameters);
		case e_power_spec::power_law: return supp*A*power_spectrum_scale_free(k, parameters);
		case e_power_spec::flat: return supp*A*flat_power_spectrum(k, parameters);
		case e_power_spec::single: return supp*A*single_power_spectrum_T(k, parameters);
		default: return supp*A*power_spectrum_T(k, parameters);
	}
}


double power_spectrum_s8(double k, void* parameters)
{
	return	k*k/(2.*PI*PI) // spherical factor
			*power_spectrum(k, *static_cast<const Cosmo_Param*>(parameters)) // P(k)
			*pow(3.*gsl_sf_bessel_j1(k*8)/(k*8), 2.); // window function R = 8 Mpc/h
}

void norm_pwr_gsl(Cosmo_Param* cosmo)
{
    /* Normalize the power spectrum */
    Integr_obj_qagiu s8(&power_spectrum_s8, 0, 0, 1e-7, 1000);
	cosmo->A = pow(cosmo->sigma8, 2)/s8(cosmo);
}

void norm_pwr_ccl(Cosmo_Param* cosmo)
{
    /* Normalize the power spectrum */
    int status = 0;
    ccl_sigma8(cosmo->cosmo, &status);
}

void norm_pwr(Cosmo_Param* cosmo)
{
    printf("Initializing given power spectrum...\n");
    if (cosmo->pwr_type_i < 4) norm_pwr_gsl(cosmo);
    else norm_pwr_ccl(cosmo);
}

double hubble_param(double a, const Cosmo_Param& cosmo)
{
    const double Om = cosmo.Omega_m;
    const double OL = cosmo.Omega_L();
    return sqrt(Om*pow(a, -3) + OL);
}

double Omega_lambda(double a, const Cosmo_Param& cosmo)
{
    const double Om = cosmo.Omega_m;
    const double OL = cosmo.Omega_L();
    return OL/(Om*pow(a, -3) + OL);
}

double growth_factor_integrand(double a, void* parameters)
{    
    return pow(a*hubble_param(a, static_cast<Cosmo_wrapper*>(parameters)->cosmo), -3);
}

double growth_factor(double a, const Cosmo_Param& cosmo)
{
    #ifdef CCL_GROWTH
    int status = 0;
    return ccl_growth_factor(cosmo.cosmo, a, &status);
    #else
    if (a == 0) return 0;
    Integr_obj_qag D(&growth_factor_integrand, 0, a, 0, 1e-12, 1000, GSL_INTEG_GAUSS61);
    Cosmo_wrapper param(cosmo);
    return hubble_param(a, cosmo)*D(&param)/cosmo.D_norm;
    #endif
}

double growth_factor(double a, void* parameters)
{
    return growth_factor(a, static_cast<Cosmo_wrapper*>(parameters)->cosmo);
}

double ln_growth_factor(double log_a, void* parameters)
{
    return log(growth_factor(exp(log_a), static_cast<Cosmo_wrapper*>(parameters)->cosmo));
}

double growth_rate(double a, const Cosmo_Param& cosmo)
{
    #ifdef CCL_GROWTH
    int status = 0;
    return ccl_growth_rate(cosmo.cosmo, a, &status);
    #else
    if (a == 0) return 1;
    Cosmo_wrapper param(cosmo);
    gsl_function F = {&ln_growth_factor, &param};
    double dDda, error;
    gsl_deriv_central(&F, log(a), 1e-12, &dDda, &error);
    return dDda;
    #endif
}

double growth_change(double a, const Cosmo_Param& cosmo)
{ // normal derivative dD/da, not logarithmic one
    Cosmo_wrapper param(cosmo);
    gsl_function F = {&growth_factor, &param};
    double dDda, error;
    if (a == 0) gsl_deriv_forward(&F, a, 1e-12, &dDda, &error);
    else gsl_deriv_central(&F, a, 1e-12, &dDda, &error);
    return dDda;
}

double lin_pow_spec(double k, const Cosmo_Param& cosmo, double a)
{
    if (cosmo.pwr_type_i < 4){
        return power_spectrum(k, cosmo);
    } else {
        int status = 0;
        return ccl_linear_matter_power(cosmo.cosmo, k*cosmo.h, a, &status)*pow(cosmo.h, 3);
    }
}

double lin_pow_spec(double k, const Cosmo_Param& cosmo)
{
    return lin_pow_spec(k, cosmo, 1.);
}

double find_pk_max(double k, void* parameters)
{
    return -lin_pow_spec(k, *static_cast<Cosmo_Param*>(parameters));
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
 * @brief:	Steffen interpolation of data [x, y]
 */

template<unsigned N>
Interp_obj::Interp_obj(const Data_Vec<double, N>& data): is_init(true)
{
    acc = gsl_interp_accel_alloc ();
    spline = gsl_spline_alloc (gsl_interp_steffen, data.size());
    gsl_spline_init (spline, data[0].data(), data[1].data(), data.size());

    x_min = data[0].front();
    x_max = data[0].back();
}

template<unsigned N>
void Interp_obj::init(const Data_Vec<double, N>& data)
{
    is_init = true;
    acc = gsl_interp_accel_alloc ();
    spline = gsl_spline_alloc (gsl_interp_steffen, data.size());
    gsl_spline_init (spline, data[0].data(), data[1].data(), data.size());

    x_min = data[0].front();
    x_max = data[0].back();
}

Interp_obj::~Interp_obj()
{
    if (is_init)
    {
        gsl_spline_free(spline);
        gsl_interp_accel_free (acc);
    }
}

double Interp_obj::operator()(double x) const{ return gsl_spline_eval(spline, x, acc); }

/**
 * @class:	Extrap_Pk
 * @brief:	linear interpolation of data [k, P(k)] within 'useful' range
            fit to primordial P_i(k) below the 'useful' range
            fit to Padé approximant R [0/3] above the 'useful' range
 */

template<unsigned N>
Extrap_Pk::Extrap_Pk(const Data_Vec<double, N>& data, const Sim_Param& sim, const unsigned m_l, const unsigned n_l,
    const unsigned m_u, const unsigned n_u):
    Interp_obj(data), cosmo(sim.cosmo)
{
    // LOWER RANGE -- fit linear power spectrum to data[m_l:n_l)
    #ifndef SWIG
    printf("Fitting amplitude of P(k) in lower range.\n");
    #endif
    k_min = data[0][m_l]; // first k in data
    fit_lin(data, m_l, n_l, A_low);

    // UPPER RANGE -- fit Ak^ns to data[m_u,n_u)
    #ifndef SWIG
    printf("Fitting amplitude of P(k) in upper range.\n");
    #endif
    k_max =  data[0][n_u-1]; // last k in data
    fit_power_law(data, m_u, n_u, A_up, n_s);
}

template<unsigned N>
Extrap_Pk::Extrap_Pk(const Data_Vec<double, N>& data, const Sim_Param& sim, const unsigned m_l, const unsigned n_u):
    Extrap_Pk(data, sim,
        // fit over first and last half of a decade
        m_l, m_l + sim.out_opt.bins_per_decade/2,
        n_u - sim.out_opt.bins_per_decade/2, n_u
    ) {}

template<unsigned N>
Extrap_Pk::Extrap_Pk(const Data_Vec<double, N>& data, const Sim_Param& sim):
    Extrap_Pk(data, sim,
        // trust simulation up to HALF k_nq
        0, get_nearest(sim.other_par.nyquist.at("particle")/2., data[0]) + 1
    ) {}

template<unsigned N>
void Extrap_Pk::fit_lin(const Data_Vec<double, N>& data, const unsigned m, const unsigned n, double &A)
{
    // FIT linear power spectrum to data[m:n)
    // fit 'P(k) = A * P_lin(k)' via A
    // for N = 2 perform non-weighted least-square fitting
    // for N = 3 use data[2] as sigma, w = 1/sigma^2
    double Pk, err;
    vector<double> Pk_res, w;
    Pk_res.reserve(n-m);
    if (N == 3){
        w.reserve(n-m);
    }
    vector<double> A_vec(n-m, 1);

    for(unsigned i = m; i < n; i++){
        Pk = lin_pow_spec(data[0][i], cosmo);
        Pk_res.push_back(data[1][i] / Pk);
        if (N == 3) w.push_back(pow(Pk/data[2][i], 2));
    }
    double A_sigma2, sumsq;

    if (!isfinite(A)) throw runtime_error("Encountered error during fitting");

    if (N == 3) gsl_fit_wmul(A_vec.data(), 1, w.data(), 1, Pk_res.data(), 1, n-m, &A, &A_sigma2, &sumsq);
    else gsl_fit_mul(A_vec.data(), 1, Pk_res.data(), 1, n-m, &A, &A_sigma2, &sumsq);

    #ifndef SWIG
    printf("\t[%sfit A = %.1e, err = %.2f\%]\n", N == 3 ? "weighted-" : "", A, 100*sqrt(A_sigma2)/A);
    #endif
}

template<unsigned N>
void Extrap_Pk::fit_power_law(const Data_Vec<double, N>& data, const unsigned m, const unsigned n, double &A, double &n_s)
{
    // FIT scale-free power spectrum to data[m,n)
    // fit 'log P(k) = log A + n_s * log k' via A, n_s
    // for N = 2 perform non-weighted least-square fitting
    // for N = 3 use data[2] as sigma, w = 1/sigma^2
    vector<double> k, Pk, w;
    double err;
    k.reserve(n-m);
    Pk.reserve(n-m);
    if (N == 3){
        w.reserve(n-m);
    }
    for(unsigned i = m; i < n; i++){
        k.push_back(log(data[0][i]));
        Pk.push_back(log(data[1][i]));
        if (N == 3)w.push_back(pow(data[1][i]/data[2][i], 2)); // weight = 1/sigma^2, approx log(1+x) = x for x rel. error
    }
    double cov00, cov01, cov11, sumsq;
    if (N == 3) gsl_fit_wlinear(k.data(), 1, w.data(), 1, Pk.data(), 1, n-m, &A, &n_s, &cov00, &cov01, &cov11, &sumsq);
    else gsl_fit_linear(k.data(), 1, Pk.data(), 1, n-m, &A, &n_s, &cov00, &cov01, &cov11, &sumsq);

    if (!isfinite(A)) throw runtime_error("Encountered error during fitting");

    A = exp(A); // log A => A
    #ifndef SWIG
    printf("\t[%sfit A = %.1e, err = %.2f\%, n_s = %.3f, err = %.2f\%, corr = %.2f\%]\n",
            N == 3 ? "weighted-" : "", A, 100*sqrt(cov00), n_s,  100*sqrt(cov11)/abs(n_s), 100*cov01/sqrt(cov00*cov11));
   #endif
}

double Extrap_Pk::operator()(double k) const
{
    if (k < k_min)
    {
        return A_low*lin_pow_spec(k, cosmo);
    }
    else if (k <= k_max) return Interp_obj::operator()(k);
    else return A_up*pow(k, n_s);
}

template <class P>
struct xi_integrand_param
{
    xi_integrand_param(double r, const P& P_k): r(r), P_k(P_k) {}
    double r;
    const P& P_k;
};

struct xi_integrand_param_lin
{
    double r;
    const Cosmo_Param& cosmo;
};

/**
 * @brief integrand for correlation function when weight-function 'sin(kr)' is used in integration
 */
template <class P>
double xi_integrand_W(double k, void* params){
    xi_integrand_param<P>* my_par = (xi_integrand_param<P>*) params;
    const double r = my_par->r;
    const P& P_k = my_par->P_k;
    return 1/(2*PI*PI)*k/r*P_k(k);
};

double xi_integrand_W_lin(double k, void* params){
    xi_integrand_param_lin* my_par = (xi_integrand_param_lin*) params;
    const double r = my_par->r;
    const double P_k = lin_pow_spec(k, my_par->cosmo);
    return 1/(2*PI*PI)*k/r*P_k;
};
/**
 * @brief integrand for correlation function when non-weighted integration is used
 */
template <class P>
double xi_integrand_G(double k, void* params){
    xi_integrand_param<P>* my_par = (xi_integrand_param<P>*) params;
    const double r = my_par->r;
    double j0;
    if (k*r < 1e-6) j0 = 1 - k*k*r*r/6.;
    else j0 = sin(k*r) / (k*r);
    const P& P_k = my_par->P_k;
    return 1/(2*PI*PI)*k*k*j0*P_k(k);
};

template <class P, class T> // both callable with 'operator()(double)'
void gen_corr_func_binned_gsl(const Sim_Param &sim, const P& P_k, Data_Vec<double, 2>* corr_func_binned, T& xi_r)
{
    const double x_min = sim.other_par.x_corr.lower;
    const double x_max = sim.other_par.x_corr.upper;

    xi_integrand_param<P> my_param(0, P_k);

    const double lin_bin = 10./sim.out_opt.points_per_10_Mpc;
    unsigned req_size = (unsigned)ceil((x_max - x_min)/lin_bin);
    corr_func_binned->resize(req_size);

	double r;
	for(unsigned i = 0; i < req_size; i++){
        r = x_min + i*lin_bin;
        my_param.r = r;
        (*corr_func_binned)[0][i] = r;
        (*corr_func_binned)[1][i] = xi_r(r, &my_param);
    }
}

template <class T>
void gen_corr_func_binned_gsl_lin(const Sim_Param &sim, double a, Data_Vec<double, 2>* corr_func_binned, T& xi_r)
{
    const double x_min = sim.other_par.x_corr.lower;
    const double x_max = sim.other_par.x_corr.upper;

    xi_integrand_param_lin my_param = {0, sim.cosmo};

    const double lin_bin = 10./sim.out_opt.points_per_10_Mpc;
    unsigned req_size = (unsigned)ceil((x_max - x_min)/lin_bin);
    corr_func_binned->resize(req_size);

    const double D2 = pow(growth_factor(a, sim.cosmo), 2);
	double r;
	for(unsigned i = 0; i < req_size; i++){
        r = x_min + i*lin_bin;
        my_param.r = r;
        (*corr_func_binned)[0][i] = r;
        (*corr_func_binned)[1][i] = D2*xi_r(r, &my_param);
    }
}

#define EPSABS 1e-9
#define EPSREL 1e-6

template<class P> // everything callable P_k(k)
void gen_corr_func_binned_gsl_qawf(const Sim_Param &sim, const P& P_k, Data_Vec<double, 2>* corr_func_binned)
{
    printf("Computing correlation function via GSL integration QAWF of extrapolated power spectrum...\n");

    Integr_obj_qawf xi_r(&xi_integrand_W<P>, 0, EPSABS,  4000, 50);
    gen_corr_func_binned_gsl(sim, P_k, corr_func_binned, xi_r);
}

void gen_corr_func_binned_gsl_qawf_lin(const Sim_Param &sim, double a, Data_Vec<double, 2>* corr_func_binned)
{
    printf("Computing correlation function via GSL integration QAWF of linear power spectrum...\n");

    Integr_obj_qawf xi_r(&xi_integrand_W_lin, 0, EPSABS,  4000, 50);
    gen_corr_func_binned_gsl_lin(sim, a, corr_func_binned, xi_r);
}

double  get_max_Pk(Sim_Param* sim)
{
    const gsl_min_fminimizer_type * T = gsl_min_fminimizer_brent;
    gsl_min_fminimizer * s = gsl_min_fminimizer_alloc (T);
    gsl_function F;
    F.function = &find_pk_max;
    F.params = &sim->cosmo;
    double k_l = sim->other_par.k_print.lower;
    double k_u = sim->other_par.k_print.upper;
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

// maybe replace function templates by class template? better to maintain, less variability
template Interp_obj::Interp_obj(const Data_Vec<double, 2>&);
template Interp_obj::Interp_obj(const Data_Vec<double, 3>&);
template void Interp_obj::init(const Data_Vec<double, 2>&);
template void Interp_obj::init(const Data_Vec<double, 3>&);
template Extrap_Pk::Extrap_Pk(const Data_Vec<double, 2>&, const Sim_Param&);
template Extrap_Pk::Extrap_Pk(const Data_Vec<double, 3>&, const Sim_Param&);
template Extrap_Pk::Extrap_Pk(const Data_Vec<double, 2>&, const Sim_Param&, const unsigned, const unsigned);
template Extrap_Pk::Extrap_Pk(const Data_Vec<double, 3>&, const Sim_Param&, const unsigned, const unsigned);

template Extrap_Pk::Extrap_Pk(const Data_Vec<double, 2>&, const Sim_Param&, const unsigned, const unsigned,
              const unsigned, const unsigned);
template Extrap_Pk::Extrap_Pk(const Data_Vec<double, 3>&, const Sim_Param&, const unsigned, const unsigned,
              const unsigned, const unsigned);
template void Extrap_Pk::fit_lin(const Data_Vec<double, 2>&, const unsigned m, const unsigned n, double& A);
template void Extrap_Pk::fit_lin(const Data_Vec<double, 3>&, const unsigned m, const unsigned n, double& A);
template void Extrap_Pk::fit_power_law(const Data_Vec<double, 2>&, const unsigned m, const unsigned n, double& A, double& n_s);
template void Extrap_Pk::fit_power_law(const Data_Vec<double, 3>&, const unsigned m, const unsigned n, double& A, double& n_s);

template void gen_corr_func_binned_gsl_qawf(const Sim_Param&, const Extrap_Pk&, Data_Vec<double, 2>*);
template void gen_corr_func_binned_gsl_qawf(const Sim_Param&, const Extrap_Pk_Nl&, Data_Vec<double, 2>*);