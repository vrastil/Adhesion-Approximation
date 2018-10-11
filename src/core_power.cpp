/**
 * @brief definitions of cosmological functions like power spectrum, growth, etc.
 * 
 * @file core_power.cpp
 * @author Michal Vrastil
 * @date 2018-06-23
 */

#include "core_power.h"
#include "templates/class_data_vec.hpp"
#include "params.hpp"

#include <ccl_defs.h>
#include <ccl_background.h>
#include <ccl_power.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fit.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_spline.h>

/*****************************//**
 * PRIVATE FUNCTIONS DEFINITIONS *
 *********************************/
namespace{
    
/**
 * @brief general integrand with callable function and one parameter
 * 
 * @tparam P callable object with 'operator()(double)'
 */
template <class P>
struct integrand_param
{
    integrand_param(double r, const P& P_k): r(r), P_k(P_k) {}
    double r;
    const P& P_k;
};

/**
 * @note GSL functions are mostly witten for double so even with different
 * PRECISION set these are computed in double precision 
 */

/**
 * @class Integr_obj
 * @brief basic integration object, wrapper for GSL integration functions
 * 
 */
class Integr_obj
{
public:
    // CONSTRUCTOR
    Integr_obj(double(*f) (double, void*),  const double a, const double b,
               const double epsabs, const double epsrel, const size_t limit):
    a(a), b(b), L(b-a), epsabs(epsabs), epsrel(epsrel),  limit(limit)
    {
        w = gsl_integration_workspace_alloc(limit);
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

/**
 * @class Integr_obj_qag
 * @brief QAG adaptive integration
 * 
 * The QAG algorithm is a simple adaptive integration procedure.
 */
class Integr_obj_qag : public Integr_obj
{
public:
    // CONSTRUCTOR
    Integr_obj_qag(double(*f) (double, void*),  const double a, const double b,
        const double epsabs, const double epsrel, size_t limit, int key):
    Integr_obj(f, a, b, epsabs, epsrel, limit), key(key) {}

    // METHODS
    double operator()(double r, void* params)
    {
        F.params = params;
        gsl_errno = gsl_integration_qag(&F, a, b, epsabs, epsrel, limit, key, w, &result, &error);
        if (gsl_errno) throw std::runtime_error("GSL integration error: " + std::string(gsl_strerror(gsl_errno)));
        else return result;
    }

    double operator()(void* params)
    {
        return this->operator()(0, params);
    }

protected:
    // VARIABLES
    int key;
};

/**
 * @class Integr_obj_qagiu
 * @brief QAGI adaptive integration on infinite intervals
 * 
 */
class Integr_obj_qagiu : public Integr_obj
{
public:
    // CONSTRUCTOR
    Integr_obj_qagiu(double(*f) (double, void*), const double a,
        const double epsabs, const double epsrel, size_t limit):
    Integr_obj(f, a, 0, epsabs, epsrel, limit) {}

    // METHODS
    double operator()(double r, void* params)
    {
        F.params = params;
        gsl_errno = gsl_integration_qagiu(&F, a, epsabs, epsrel, limit, w, &result, &error);
        if (gsl_errno) throw std::runtime_error("GSL integration error: " + std::string(gsl_strerror(gsl_errno)));
        else return result;
    }

    double operator()(void* params)
    {
        return this->operator()(0, params);
    }
};

/**
 * @class Integr_obj_qawo
 * @brief QAWO adaptive integration for oscillatory functions
 * 
 * The QAWO algorithm is designed for integrands with an oscillatory factor,
 * \f$\sin(\omega x)\f$ or \f$\cos(\omega x)\f$.
 */
class Integr_obj_qawo : public Integr_obj
{
public:
    // CONSTRUCTOR
    Integr_obj_qawo(double(*f) (double, void*), const double a, const double b,
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
        gsl_errno = gsl_integration_qawo(&F, a, epsabs, epsrel, limit, w, wf, &result, &error);
        if (gsl_errno) throw std::runtime_error("GSL integration error: " + std::string(gsl_strerror(gsl_errno)));
        else return result;
    }
protected:
    gsl_integration_qawo_table* wf;
};

/**
 * @class Integr_obj_qawf
 * @brief QAWF adaptive integration for Fourier integrals
 * 
 */
class Integr_obj_qawf : public Integr_obj_qawo
{
public:
    // CONSTRUCTOR
    Integr_obj_qawf(double(*f) (double, void*), const double a,
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
        gsl_errno = gsl_integration_qawf(&F, a, epsabs, limit, w, wc, wf, &result, &error);
        if (gsl_errno) throw std::runtime_error("GSL integration error: " + std::string(gsl_strerror(gsl_errno)));
        else return result;
    }
protected:
    gsl_integration_workspace* wc;
};

/**
 * @class:	ODE_Solver
 * @brief:	Explicit embedded Runge-Kutta Prince-Dormand (8, 9) method.
 */

class ODE_Solver
{
public:
    ODE_Solver(int (* function) (double t, const double y[], double dydt[], void * params), size_t dim, void* params,
               const double hstart = 1e-6, const double epsabs = 1e-6, const double epsrel = 0):
        sys({function, NULL, dim, params}),
        d(gsl_odeiv2_driver_alloc_y_new(&sys, gsl_odeiv2_step_rk8pd, hstart, epsabs, epsrel)) {}

    ~ODE_Solver()
    {
        gsl_odeiv2_driver_free (d);
    }

    void update(double &t, const double t1, double y[])
    {
        status = gsl_odeiv2_driver_apply (d, &t, t1, y);
        if (status) throw std::runtime_error("GSL ODE error: " + std::string(gsl_strerror(status)));
    }

    int status;
    gsl_odeiv2_system sys;
    gsl_odeiv2_driver* d;
};

template<typename T>
T hubble_param(T a, const Cosmo_Param& cosmo)
{   // hubble normalize to H(a = 1) == 1
    return sqrt(cosmo.Omega_m/pow(a, 3) + cosmo.Omega_L());
}

double growth_factor_integrand(double a, void* params)
{    
    return pow(a*hubble_param(a, *static_cast<const Cosmo_Param*>(params)), -3);
}

double growth_factor(double a, void* params)
{
    return growth_factor(a, *static_cast<const Cosmo_Param*>(params));
}

double ln_growth_factor(double log_a, void* parameters)
{
    return log(growth_factor(exp(log_a), parameters));
}

template<typename T>
int get_nearest(const T val, const std::vector<T>& vec)
{   // assume data in 'vec' are ordered, vec[0] < vec[1] < ...
    return lower_bound(vec.begin(), vec.end(), val) - vec.begin();
}

/**
 * @brief integrand for correlation function when weight-function 'sin(kr)' is used in integration
 * 
 * @tparam P callable object with 'operator()(double)'
 * @param k 
 * @param params 
 * @return double 
 */
template <class P>
double xi_integrand_W(double k, void* params){
    integrand_param<P>* my_par = (integrand_param<P>*) params;
    const double r = my_par->r;
    const P& P_k = my_par->P_k;
    return 1/(2*PI*PI)*k/r*P_k(k);
};

/**
 * @brief integrand for correlation function when non-weighted integration is used
 * 
 * @tparam P callable object with 'operator()(double)'
 * @param k 
 * @param params 
 * @return double 
 */
template <class P>
double xi_integrand_G(double k, void* params){
    integrand_param<P>* my_par = (integrand_param<P>*) params;
    const double kr = k*my_par->r;
    double j0 = kr < 1e-6 ? 1 - kr*kr/6. : sin(kr) / (kr);
    const P& P_k = my_par->P_k;
    return 1/(2*PI*PI)*k*k*j0*P_k(k);
};

/**
 * @brief integrand for amplitude of mass fluctuaton when non-weighted integration is used
 * 
 * @tparam P 
 * @param k 
 * @param params 
 * @return double 
 */
template <class P>
double sigma_integrand_G(double k, void* params){
    integrand_param<P>* my_par = (integrand_param<P>*) params;
    const double kr = k*my_par->r;
    const double w = kr < 0.1 ?
        1.-0.1*kr*kr+0.003571429*kr*kr*kr*kr
        -6.61376E-5*kr*kr*kr*kr*kr*kr
        +7.51563E-7*kr*kr*kr*kr*kr*kr*kr*kr
        :
        3.*(sin(kr) - kr*cos(kr))/(kr*kr*kr);

    const P& P_k = my_par->P_k;
    return 1/(2*PI*PI)*k*k*w*w*P_k(k);
};

/**
 * @brief compute given function in linear range and store results
 * 
 * @tparam P callable object with 'operator()(double)'
 * @tparam T callable object with 'operator()(double)'
 * @param x_min starting radius of function evaluation
 * @param x_max last radius of function evaluation
 * @param lin_bin linear step between function evaluations
 * @param P_k power spectrum (callable)
 * @param fce_binned object to store binned function
 * @param fce_r callable function
 */

template <class P, class T> 
void gen_fce_r_binned_gsl(const double x_min, const double x_max, const double lin_bin, const P& P_k, Data_Vec<FTYPE_t, 2>& fce_binned, T& fce_r)
{
    const size_t req_size = (size_t)ceil((x_max - x_min)/lin_bin);
    
    fce_binned.resize(req_size);
    integrand_param<P> my_param(0, P_k);

	for(size_t i = 0; i < req_size; i++){
        double r = x_min + i*lin_bin;
        my_param.r = r;
        fce_binned[0][i] = r;
        fce_binned[1][i] = fce_r(r, &my_param);
    }
}

/**
 * @brief compute correlation function and store results
 * 
 * @tparam P callable object with 'operator()(double)'
 * @tparam T callable object with 'operator()(double)'
 * @param sim simulation parameters
 * @param P_k power spectrum (callable)
 * @param corr_func_binned object to store binned correlation function
 * @param xi_r correlation function at given r (callable)
 */

template <class P, class T> // both callable with 'operator()(double)'
void gen_corr_func_binned_gsl(const Sim_Param &sim, const P& P_k, Data_Vec<FTYPE_t, 2>& corr_func_binned, T& xi_r)
{
    const double x_min = sim.other_par.x_corr.lower;
    const double x_max = sim.other_par.x_corr.upper;
    const double lin_bin = 10./sim.out_opt.points_per_10_Mpc;
    gen_fce_r_binned_gsl(x_min, x_max, lin_bin, P_k, corr_func_binned, xi_r);
}

/**
 * @brief compute amplitude of density fluctuation and store results
 * 
 * @tparam P callable object with 'operator()(double)'
 * @tparam T callable object with 'operator()(double)'
 * @param sim simulation parameters
 * @param P_k power spectrum (callable)
 * @param sigma_binned object to store binned amplitude of density fluctuation
 * @param sigma_r amplitude of density fluctuation at given r (callable)
 */

template <class P, class T> // both callable with 'operator()(double)'
void gen_sigma_func_binned_gsl(const Sim_Param &sim, const P& P_k, Data_Vec<FTYPE_t, 2>& sigma_binned, T& sigma_r)
{
    const double x_min = sim.other_par.x_corr.lower;
    const double x_max = sim.other_par.x_corr.upper;
    const double lin_bin = 10./sim.out_opt.points_per_10_Mpc;
    gen_fce_r_binned_gsl(x_min, x_max, lin_bin, P_k, sigma_binned, sigma_r);
}

constexpr double GSL_EPSABS = 1e-9; ///< absolute error
constexpr size_t GSL_LIMIT = 4000; ///< max. number of subintervals for adaptive integration
constexpr size_t GSL_N = 50; ///< max. number of bisections of integration interval (QAWO / QAWF)

} ///< end of anonymous namespace (private definitions)

/****************************//**
 * PUBLIC FUNCTIONS DEFINITIONS *
 ********************************/

void norm_pwr(Cosmo_Param& cosmo)
{
    #ifndef LESSINFO
    std::cout << "Initializing CCL power spectrum...\n";
    #endif
    int status = 0;
    ccl_sigma8(cosmo.cosmo, &status);
    if (status) throw std::runtime_error(cosmo.cosmo->status_message);
}

FTYPE_t norm_growth_factor(const Cosmo_Param& cosmo)
{
    Integr_obj_qag D(&growth_factor_integrand, 0, 1, 0, 1e-12, 1000, GSL_INTEG_GAUSS61);
    return D(static_cast<void*>(cosmo));
}

FTYPE_t growth_factor(FTYPE_t a, const Cosmo_Param& cosmo)
{
    // D(0) == 0; D(1) == 1, do not check (not often, better performance)
    if (!a) return 0;

    // try ccl range
    int status = 0;
    FTYPE_t D_ccl = ccl_growth_factor(cosmo.cosmo, a, &status);
    if (!status) return D_ccl;

    // integrate outside range
    Integr_obj_qag D(&growth_factor_integrand, 0, a, 0, 1e-12, 1000, GSL_INTEG_GAUSS61);
    return hubble_param(a, cosmo)*D(static_cast<void*>(cosmo))/cosmo.D_norm;
}

FTYPE_t growth_rate(FTYPE_t a, const Cosmo_Param& cosmo)
{
    // f(0) == 1
    if (!a) return 1;

    // try ccl range
    int status = 0;
    double f = ccl_growth_rate(cosmo.cosmo, a, &status);
    if (!status) return f;
    
    // logarithmic derivative outside range
    gsl_function F = {&ln_growth_factor, static_cast<void*>(cosmo)};
    double error;
    status = gsl_deriv_central(&F, log(a), 1e-12, &f, &error);
    if (!status) return f;
    else throw std::runtime_error("GSL ODE error: " + std::string(gsl_strerror(status)));
}

FTYPE_t growth_change(FTYPE_t a, const Cosmo_Param& cosmo)
{
    if (a){
        // compute with growth rate
        return growth_rate(a, cosmo)*growth_factor(a, cosmo)/a;
    }
    else{
        gsl_function F = {&growth_factor, static_cast<void*>(cosmo)};
        double dDda, error;
        int status = gsl_deriv_forward(&F, a, 1e-12, &dDda, &error);
        if (!status) return dDda;
        else throw std::runtime_error("GSL ODE error: " + std::string(gsl_strerror(status)));
    }
}

FTYPE_t Omega_lambda(FTYPE_t a, const Cosmo_Param& cosmo)
{
    // try ccl range
    int status = 0;
    FTYPE_t OL = ccl_omega_x(cosmo.cosmo, a, ccl_species_l_label, &status);
    if(!status) return OL;

    // compute outside range
    OL = cosmo.Omega_L();
    return OL/(cosmo.Omega_m/pow(a, 3) + OL);
}

FTYPE_t lin_pow_spec(FTYPE_t a, FTYPE_t k, const Cosmo_Param& cosmo)
{
    int status = 0;
    FTYPE_t pk = ccl_linear_matter_power(cosmo.cosmo, k*cosmo.h, a, &status)*pow(cosmo.h, 3);
    if (!status) return pk;
    status = 0;
    pk = ccl_linear_matter_power(cosmo.cosmo, k*cosmo.h, 1, &status)*pow(cosmo.h, 3);
    FTYPE_t D = growth_factor(a, cosmo);
    if (!status) return D*D*pk;
    else throw std::runtime_error(cosmo.cosmo->status_message);
}

FTYPE_t non_lin_pow_spec(FTYPE_t a, FTYPE_t k, const Cosmo_Param& cosmo)
{
    if (a < 1e-3) return lin_pow_spec(a, k, cosmo);
    int status = 0;
    FTYPE_t pk = ccl_nonlin_matter_power(cosmo.cosmo, k*cosmo.h, a, &status)*pow(cosmo.h, 3);
    if (!status) return pk;
    else throw std::runtime_error(cosmo.cosmo->status_message);
}

template <typename T, size_t N>
void Interp_obj::init(const Data_Vec<T, N>& data)
{
    is_init = true;
    acc = gsl_interp_accel_alloc ();
    spline = gsl_spline_alloc (gsl_interp_steffen, data.size());

    std::vector<double> tmp_x(data[0].begin(), data[0].end());
    std::vector<double> tmp_y(data[1].begin(), data[1].end());

    gsl_spline_init (spline, tmp_x.data(), tmp_y.data(), data.size());

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

template <typename T, size_t N>
Extrap_Pk<T, N>::Extrap_Pk(const Data_Vec<T, N>& data, const Sim_Param& sim, const size_t m_l, const size_t n_l,
    const size_t m_u, const size_t n_u):
    cosmo(sim.cosmo)
{
    this->init(data);//< initialize Interp_obj
    // LOWER RANGE -- fit linear power spectrum to data[m_l:n_l)
    #ifndef LESSINFO
    printf("Fitting amplitude of P(k) in lower range.\n");
    #endif
    k_min = data[0][m_l]; // first k in data
    fit_lin(data, m_l, n_l, A_low);

    // UPPER RANGE -- fit Ak^ns to data[m_u,n_u)
    #ifndef LESSINFO
    printf("Fitting amplitude of P(k) in upper range.\n");
    #endif
    k_max =  data[0][n_u-1]; // last k in data
    fit_power_law(data, m_u, n_u, A_up, n_s);
}

template <typename T, size_t N>
Extrap_Pk<T, N>::Extrap_Pk(const Data_Vec<T, N>& data, const Sim_Param& sim, const size_t m_l, const size_t n_u):
    Extrap_Pk(data, sim,
        // fit over first and last half of a decade
        m_l, m_l + sim.out_opt.bins_per_decade/2,
        n_u - sim.out_opt.bins_per_decade/2, n_u
    ) {}

template <typename T, size_t N>
Extrap_Pk<T, N>::Extrap_Pk(const Data_Vec<T, N>& data, const Sim_Param& sim):
    Extrap_Pk(data, sim,
        // trust simulation up to HALF k_nq
        0, get_nearest(sim.other_par.nyquist.at("particle")/2, data[0]) + 1
    ) {}

template <typename T, size_t N>
void Extrap_Pk<T, N>::fit_lin(const Data_Vec<T, N>& data, const size_t m, const size_t n, double &A)
{
    // FIT linear power spectrum to data[m:n)
    // fit 'P(k) = A * P_lin(k)' via A
    // for N = 2 perform non-weighted least-square fitting
    // for N = 3 use data[2] as sigma, w = 1/sigma^2
    double Pk;
    std::vector<double> Pk_res, w;
    Pk_res.reserve(n-m);
    if (N == 3){
        w.reserve(n-m);
    }
    std::vector<double> A_vec(n-m, 1);

    for(size_t i = m; i < n; i++){
        Pk = lin_pow_spec(1, data[0][i], cosmo);
        Pk_res.push_back(data[1][i] / Pk);
        if (N == 3) w.push_back(pow2(Pk/data[2][i]));
    }
    double A_sigma2, sumsq;

    int gsl_errno;
    if (N == 3) gsl_errno = gsl_fit_wmul(A_vec.data(), 1, w.data(), 1, Pk_res.data(), 1, n-m, &A, &A_sigma2, &sumsq);
    else gsl_errno = gsl_fit_mul(A_vec.data(), 1, Pk_res.data(), 1, n-m, &A, &A_sigma2, &sumsq);
    if (gsl_errno) throw std::runtime_error("GSL integration error: " + std::string(gsl_strerror(gsl_errno)));

    #ifndef LESSINFO
    printf("\t[%sfit A = %.1e, err = %.2f%%]\n", N == 3 ? "weighted-" : "", A, 100*sqrt(A_sigma2)/A);
    #endif
}

template <typename T, size_t N>
void Extrap_Pk<T, N>::fit_power_law(const Data_Vec<T, N>& data, const size_t m, const size_t n, double &A, double &n_s)
{
    // FIT scale-free power spectrum to data[m,n)
    // fit 'log P(k) = log A + n_s * log k' via A, n_s
    // for N = 2 perform non-weighted least-square fitting
    // for N = 3 use data[2] as sigma, w = 1/sigma^2
    std::vector<double> k, Pk, w; //< need double for GSL
    k.reserve(n-m);
    Pk.reserve(n-m);
    if (N == 3){
        w.reserve(n-m);
    }
    for(size_t i = m; i < n; i++){
        k.push_back(log(data[0][i]));
        Pk.push_back(log(data[1][i]));
        if (N == 3)w.push_back(pow2(data[1][i]/data[2][i])); // weight = 1/sigma^2, approx log(1+x) = x for x rel. error
    }
    double cov00, cov01, cov11, sumsq;
    int gsl_errno;
    if (N == 3) gsl_errno = gsl_fit_wlinear(k.data(), 1, w.data(), 1, Pk.data(), 1, n-m, &A, &n_s, &cov00, &cov01, &cov11, &sumsq);
    else gsl_errno = gsl_fit_linear(k.data(), 1, Pk.data(), 1, n-m, &A, &n_s, &cov00, &cov01, &cov11, &sumsq);
    if (gsl_errno) throw std::runtime_error("GSL integration error: " + std::string(gsl_strerror(gsl_errno)));

    A = exp(A); // log A => A
    #ifndef LESSINFO
    printf("\t[%sfit A = %.1e, err = %.2f%%, n_s = %.3f, err = %.2f%%, corr = %.2f%%]\n",
            N == 3 ? "weighted-" : "", A, 100*sqrt(cov00), n_s,  100*sqrt(cov11)/std::abs(n_s), 100*cov01/sqrt(cov00*cov11));
    #endif
}

template <typename T, size_t N>
double Extrap_Pk<T, N>::operator()(double k) const
{
    if (k < k_min)
    {
        return A_low*lin_pow_spec(1, k, cosmo);
    }
    else if (k <= k_max) return Interp_obj::operator()(k);
    else return A_up*pow(k, n_s);
}

template <typename T, size_t N>
Extrap_Pk_Nl<T, N>::Extrap_Pk_Nl(const Data_Vec<T, N>& data, const Sim_Param &sim, T A_nl, T a_eff):
    Extrap_Pk<T, N>(data, sim), A_nl(A_nl), a_eff(a_eff), k_split(this->k_max) {}

template <typename T, size_t N>
double Extrap_Pk_Nl<T, N>::operator()(double k) const {
    if (k < k_split) return Extrap_Pk<T, N>::operator()(k);
    else return (1-A_nl)*lin_pow_spec(a_eff, k, this->cosmo) + A_nl*non_lin_pow_spec(a_eff, k, this->cosmo);
}

template<class P>
void gen_corr_func_binned_gsl_qawf(const Sim_Param &sim, const P& P_k, Data_Vec<FTYPE_t, 2>& corr_func_binned)
{
    #ifndef LESSINFO
    printf("Computing correlation function via GSL integration QAWF...\n");
    #endif
    Integr_obj_qawf xi_r(&xi_integrand_W<P>, 0, GSL_EPSABS,  GSL_LIMIT, GSL_N);
    gen_corr_func_binned_gsl(sim, P_k, corr_func_binned, xi_r);
}

void gen_corr_func_binned_gsl_qawf_lin(const Sim_Param &sim, FTYPE_t a, Data_Vec<FTYPE_t, 2>& corr_func_binned)
{
    auto P_k = [&](FTYPE_t k){ return lin_pow_spec(a, k, sim.cosmo); };
    gen_corr_func_binned_gsl_qawf(sim, P_k, corr_func_binned);
}

void gen_corr_func_binned_gsl_qawf_nl(const Sim_Param &sim, FTYPE_t a, Data_Vec<FTYPE_t, 2>& corr_func_binned)
{
    auto P_k = [&](FTYPE_t k){ return non_lin_pow_spec(a, k, sim.cosmo); };
    gen_corr_func_binned_gsl_qawf(sim, P_k, corr_func_binned);
}

template<class P>
void gen_sigma_binned_gsl_qawf(const Sim_Param &sim, const P& P_k, Data_Vec<FTYPE_t, 2>& sigma_binned)
{
    #ifndef LESSINFO
    printf("Computing mass fluctuations via GSL integration QAWF...\n");
    #endif
    Integr_obj_qagiu sigma_r(&sigma_integrand_G<P>, 0, GSL_EPSABS,  GSL_LIMIT, GSL_N);
    gen_sigma_func_binned_gsl(sim, P_k, sigma_binned, sigma_r);
}

void gen_sigma_func_binned_gsl_qawf_lin(const Sim_Param &sim, FTYPE_t a, Data_Vec<FTYPE_t, 2>& sigma_binned)
{
    auto P_k = [&](FTYPE_t k){ return lin_pow_spec(a, k, sim.cosmo); };
    gen_sigma_binned_gsl_qawf(sim, P_k, sigma_binned);
}

void gen_sigma_func_binned_gsl_qawf_nl(const Sim_Param &sim, FTYPE_t a, Data_Vec<FTYPE_t, 2>& sigma_binned)
{
    auto P_k = [&](FTYPE_t k){ return non_lin_pow_spec(a, k, sim.cosmo); };
    gen_sigma_binned_gsl_qawf(sim, P_k, sigma_binned);
}

/**********************//**
 * EXPLICIT INSTANTIATION *
 **************************/

template void Interp_obj::init(const Data_Vec<FTYPE_t, 2>& data);
template void Interp_obj::init(const Data_Vec<FTYPE_t, 3>& data);

template class Extrap_Pk<FTYPE_t, 2>;
template class Extrap_Pk<FTYPE_t, 3>;
template class Extrap_Pk_Nl<FTYPE_t, 2>;
template class Extrap_Pk_Nl<FTYPE_t, 3>;

template void gen_corr_func_binned_gsl_qawf(const Sim_Param&, const Extrap_Pk<FTYPE_t, 2>&, Data_Vec<FTYPE_t, 2>&);
template void gen_corr_func_binned_gsl_qawf(const Sim_Param&, const Extrap_Pk<FTYPE_t, 3>&, Data_Vec<FTYPE_t, 2>&);
template void gen_corr_func_binned_gsl_qawf(const Sim_Param&, const Extrap_Pk_Nl<FTYPE_t, 2>&, Data_Vec<FTYPE_t, 2>&);
template void gen_corr_func_binned_gsl_qawf(const Sim_Param&, const Extrap_Pk_Nl<FTYPE_t, 3>&, Data_Vec<FTYPE_t, 2>&);

template void gen_sigma_binned_gsl_qawf(const Sim_Param&, const Extrap_Pk<FTYPE_t, 2>&, Data_Vec<FTYPE_t, 2>&);
template void gen_sigma_binned_gsl_qawf(const Sim_Param&, const Extrap_Pk<FTYPE_t, 3>&, Data_Vec<FTYPE_t, 2>&);
template void gen_sigma_binned_gsl_qawf(const Sim_Param&, const Extrap_Pk_Nl<FTYPE_t, 2>&, Data_Vec<FTYPE_t, 2>&);
template void gen_sigma_binned_gsl_qawf(const Sim_Param&, const Extrap_Pk_Nl<FTYPE_t, 3>&, Data_Vec<FTYPE_t, 2>&);
