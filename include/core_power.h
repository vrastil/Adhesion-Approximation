/**
 * @brief handle cosmological functions like power spectrum, growth, etc.
 * 
 * @file core_power.h
 * @author Michal Vrastil
 * @date 2018-06-21
 */

#pragma once
#include "precision.hpp"
#include <gsl/gsl_spline.h>

/*******************//**
* FORWARD DECLARATIONS *
************************/

class Cosmo_Param; ///< declaration in params.hpp
class Sim_Param; ///< declaration in params.hpp
template <typename T, unsigned N> class Data_Vec; ///< declaration in class_data_vec.hpp

/**
 * @brief Initialize CCL power spectrum
 * 
 * @param cosmo cosmological parameters
 */
void norm_pwr(Cosmo_Param& cosmo);

/**
 * @brief when computing growth factor outside CCL range we need to normalize the growth factor; \f$D(a=1)\equiv1\f$
 * 
 * @param cosmo cosmological parameters
 * @return FTYPE_t unnormalized growth factor at 'a = 1'
 */
FTYPE_t norm_growth_factor(const Cosmo_Param& cosmo);

/**
 * @brief 
 * 
 * @param a 
 * @param cosmo cosmological parameters
 * @return FTYPE_t 
 */
FTYPE_t growth_factor(FTYPE_t a, const Cosmo_Param& cosmo);

/**
 * @brief 
 * 
 * @param a 
 * @param cosmo cosmological parameters
 * @return FTYPE_t 
 */
FTYPE_t growth_rate(FTYPE_t a, const Cosmo_Param& cosmo);

/**
 * @brief 
 * 
 * @param a 
 * @param cosmo cosmological parameters
 * @return FTYPE_t 
 */
FTYPE_t growth_change(FTYPE_t a, const Cosmo_Param& cosmo);

/**
 * @brief 
 * 
 * @param a 
 * @param cosmo cosmological parameters
 * @return FTYPE_t 
 */
FTYPE_t Omega_lambda(FTYPE_t a, const Cosmo_Param& cosmo);

/**
 * @brief 
 * 
 * @param a 
 * @param k 
 * @param cosmo cosmological parameters
 * @return FTYPE_t 
 */
FTYPE_t lin_pow_spec(FTYPE_t a, FTYPE_t k, const Cosmo_Param& cosmo);

/**
 * @brief 
 * 
 * @param a 
 * @param k 
 * @param cosmo cosmological parameters
 * @return FTYPE_t 
 */
FTYPE_t non_lin_pow_spec(FTYPE_t a, FTYPE_t k, const Cosmo_Param& cosmo);

/**
 * @class:	Interp_obj
 * @brief:	linear interpolation (Steffen) of data [x, y]
 */
class Interp_obj
{
public:
    Interp_obj(): is_init(false) {}
    ~Interp_obj();
    double operator()(double x) const;
    template <typename T, unsigned N>
    void init(const Data_Vec<T, N>& data);
    double x_min, x_max;

private:
    bool is_init;
    gsl_spline* spline;
    gsl_interp_accel* acc;
};

/**
 * @class:	Extrap_Pk
 * @brief:	linear interpolation of data [k, P(k)] within 'useful' range
 *          fit to primordial P_i(k) below the 'useful' range
 *          fit to Padé approximant R [0/3] above the 'useful' range
 *
 *  Steffen interpolation of data [k, P(k)] within range k_min, k_max
 *  fit to primordial P_i(k) below this range, fit A*k^ns above
 */
template <typename T, unsigned N>
class Extrap_Pk : public Interp_obj
{
public:
    Extrap_Pk(const Data_Vec<T, N>& data, const Sim_Param& sim);
    Extrap_Pk(const Data_Vec<T, N>& data, const Sim_Param& sim, const unsigned m_l, const unsigned n_u);
    Extrap_Pk(const Data_Vec<T, N>& data, const Sim_Param& sim, const unsigned m_l, const unsigned n_l,
              const unsigned m_u, const unsigned n_u);
    double operator()(double k) const;

    void fit_lin(const Data_Vec<T, N>& data, const unsigned m, const unsigned n, double& A);
    void fit_power_law(const Data_Vec<T, N>& data, const unsigned m, const unsigned n, double& A, double& n_s);

    double A_low; ///< amplitude of linear power in lower range
    const Cosmo_Param& cosmo;
    double A_up, n_s; ///< scale-free power spectrum in upper range
    T k_min, k_max; ///< interpolation range
};

/**
 * @class:	Extrap_Pk_Nl
 * @brief:	creates Extrapolate object (linear power spectrum) from data and store non-linear parameters
            call 'operator()(k)' based on k_split (upper range of the linear)
 */
template <typename T, unsigned N>
class Extrap_Pk_Nl : public Extrap_Pk<T, N>
{
public:
    Extrap_Pk_Nl(const Data_Vec<T, N>& data, const Sim_Param &sim, T A_nl, T a_eff);
    const T A_nl, a_eff, k_split;
    double operator()(double k) const;
};

/**
 * @brief compute correlation function and store results
 * 
 * @tparam P callable object with 'operator()(double)'
 * @param sim simulation parameters
 * @param P_k power spectrum (callable)
 * @param corr_func_binned object to store binned correlation function
 */
template<class P>
void gen_corr_func_binned_gsl_qawf(const Sim_Param &sim, const P& P_k, Data_Vec<FTYPE_t, 2>& corr_func_binned);

/**
 * @brief compute linear correlation function and store results
 * 
 * @param sim simulation parameters
 * @param a scale factor
 * @param corr_func_binned object to store binned correlation function
 */
void gen_corr_func_binned_gsl_qawf_lin(const Sim_Param &sim, FTYPE_t a, Data_Vec<FTYPE_t, 2>& corr_func_binned);

/**
 * @brief compute non-linear correlation function and store results
 * 
 * @param sim simulation parameters
 * @param a scale factor
 * @param corr_func_binned object to store binned correlation function
 */
void gen_corr_func_binned_gsl_qawf_nl(const Sim_Param &sim, FTYPE_t a, Data_Vec<FTYPE_t, 2>& corr_func_binned);

/**
 * @brief compute amplitude of density fluctuation and store results
 * 
 * @tparam P callable object with 'operator()(double)'
 * @param sim simulation parameters
 * @param P_k power spectrum (callable)
 * @param sigma_binned object to store binned correlation function
 */
template<class P>
void gen_sigma_binned_gsl_qawf(const Sim_Param &sim, const P& P_k, Data_Vec<FTYPE_t, 2>& sigma_binned);

/**
 * @brief compute linear amplitude of density fluctuation and store results
 * 
 * @param sim simulation parameters
 * @param a scale factor
 * @param sigma_binned object to store binned correlation function
 */
void gen_sigma_func_binned_gsl_qawf_lin(const Sim_Param &sim, FTYPE_t a, Data_Vec<FTYPE_t, 2>& sigma_binned);

/**
 * @brief compute non-linear amplitude of density fluctuation and store results
 * 
 * @param sim simulation parameters
 * @param a scale factor
 * @param sigma_binned object to store binned correlation function
 */
void gen_sigma_func_binned_gsl_qawf_nl(const Sim_Param &sim, FTYPE_t a, Data_Vec<FTYPE_t, 2>& sigma_binned);
