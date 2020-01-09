/**
 * @brief various simulation parameters
 * 
 * @file params.hpp
 * @author Michal Vrastil
 * @date 2018-07-11
 */

#pragma once
#include "stdafx.h"
#include <ccl_defs.h>
#include <ccl_config.h>
#include <ccl_core.h>
#include <map>
#include "precision.hpp"

/**
 * @brief cosmological & CCL parameters
 * @class Cosmo_Param
 * 
 */
class Cosmo_Param
{
public:
    // CONSTRUCTORS, DESTRUCTOR
    void init(); ///< lazy constructor
    Cosmo_Param();
    ~Cosmo_Param();

    // CCL VARIABLES
    ccl_configuration config;
    ccl_cosmology* cosmo;

    // COSMOLOGY (flat LCDM)
    FTYPE_t ns, k2_G, sigma8;
    FTYPE_t Omega_m, Omega_b, H0, h;
    FTYPE_t Omega_c() const { return Omega_m - Omega_b; }
    FTYPE_t Omega_L() const { return 1 - Omega_m; }

    // TRUNCATION OF INITIAL POWER SPECTRUM
    bool truncated_pk = false;

    // PRECOMPUTED VALUES
    FTYPE_t D_norm;

    // DEALING WITH GSL 'void* param'
    explicit operator void*() const;
};

/**
 * @brief simulation box options
 * @struct Box_Opt
 * 
 */
struct Box_Opt {
    void init(const Cosmo_Param& cosmo);
    /* cmd args */
    size_t par_num_1d, mesh_num, mesh_num_pwr;
    FTYPE_t box_size;
    /* derived param*/
    size_t par_num, Ng, Ng_pwr;
    FTYPE_t mass_p_log; ///< logarithm of particle mass in \f$M_\odot\f$
};

/**
 * @brief integration options
 * @struct Integ_Opt
 * 
 */
struct Integ_Opt {
    void init();
    FTYPE_t z_in, z_out, db; ///< cmd args
    FTYPE_t a_in, a_out; ///< derived parameters
};


/**
 * @brief output options
 * @struct Out_Opt
 * 
 */
struct Out_Opt {
    void init();
    /* cmd args */
    size_t print_every, bins_per_decade, points_per_10_Mpc;
    std::string out_dir; //< where to save output of the simulation
    bool print_par_pos, print_dens, print_pwr, print_extrap_pwr, print_corr, print_vel_pwr;
    /* derived param*/
    bool get_rho, get_pwr, get_pk_extrap;
};


/**
 * @brief specify which approximations are run
 * @struct Comp_App
 * 
 */
struct Comp_App {
    /* cmd args */
    bool ZA, TZA, FF, FP, AA, FP_pp; //< approximations
    bool PM; //< PM codes
    bool chi, chi_ff; //< modified gravities
    void reset();
    bool is_ready();
};


/**
 * @brief approximations options
 * @struct App_Opt
 * 
 */
struct App_Opt {
    void init(const Box_Opt&);
    /* cmd args */
    FTYPE_t nu, rs;
    /* derived param*/
    FTYPE_t Hc, a, nu_dim;
    size_t M;
};


/**
 * @brief run options
 * @struct Run_Opt
 * 
 */
struct Run_Opt {
    void init();
    void reset();
    bool is_ready();
    bool simulate();
    /* cmd args */
    size_t nt, nt_fftw, mlt_runs;
    size_t seed;
    bool pair;
    /* other*/
    bool phase = true;
};

/**
 * @brief testing options
 * @struct Test_Opt
 * 
 */
struct Test_Opt {
    /* cmd args */
    FTYPE_t R_sphere, rho_sphere;
    size_t fine_sweeps, coarse_sweeps, max_steps, step_per_iter;
    size_t N_grid, N_min;
    bool verbose;

    /* derived param*/
    FTYPE_t rho_b;
};

// define Range outside because of SWIG
/**
 * @brief lower and upper boundaries
 * @struct Range
 * 
 */
struct Range { FTYPE_t lower, upper; };

/**
 * @brief other parameters
 * @struct Other_par
 * 
 */
struct Other_par {
    void init(const Box_Opt&);
    // k-range where to use (linear) interpolation and k-range in which print 'pwr_spec_extrap_*'
    ///range in which compute the correlation function 
    Range k_print, x_corr;
    std::map<std::string,FTYPE_t> nyquist; //< Nyquist frequencies of potential mesh, analyses mesh and particle separation
};

/**
 * @brief chameleon options
 * @struct Chi_Opt
 * 
 */
struct Chi_Opt {
    /* cmd args */
    FTYPE_t beta, n, phi;
    bool linear;
};

/**
 * @class:	Sim_Param
 * @brief:	class storing simulation parameters
 */
class Sim_Param
{
public:
    // CONSTRUCTOR
    Sim_Param(){}; //< from simpy and partial initialization
    Sim_Param(int ac, const char* const av[]); //< load from command line arguments
    Sim_Param(std::string file_name); //< load from sim_param.json file

    // VARIABLES
    Box_Opt box_opt;
    Integ_Opt integ_opt;
    Out_Opt out_opt;
    Comp_App comp_app;
    Cosmo_Param cosmo;
    App_Opt app_opt;
    Run_Opt run_opt;
    Other_par other_par;
    Chi_Opt chi_opt;
    Test_Opt test_opt;

	// METHODS
    void print_info(std::string out, std::string app) const;
	void print_info() const;
	FTYPE_t x_0() const{return box_opt.box_size/box_opt.mesh_num;}
    FTYPE_t x_0_pwr() const{return box_opt.box_size/box_opt.mesh_num_pwr;}
    bool simulate() { return run_opt.simulate(); }
    void reset();
    bool is_ready();
};