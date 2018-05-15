/**
 * @file:	params.hpp
 * @brief:	various parameter classes declarations
 */

#pragma once
#include "stdafx.h"
#include <ccl_config.h>
#include <ccl_core.h>
#include <json.hpp>
#include <map>
#include "precision.hpp"

/* COSMOLOGICAL & CCL PARAMETERS*/
class Cosmo_Param
{
public:
    // CONSTRUCTORS, DESTRUCTOR
    void init(); //< lazy constructor
    Cosmo_Param();
    ~Cosmo_Param();

    // CCL VARIABLES
    ccl_configuration config;
    ccl_cosmology* cosmo;

    // COSMOLOGY (flat LCDM)
    FTYPE_t A = 1, ns, k2_G, sigma8;
    FTYPE_t Omega_m, Omega_b, H0, h;
    FTYPE_t Omega_c() const { return Omega_m - Omega_b; }
    FTYPE_t Omega_L() const { return 1 - Omega_m; }

    // PRECOMPUTED VALUES
    FTYPE_t D_norm;

    // DEALING WITH GSL 'void* param'
    explicit operator void*() const;
};

/* SIMULATION BOX*/
struct Box_Opt {
    void init();
    /* cmd args */
    unsigned par_num_1d, mesh_num, mesh_num_pwr;
    FTYPE_t box_size;
    /* derived param*/
    unsigned par_num, Ng, Ng_pwr;
};


/* INTEGRATION */
struct Integ_Opt {
    void init();
    /* cmd args */
    FTYPE_t z_in, z_out, db;
    /* derived param*/
    FTYPE_t b_in, b_out;
};


/* OUTPUT */
struct Out_Opt {
    void init();
    /* cmd args */
    unsigned print_every, bins_per_decade, points_per_10_Mpc;
    std::vector<FTYPE_t> print_z; //< for which redshifts print output on top of print_every (optional)
    std::string out_dir; //< where to save output of the simulation
    bool print_par_pos, print_dens, print_pwr, print_extrap_pwr, print_corr, print_vel_pwr;
    /* derived param*/
    bool get_rho, get_pwr, get_pk_extrap;
};


/* APPROXIMATIONS */
struct Comp_App {
    /* cmd args */
    bool ZA, FF, FP, AA, FP_pp; //< approximations
    bool chi; //< modified gravities
};


/* APPROXIMATIONS */
struct App_Opt {
    void init(const Box_Opt&);
    /* cmd args */
    FTYPE_t nu, rs;
    /* derived param*/
    FTYPE_t Hc, a, nu_dim;
    unsigned M;
};


/* RUN */
struct Run_Opt {
    void init();
    bool simulate();
    /* cmd args */
    unsigned nt, mlt_runs;
    unsigned long seed;
    bool pair;        
    /* other*/
    bool phase;
};

/* TEST */
struct Test_Opt {
    FTYPE_t R_sphere, rho_sphere;
    unsigned N_grid;
};

// define Range outside because of SWIG
struct Range { FTYPE_t lower, upper; };

/* OTHER PARAMETERS */
struct Other_par {
    void init(const Box_Opt&);
    // k-range where to use (linear) interpolation and k-range in which print 'pwr_spec_extrap_*'
    ///range in which compute the correlation function 
    Range k_print, x_corr;
    std::map<std::string,FTYPE_t> nyquist; //< Nyquist frequencies of potential mesh, analyses mesh and particle separation
};

/* CHAMELEON */
struct Chi_Opt {
    /* cmd args */
    FTYPE_t beta, n, phi;
};

/**
 * @class:	Sim_Param
 * @brief:	class storing simulation parameters
 */

class Sim_Param
{
public:
    // CONSTRUCTOR
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
};

// interaction with json files

void to_json(nlohmann::json&, const Cosmo_Param&);
void to_json(nlohmann::json&, const Box_Opt&);
void to_json(nlohmann::json&, const Integ_Opt&);
void to_json(nlohmann::json&, const App_Opt&);
void to_json(nlohmann::json&, const Run_Opt&);
void to_json(nlohmann::json&, const Out_Opt&);

void from_json(const nlohmann::json&, Box_Opt&);
void from_json(const nlohmann::json&, Integ_Opt&);
void from_json(const nlohmann::json&, App_Opt&);
void from_json(const nlohmann::json&, Run_Opt&);
void from_json(const nlohmann::json&, Out_Opt&);
void from_json(const nlohmann::json&, Cosmo_Param&);