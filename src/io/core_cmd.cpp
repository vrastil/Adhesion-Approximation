/**
 * @brief command line arguments manipulation
 * 
 * @file core_cmd.cpp
 * @author Michal Vrastil
 * @date 2018-07-11
 */

#include "stdafx.h"
#include <boost/any.hpp>
#include <boost/program_options.hpp>
#include <iomanip>
#include <fstream>
#include "params.hpp"
#include "core_cmd.h"
#include "config.h"

namespace po = boost::program_options;

/****************************//**
 * PUBLIC FUNCTIONS DEFINITIONS *
 ********************************/

void handle_cmd_line(int ac, const char* const av[], Sim_Param& sim){
    std::string config_file;
    bool no_run = false; // do not run approximations if '--help= or '--version' cmd option
    unsigned int trans_func_cmd, matter_pwr_cmd, baryons_pwr_cmd, mass_func_cmd;
    FTYPE_t k2_G;
    // options ONLY on command line
    po::options_description generic("Generic options");
    generic.add_options()
        ("help,h", "produce this help message")			
        ("config,c", po::value<std::string>(&config_file)->default_value("INPUT.cfg"), "configuration file name (optional)")
    	("version,v", "print current version of the program")
        ;
        
    // options both on command line	and in configuration file
    po::options_description config_mesh("Simulation box options");
    config_mesh.add_options()
        ("mesh_num,m", po::value<size_t>(&sim.box_opt.mesh_num)->default_value(128), "number of mesh cells per dimension (potential)")
        ("mesh_num_pwr,M", po::value<size_t>(&sim.box_opt.mesh_num_pwr)->default_value(256), "number of mesh cells per dimension (power spectrum)")
        ("par_num,p", po::value<size_t>(&sim.box_opt.par_num_1d)->default_value(128), "number of particles per dimension")
        ("box_size,L", po::value<FTYPE_t>(&sim.box_opt.box_size)->default_value(512, "512"), "box size in units of Mpc/h")
        ;
        
    po::options_description config_integ("Integration options");
    config_integ.add_options()
        ("redshift,z", po::value<FTYPE_t>(&sim.integ_opt.z_in)->default_value(200.), "redshift at the start of the simulation")
        ("redshift_0,Z", po::value<FTYPE_t>(&sim.integ_opt.z_out)->default_value(10.), "redshift at the end of the simulation")
        ("time_step,a", po::value<FTYPE_t>(&sim.integ_opt.db)->default_value(0.1, "0.1"), "dimensionless time-step (scale factor)")
        ;
    
    po::options_description config_output("Output options");
    config_output.add_options()
        ("print_every", po::value<size_t>(&sim.out_opt.print_every)->default_value(1, "1"), "save particle positions and power spectrum "
                                                                                        "every n-th step, set 0 for no printing")
        ("pwr_bins", po::value<size_t>(&sim.out_opt.bins_per_decade)->default_value(30), "number of bins per decade in power spectrum")
        ("corr_pt", po::value<size_t>(&sim.out_opt.points_per_10_Mpc)->default_value(10), "number of points per 10 Mpc in correlation function")
        ("out_dir,o", po::value<std::string>(&sim.out_opt.out_dir)->default_value("output/"), "output folder name")
        ("print_par_pos", po::value<bool>(&sim.out_opt.print_par_pos)->default_value(false), "print particles positions")
        ("print_dens", po::value<bool>(&sim.out_opt.print_dens)->default_value(false), "print density map and histogram")
        ("print_pwr", po::value<bool>(&sim.out_opt.print_pwr)->default_value(false), "print power spectrum")
        ("print_extrap_pwr", po::value<bool>(&sim.out_opt.print_extrap_pwr)->default_value(false), "print extrapolated power spectrum")
        ("print_corr", po::value<bool>(&sim.out_opt.print_corr)->default_value(false), "print correlation function")
        ("print_vel_pwr", po::value<bool>(&sim.out_opt.print_vel_pwr)->default_value(false), "print velocity power spectrum")
        ;
    
    po::options_description config_app("Approximations");
    config_app.add_options()
        ("comp_ZA", po::value<bool>(&sim.comp_app.ZA)->default_value(false), "compute Zeldovich approximation")
        ("comp_TZA", po::value<bool>(&sim.comp_app.TZA)->default_value(false), "compute Truncated Zeldovich approximation")
        ("comp_FF", po::value<bool>(&sim.comp_app.FF)->default_value(false), "compute Frozen-flow approximation")
        ("comp_FP", po::value<bool>(&sim.comp_app.FP)->default_value(false), "compute Frozen-potential approximation")
        ("comp_PM", po::value<bool>(&sim.comp_app.PM)->default_value(false), "compute Particle-mesh approximation")
        ("comp_AA", po::value<bool>(&sim.comp_app.AA)->default_value(false), "compute Adhesion approximation")
        ("comp_FP_pp", po::value<bool>(&sim.comp_app.FP_pp)->default_value(false),
                            "compute Frozen-potential approximation (particle-particle interaction)")
        ;
        
    po::options_description config_power("Cosmological parameters");
    config_power.add_options()
        ("Omega_b,B", po::value<FTYPE_t>(&sim.cosmo.Omega_b)->default_value(0.05, "0.05"), "density of baryons relative to the critical density")
        ("Omega_m,C", po::value<FTYPE_t>(&sim.cosmo.Omega_m)->default_value(0.25, "0.25"), "density of CDM relative to the critical density")
        ("Hubble,H", po::value<FTYPE_t>(&sim.cosmo.H0)->default_value(67, "67"), "Hubble constant in units of km/s/Mpc")
        ("transfer_function", po::value<unsigned int>(&trans_func_cmd)->default_value(3), "transfer function type")
        ("matter_power_spectrum", po::value<unsigned int>(&matter_pwr_cmd)->default_value(1), "matter power spectrum type")
        ("baryons_power_spectrum", po::value<unsigned int>(&baryons_pwr_cmd)->default_value(0), "baryons power spectrum type")
        ("mass_function", po::value<unsigned int>(&mass_func_cmd)->default_value(2), "mass function type")
        ("n_s,n", po::value<FTYPE_t>(&sim.cosmo.ns)->default_value(1.), "spectral index of the scale-free power spectrum")
        ("sigma8,s", po::value<FTYPE_t>(&sim.cosmo.sigma8)->default_value(1.), "normalization of the power spectrum at R = 8 Mpc/h")
        ("smoothing_k,k", po::value<FTYPE_t>(&k2_G)->default_value(0.),
                            "smoothing wavenumber of TZA in units of h/Mpc, set 0 for ZA")
        ;

    po::options_description config_run("Run options");
    config_run.add_options()
        ("num_thread,t", po::value<size_t>(&sim.run_opt.nt)->default_value(0), "number of threads the program will use, set 0 for max. available")
        ("num_thread_fftw", po::value<size_t>(&sim.run_opt.nt_fftw)->default_value(0), "number of threads the program will use for FFTW, set 0 for max. available")
        ("seed", po::value<size_t>(&sim.run_opt.seed)->default_value(0), "seed to random number generator, use 0 for random")
        ("pair", po::value<bool>(&sim.run_opt.pair)->default_value(false), "if true run two simulations with opposite phases of random field")
        ("mlt_runs", po::value<size_t>(&sim.run_opt.mlt_runs)->default_value(1), "how many runs should be simulated (only if seed = 0)")
        ;

    po::options_description mod_grav("Modified Gravities");
    mod_grav.add_options()
        ("comp_chi", po::value<bool>(&sim.comp_app.chi)->default_value(false), "compute chameleon gravity (frozen-potential)")
        ("comp_chi_ff", po::value<bool>(&sim.comp_app.chi_ff)->default_value(false), "compute chameleon gravity (frozen-flow)")
        ;
    
    po::options_description config_cham("Chameleon parameters");
    config_cham.add_options()
        ("chi_beta", po::value<FTYPE_t>(&sim.chi_opt.beta)->default_value(1/sqrt(6), "(1/6)^1/2"), "coupling constant")
        ("chi_n", po::value<FTYPE_t>(&sim.chi_opt.n)->default_value(0.5, "1/2"), "chameleon power-law potential exponent,\n0 < n < 1")
        ("chi_phi", po::value<FTYPE_t>(&sim.chi_opt.phi)->default_value(1E-6, "1E-6"), "screening potential")
        ("comp_chi_lin", po::value<bool>(&sim.chi_opt.linear)->default_value(false), "use only linear prediction in k-space")
        ;  
    mod_grav.add(config_cham);
    
    po::options_description config_test("Test options");
    config_test.add_options()
        ("R_sphere", po::value<FTYPE_t>(&sim.test_opt.R_sphere)->default_value(1., "3.0"), "radius of a sphere sitting in a vacuum")
        ("rho_sphere", po::value<FTYPE_t>(&sim.test_opt.rho_sphere)->default_value(2.7, "1.0"), "density  of a sphere sitting in a vacuum")
        ("N_grid", po::value<size_t>(&sim.test_opt.N_grid)->default_value(64, "64"), "finest grid size")
        ("N_min", po::value<size_t>(&sim.test_opt.N_min)->default_value(2, "2"), "coarsest grid size")
        ("step_per_iter", po::value<size_t>(&sim.test_opt.step_per_iter)->default_value(5, "5"), "print chameleon results after X steps")
        ("fine_sweeps", po::value<size_t>(&sim.test_opt.fine_sweeps)->default_value(3, "3"), "number of sweeps on fine grid")
        ("coarse_sweeps", po::value<size_t>(&sim.test_opt.coarse_sweeps)->default_value(3, "3"), "number of sweeps on coarse grids")
        ("max_steps", po::value<size_t>(&sim.test_opt.max_steps)->default_value(1, "1"), "max number of V-cycles")
        ("verbose", po::value<bool>(&sim.test_opt.verbose)->default_value(true, "true"), "verbosity")
        ;


    po::options_description cmdline_options("\n" PROJECT " v" PROJECT_VERSION);//< store all normal parameters
    cmdline_options.add(generic).add(config_mesh).add(config_power).add(config_integ);
    cmdline_options.add(config_output).add(config_app).add(config_run);
    cmdline_options.add(mod_grav);

    config_test.add(cmdline_options);//< union of test parameters and normaln ones
    po::variables_map vm;
    //po::store(po::command_line_parser(ac, av).options(config_test).run(), vm);
    po::store(po::parse_command_line(ac, av, config_test), vm);
    po::notify(vm);

    if (vm.count("help")) {
        std::cout << std::setprecision(3) << cmdline_options;
        no_run = true;
    }
    else if (vm.count("version")) {
        std::cout << PROJECT " v" PROJECT_VERSION "\n";
        no_run = true;
    }
    else if (vm.count("config")) {
        std::ifstream ifs(config_file.c_str());
        if (ifs){
        BOOST_LOG_TRIVIAL(debug) << "Using configuration options defined in file " << config_file << " (given command line options have higher priority).";
            po::store(po::parse_config_file(ifs, config_test), vm);
            // po::store(po::parse_config_file(ifs, cmdline_options), vm);
            notify(vm);
        } else{
            BOOST_LOG_TRIVIAL(warning) << "Cannot open config file '" << config_file << "'. Using default and command lines values.";
        }
    }

    if (no_run) sim.reset();

    // !!!>>> THIS NEEDS TO BE AFTER ALL CALLS TO NOTIFY() <<<!!!
    sim.cosmo.config.transfer_function_method = static_cast<transfer_function_t>(trans_func_cmd);
    sim.cosmo.config.matter_power_spectrum_method = static_cast<matter_power_spectrum_t>(matter_pwr_cmd);
    sim.cosmo.config.baryons_power_spectrum_method = static_cast<baryons_power_spectrum_t>(baryons_pwr_cmd);
    sim.cosmo.config.mass_function_method = static_cast<mass_function_t>(mass_func_cmd);
    sim.cosmo.k2_G = k2_G*k2_G;
    sim.comp_app.TZA &= bool(sim.cosmo.k2_G);
}