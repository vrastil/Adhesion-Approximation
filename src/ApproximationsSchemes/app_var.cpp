/**
 * @brief classes handling approximations data
 * 
 * @file app_var.cpp
 * @author Michal Vrastil
 * @date 2018-07-11
 */
#include "app_var.hpp"
#include "core_out.h"
#include "core_app.h"
#include "core_mesh.h"

#include <algorithm>
#include <iomanip>
#include <sstream>

namespace {
const char *humanSize(uint64_t bytes){
	char const *suffix[] = {"B", "KB", "MB", "GB", "TB"};
	char length = sizeof(suffix) / sizeof(suffix[0]);

	size_t i = 0;
	double dblBytes = bytes;

	if (bytes > 1024) {
		for (i = 0; (bytes / 1024) > 0 && i<length-1; i++, bytes /= 1024)
			dblBytes = bytes / 1024.0;
	}

	static char output[200];
	sprintf(output, "%.02lf %s", dblBytes, suffix[i]);
	return output;
}

void print_mem(uint64_t memory_alloc)
{
    BOOST_LOG_TRIVIAL(debug) << "Allocated " << humanSize(memory_alloc) << " of memory.";
}

/**
 * @class:	Tracking
 * @brief:	class storing info about tracked particles
 */

class Tracking
{
public:
	// CONSTRUCTOR
	Tracking(size_t sqr_num_track_par, size_t par_num_per_dim)
    {
        set_par_ids(sqr_num_track_par, par_num_per_dim);
    }

    // METHODS
    size_t get_num_track_par() const
    {
        return par_ids.size();
    }

    size_t get_num_steps() const
    {
        return par_pos.size();
    }

	template <class T> void update_track_par(const std::vector<T>& particles)
    {
        std::vector<Particle_x<FTYPE_t>> par_pos_step;
        par_pos_step.reserve(par_ids.size());
        for (size_t i=0; i < par_ids.size(); i++){
            par_pos_step.emplace_back(particles[par_ids[i]].position);
        }
        par_pos.push_back(par_pos_step);
    }

    void print_track_par(const Sim_Param &sim, std::string out_dir, std::string suffix) const
    {
        out_dir += "par_cut/";
        std::string file_name = out_dir + "track_par_pos" + suffix + ".dat";
        Ofstream File(file_name);

        FTYPE_t x,y,z;
        const FTYPE_t x_0 = sim.x_0();
        BOOST_LOG_TRIVIAL(debug) << "Writing positons of " << par_ids.size() << " tracked particles into file " << file_name;
        File << "# This file contains positions of particles in units [Mpc/h].\n"
                "# x [Mpc/h]\tz [Mpc/h]";
        for (size_t i=0; i < par_ids.size(); i++){
            for (const auto& par_pos_step : par_pos){
                x = par_pos_step[i].position[0];
                y = par_pos_step[i].position[1];
                z = par_pos_step[i].position[2];
                File << x*x_0 << "\t" << z*x_0 << "\t" << y*x_0 << "\n";
            }
            File << "\n\n";
        }
    }

private:
    void set_par_ids(size_t sqr_num_track_par, size_t par_num_per_dim)
    {
        size_t num_track_par = sqr_num_track_par*sqr_num_track_par;
        par_ids.reserve(num_track_par);
        size_t x, y, z;
        FTYPE_t s;
        y = par_num_per_dim / 2; // middle of the cube
        s = par_num_per_dim / FTYPE_t(4*(sqr_num_track_par+1)); // quarter of the cube
        for (size_t i=1; i<=sqr_num_track_par;i++)
        {
            z = (size_t)(s*i);
            for (size_t j=1; j<=sqr_num_track_par;j++)
            {
                x = (size_t)(s*j);
                par_ids.push_back(x*par_num_per_dim*par_num_per_dim+y*par_num_per_dim+z);
            }
        }
    }
	
	// VARIABLES
	std::vector<size_t> par_ids;
	std::vector<std::vector<Particle_x<FTYPE_t>>> par_pos;
};

//  ******************************
}// * END OF ANONYMOUS NAMESPACE *
//  ******************************

/**
 * @brief implementation of class App_Var<T>
 * @class Impl
 * 
 */
template <class T> 
class App_Var<T>::Impl
{
public:
    // CONSTRUCTOR
    Impl(const Sim_Param &sim, const std::string& app_short, const std::string& app_long):
        step(0), print_every(sim.out_opt.print_every),
        app_str(app_short), app_long(app_long), z_suffix_const("_" + app_short + "_"), out_dir_app(std_out_dir(app_short + "_run/", sim)),
        track(4, sim.box_opt.par_num_1d),
        a(sim.integ_opt.a_in), a_out(sim.integ_opt.a_out), da(sim.integ_opt.db),
        is_init_pwr_spec_0(false), is_init_vel_pwr_spec_0(false)
    {}

    // ALLOCATE MEMORY
    uint64_t alloc_mesh_vec(App_Var<T>& APP)
    {
        APP.app_field.reserve(3);
        APP.power_aux.reserve(3);
        for(size_t i = 0; i < 3; i++){
            APP.app_field.emplace_back(APP.sim.box_opt.mesh_num);
            APP.power_aux.emplace_back(APP.sim.box_opt.mesh_num_pwr);
        }

        return sizeof(FTYPE_t)*(APP.app_field[0].length*APP.app_field.size()+APP.power_aux[0].length*APP.power_aux.size());
    }

    uint64_t alloc_bin_spec(App_Var<T>& APP)
    {
        size_t bin_num = (size_t)ceil(log10(APP.sim.box_opt.mesh_num_pwr)*APP.sim.out_opt.bins_per_decade);
        APP.pwr_spec_binned.reserve(bin_num);
        APP.pwr_spec_binned_0.reserve(bin_num);
        APP.vel_pwr_spec_binned_0.reserve(bin_num);

        return sizeof(FTYPE_t)*APP.pwr_spec_binned.dim()*APP.pwr_spec_binned.size()*3;
    }

    uint64_t alloc_bin_corr(App_Var<T>& APP)
    {
        size_t bin_num =  (size_t)ceil((APP.sim.other_par.x_corr.upper - APP.sim.other_par.x_corr.lower)/ 10. * APP.sim.out_opt.points_per_10_Mpc);
        APP.corr_func_binned.reserve(bin_num);

        return sizeof(FTYPE_t)*APP.corr_func_binned.dim()*APP.corr_func_binned.size();
    }

    uint64_t alloc_particles(App_Var<T>& APP)
    {
        APP.particles.resize(APP.sim.box_opt.par_num); //< use resize instead of reserve for initialization of particles
        return sizeof(T)*APP.sim.box_opt.par_num;
    }

    void fftw_prep(App_Var<T>& APP)
    {
        const Sim_Param& sim = APP.sim; // get rid of 'APP.sim'

        if (!FFTW_PLAN_OMP_INIT()){
            throw std::runtime_error("Errors during multi-thread initialization");
        }
        FFTW_PLAN_OMP(sim.run_opt.nt_fftw);

        
        
        const size_t N_pot = sim.box_opt.mesh_num;
        const size_t N_pwr = sim.box_opt.mesh_num_pwr;

        APP.p_F = FFTW_PLAN_R2C(N_pot, N_pot, N_pot, APP.app_field[0].real(), APP.app_field[0].complex(), FFTW_ESTIMATE);
        APP.p_B = FFTW_PLAN_C2R(N_pot, N_pot, N_pot, APP.app_field[0].complex(), APP.app_field[0].real(), FFTW_ESTIMATE);
        APP.p_F_pwr = FFTW_PLAN_R2C(N_pwr, N_pwr, N_pwr, APP.power_aux[0].real(), APP.power_aux[0].complex(), FFTW_ESTIMATE);
        APP.p_B_pwr = FFTW_PLAN_C2R(N_pwr, N_pwr, N_pwr, APP.power_aux[0].complex(), APP.power_aux[0].real(), FFTW_ESTIMATE);
    }

    /*********************
    * INITIAL CONDITIONS *
    *********************/

    void set_init_pos(App_Var<T>& APP)
    {   
        set_pert_pos(APP.sim, APP.sim.integ_opt.a_in, APP.particles, APP.app_field);
    }

    void set_init_cond(App_Var<T>& APP)
    {
        /* Generating the right density distribution in k-space */	
        gen_rho_dist_k(APP.sim, APP.app_field[0], APP.p_F);

        /* Print input power spectrum (one realisation), before Zel`dovich push */
        if (print_every && APP.sim.out_opt.print_pwr) print_input_realisation(APP);
        
        /* Computing initial potential in k-space */
        gen_pot_k(APP.app_field[0], APP.power_aux[0]);
        
        /* Computing displacement in k-space */
        gen_displ_k(APP.app_field, APP.power_aux[0]);
        
        /* Computing displacement in q-space */
        BOOST_LOG_TRIVIAL(debug) << "Computing displacement in q-space...";
        fftw_execute_dft_c2r_triple(APP.p_B, APP.app_field);

        /* Set initial positions of particles (with or without velocities)  */
        set_init_pos(APP);
    }

    // PUBLIC PRINTING
    const std::string app_str, app_long, z_suffix_const, out_dir_app;

    void print_sim_name() const
    {
        std::string app_long_upper(app_long);
        std::transform(app_long_upper.begin(), app_long_upper.end(), app_long_upper.begin(), ::toupper);
        BOOST_LOG_TRIVIAL(info) << "Starting: " << app_long_upper;
    }

    void print_end()
    {
        BOOST_LOG_TRIVIAL(info) << app_long << " ended successfully.";
    }

    std::string z_suffix() const
    {
        std::stringstream z_suffix_num;
        z_suffix_num << std::fixed << std::setprecision(2) << z();
        return z_suffix_const + "z" + z_suffix_num.str();
    }
    
    void print_input_realisation(App_Var<T>& APP)
    {
        /* Print input power spectrum (one realisation), before Zel`dovich push */
        pwr_spec_k_init(APP.app_field[0], APP.power_aux[0]);
        gen_pow_spec_binned_init(APP.sim, APP.power_aux[0], APP.app_field[0].length/2, APP.pwr_spec_binned_0);
        pwr_spec_input.init(APP.pwr_spec_binned_0);
        print_pow_spec(APP.pwr_spec_binned_0, out_dir_app, z_suffix_const + "init");
    }

    void print_output(App_Var<T>& APP)
    {
        const Out_Opt& out_opt = APP.sim.out_opt;

        /* Printing positions */
        if (out_opt.print_par_pos) print_position(APP);

        /* Get discrete density from particles */
        if (out_opt.get_rho) get_rho_from_par(APP.particles, APP.power_aux[0], APP.sim);
        
        /* Printing density */
        if (out_opt.print_dens) print_density(APP);

        /* Compute power spectrum and bin it */
        if (out_opt.get_pwr) get_binned_power_spec(APP);

        /* Printing power spectrum */
        if (out_opt.print_pwr) print_power_spec(APP);

        /* Extrapolate power spectrum beyond range of simulation box */
        if (out_opt.get_pk_extrap){
            const Extrap_Pk<FTYPE_t, 2> P_k(APP.pwr_spec_binned, APP.sim);
            /* Print extrapolated power spectrum */
            if (out_opt.print_extrap_pwr) print_extrap_pwr(APP, P_k);
            /* Printing correlation function */
            if (out_opt.print_corr) print_corr(APP, P_k);
        }

        /* Velocity power spectrum */
        if (out_opt.print_vel_pwr && get_vel_from_par(APP.particles, APP.power_aux, APP.sim)) print_vel_pwr(APP);
    }

    // CREATE WORKING DIRECTORY STRUCTURE
    void create_work_dir(const Out_Opt& out_opt)
    {
        create_dir(out_dir_app + "log/");
        if (print_every)
        {
            if (out_opt.print_corr) create_dir(out_dir_app + "corr_func/");
            if (out_opt.print_par_pos) create_dir(out_dir_app + "par_cut/");
            if (out_opt.print_pwr) 
            {
                create_dir(out_dir_app + "pwr_diff/");
                create_dir(out_dir_app + "pwr_spec/");
            }
            if (out_opt.print_dens)
            {
                create_dir(out_dir_app + "rho_bin/");
                create_dir(out_dir_app + "rho_map/");
            }
            
            if (out_opt.print_vel_pwr)
            {
                create_dir(out_dir_app + "vel_pwr_diff/");
                create_dir(out_dir_app + "vel_pwr_spec/");
            }
        }
    }

    // INTEGRATION
    FTYPE_t a, a_out, da;

    void integration(App_Var<T>& APP)
    {
        print_init(APP); // WARNING: power_aux[0] is modified
        while(integrate())
        {
            BOOST_LOG_TRIVIAL(info) << "Starting computing step with z = " << z() << " (a = " << a << ")";
            APP.upd_pos();
            track.update_track_par(APP.particles);
            if (printing()) APP.print_output();
            upd_time();
        }
        print_info(APP.sim);
    }

private:
    // PRIVATE PRINTING
    unsigned int print_every, step;
    Tracking track;
    Interp_obj pwr_spec_input;

    bool printing() const
    {
        return print_every ? ((step % print_every) == 0) or (a == a_out) : false ;
    }

    void print_info(const Sim_Param& sim) const
    {
        sim.print_info(out_dir_app, app_str);
    }

    void print_position(const App_Var<T>& APP) const
    {/* Printing positions */
        print_par_pos_cut_small(APP.particles, APP.sim, out_dir_app, z_suffix());
        track.print_track_par(APP.sim, out_dir_app, z_suffix());
    }

    void print_density(App_Var<T>& APP) const
    {/* Printing density */
        gen_dens_binned(APP.power_aux[0], APP.dens_binned, APP.sim);    
        print_rho_map(APP.power_aux[0], APP.sim, out_dir_app, z_suffix());
        print_dens_bin(APP.dens_binned, out_dir_app, z_suffix());
    }

    void get_binned_power_spec(App_Var<T>& APP) const
    {/* Compute power spectrum and bin it */
        fftw_execute_dft_r2c(APP.p_F_pwr, APP.power_aux[0]);
        pwr_spec_k(APP.power_aux[0], APP.power_aux[0]);
        gen_pow_spec_binned(APP.sim, APP.power_aux[0], APP.pwr_spec_binned);
    }

    void print_power_spec(App_Var<T>& APP)
    {/* Printing power spectrum */
        print_pow_spec(APP.pwr_spec_binned, out_dir_app, "_par" + z_suffix());
        if (!is_init_pwr_spec_0){
            APP.pwr_spec_binned_0 = APP.pwr_spec_binned;
            D_init = growth_factor(a, APP.sim.cosmo);
            is_init_pwr_spec_0 = true;
        }
        FTYPE_t D_now = growth_factor(a, APP.sim.cosmo);
        print_pow_spec_diff(APP.pwr_spec_binned, APP.pwr_spec_binned_0, D_now / D_init, out_dir_app, "_par" + z_suffix());
        print_pow_spec_diff(APP.pwr_spec_binned, pwr_spec_input, D_now, out_dir_app, "_input" + z_suffix());
        print_pow_spec_diff(APP.pwr_spec_binned, APP.pwr_spec_binned_0, pwr_spec_input, D_now, D_init,
                            out_dir_app, "_hybrid" + z_suffix());
    }


    void print_extrap_pwr(App_Var<T>& APP, const Extrap_Pk<FTYPE_t, 2>& P_k) const
    {/* Print extrapolated power spectrum */
        gen_pow_spec_binned_from_extrap(APP.sim, P_k, APP.pwr_spec_binned);
        print_pow_spec(APP.pwr_spec_binned, out_dir_app, "_extrap" + z_suffix());
    }

    void print_corr(App_Var<T>& APP, const Extrap_Pk<FTYPE_t, 2>& P_k) const
    {/* Printing correlation function */
        gen_corr_func_binned_gsl_qawf(APP.sim, P_k, APP.corr_func_binned);
        print_corr_func(APP.corr_func_binned, out_dir_app, "_gsl_qawf_par" + z_suffix());
        gen_corr_func_binned_gsl_qawf_lin(APP.sim, a, APP.corr_func_binned);
        print_corr_func(APP.corr_func_binned, out_dir_app, "_gsl_qawf_lin" + z_suffix());
    }

    void print_vel_pwr(App_Var<T>& APP)
    {/* Print velocity power spectrum */
        fftw_execute_dft_r2c_triple(APP.p_F_pwr, APP.power_aux);
        vel_pwr_spec_k(APP.power_aux, APP.power_aux[0]);
        gen_pow_spec_binned(APP.sim, APP.power_aux[0], APP.pwr_spec_binned);
        print_vel_pow_spec(APP.pwr_spec_binned, out_dir_app, z_suffix());
        if (!is_init_vel_pwr_spec_0){
            APP.vel_pwr_spec_binned_0 = APP.pwr_spec_binned;
            is_init_vel_pwr_spec_0 = true;
            dDda_init = growth_change(a, APP.sim.cosmo);
        }
        print_vel_pow_spec_diff(APP.pwr_spec_binned, APP.vel_pwr_spec_binned_0, growth_change(a, APP.sim.cosmo) / dDda_init, out_dir_app, z_suffix());
    }

    void print_init(App_Var<T>& APP)
    {
        /* Setting initial (binned) power spectrum, WARNING: power_aux[0] is modified */
        track.update_track_par(APP.particles);
        if (print_every) APP.print_output();
        upd_time();
    }

    // INTEGRATION
    FTYPE_t D_init, dDda_init;
    bool is_init_pwr_spec_0, is_init_vel_pwr_spec_0;

    bool integrate() const
    {
        return (a <= a_out) && (da > 0);
    }

    void upd_time()
    {
        step++;
        if ((a_out - a) < da) da = a_out - a;
        a += da;
    }


    FTYPE_t z() const
    {
        return 1/a - 1;
    }    
};

template <class T> 
App_Var<T>::App_Var(const Sim_Param &sim, const std::string& app_short, const std::string& app_long):
	m_impl(new Impl(sim, app_short, app_long)), sim(sim), dens_binned(500)
{
    // EFFICIENTLY ALLOCATE MEMORY
    memory_alloc = m_impl->alloc_mesh_vec(*this); // app_field, power_aux
    memory_alloc += m_impl->alloc_bin_spec(*this); // pwr_spec_binned, pwr_spec_binned_0, vel_pwr_spec_binned_0
    memory_alloc += m_impl->alloc_bin_corr(*this); // corr_func_binned
    memory_alloc += m_impl->alloc_particles(*this); // particles
    
    // CREAT SUBDIR STRUCTURE
    m_impl->create_work_dir(sim.out_opt);

	// FFTW PREPARATION
    m_impl->fftw_prep(*this); // fftw omp, fftw plans
}

template <class T> 
App_Var<T>::~App_Var()
{	// FFTW CLEANUP
	FFTW_DEST_PLAN(p_F);
    FFTW_DEST_PLAN(p_B);
    FFTW_DEST_PLAN(p_F_pwr);
	FFTW_DEST_PLAN(p_B_pwr);
	FFTW_PLAN_OMP_CLEAN();
}

template <class T> 
void App_Var<T>::run_simulation()
{
    // print simulation name
    m_impl->print_sim_name();

    // print memory usage
    print_mem(memory_alloc);

    // set initial conditions
    m_impl->set_init_cond(*this);

    // CIC correction of potential (if not overriden)
    pot_corr(app_field, power_aux[0]);

    // integration
    m_impl->integration(*this);

    // end of simulation
    m_impl->print_end();
}

template <class T> 
void App_Var<T>::update_cosmo(Cosmo_Param& cosmo)
{
    cosmo.truncated_pk = false;
}

template <class T>
void App_Var<T>::pot_corr(std::vector<Mesh>& vel_field, Mesh& pot_k)
{
    /* Computing displacement in k-space with CIC opt */
    gen_displ_k_cic(vel_field, pot_k);

    /* Computing force in q-space */
    BOOST_LOG_TRIVIAL(debug) << "Computing force in q-space...";
    fftw_execute_dft_c2r_triple(p_B, vel_field);
}

template <class T> 
void App_Var<T>::print_output()
{
    m_impl->print_output(*this);
}

template <class T> 
FTYPE_t App_Var<T>::a()
{
    return m_impl->a;
}

template <class T> 
FTYPE_t App_Var<T>::a_half()
{
    return a() - da()/2.;
}

template <class T> 
FTYPE_t App_Var<T>::da()
{
    return m_impl->da;
}

template <class T> 
std::string App_Var<T>::get_out_dir() const
{
    return m_impl->out_dir_app;
}

template <class T> 
std::string App_Var<T>::get_z_suffix() const
{
    return m_impl->z_suffix();
}

template class App_Var<Particle_x<FTYPE_t>>;
template class App_Var<Particle_v<FTYPE_t>>;