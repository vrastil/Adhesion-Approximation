#include "approximations.hpp"
#include "core_app.h"
#include "core_mesh.h"
#include "core_out.h"
#include "chameleon.hpp"

using namespace std;

namespace {

/*********************
* INITIAL CONDITIONS *
*********************/

void init_cond(App_Var<Particle_x<FTYPE_t>>& APP)
{
    printf("\nSetting initial positions of particles...\n");
    set_pert_pos(APP.sim, APP.sim.integ_opt.b_in, APP.particles, APP.app_field);
}

void init_cond(App_Var<Particle_v<FTYPE_t>>& APP)
{
    printf("\nSetting initial positions and velocitis of particles...\n");
	set_pert_pos_w_vel(APP.sim, APP.sim.integ_opt.b_in, APP.particles, APP.app_field);
}

template<class T>
void init_pot_w_s2(T& APP)
{
    /* Computing displacement in k-space with S2 shaped particles */
	gen_displ_k_S2(APP.app_field, APP.power_aux[0], APP.sim.app_opt.a);
    
    /* Computing force in q-space */
    printf("Computing force in q-space...\n");
    fftw_execute_dft_c2r_triple(APP.p_B, APP.app_field);
}

template<class T>
void init_pot_w_cic(T& APP)
{
    /* Computing displacement in k-space with CIC opt */
    gen_displ_k_cic(APP.app_field, APP.power_aux[0]);

    /* Computing force in q-space */
    printf("Computing force in q-space...\n");
    fftw_execute_dft_c2r_triple(APP.p_B, APP.app_field);
}

template<class T>
void print_init(T& APP)
{
    /* Setting initial (binned) power spectrum, WARNING: power_aux[0] is modified */
    APP.track.update_track_par(APP.particles);
    if (APP.print_every) APP.print_output();
    APP.upd_time();
}

template<class T>
void print_input_realisation(T& APP)
{
    /* Print input power spectrum (one realisation), before Zel`dovich push */
    pwr_spec_k_init(APP.app_field[0], APP.power_aux[0]);
    gen_pow_spec_binned_init(APP.sim, APP.power_aux[0], APP.app_field[0].length/2, APP.pwr_spec_binned_0);
    APP.pwr_spec_input.init(APP.pwr_spec_binned_0);
    print_pow_spec(APP.pwr_spec_binned_0, APP.out_dir_app, APP.z_suffix_const + "init");
}

/***************************************
* STANDARD PREPARATION FOR INTEGRATIOM *
***************************************/

template<class T>
void standard_preparation(T& APP, function<void()> save_rho_k = [](){;})
{
    /* Generating the right density distribution in k-space */	
    gen_rho_dist_k(APP.sim, APP.app_field[0], APP.p_F);

    /* Save initial density field */
    save_rho_k();

    /* Print input power spectrum (one realisation), before Zel`dovich push */
    if (APP.print_every) print_input_realisation(APP);
    
	/* Computing initial potential in k-space */
	gen_pot_k(APP.app_field[0], APP.power_aux[0]);
	
	/* Computing displacement in k-space */
	gen_displ_k(APP.app_field, APP.power_aux[0]);
    
    /* Computing displacement in q-space */
    printf("Computing displacement in q-space...\n");
    fftw_execute_dft_c2r_triple(APP.p_B, APP.app_field);
}

/* preparation for chameleon gravity */
void standard_preparation(App_Var_chi& APP)
{
    auto save_rho_k = [&](){ APP.save_init_drho_k(APP.app_field[0], APP.app_field[1]); };
    standard_preparation<App_Var_chi>(APP, save_rho_k);
}

/**************
* INTEGRATION *
**************/

template<class T>
void integration(T& APP, function<void(T&)> upd_pos)
{
    print_init(APP); // WARNING: power_aux[0] is modified
	while(APP.integrate())
	{
		printf("\nStarting computing step with z = %.2f (b = %.3f)\n", APP.z(), APP.b);
		upd_pos(APP);
        APP.track.update_track_par(APP.particles);
		if (APP.printing()) APP.print_output();
		APP.upd_time();
    }
    APP.print_info();
}

void stream_step(const Sim_Param &sim, const FTYPE_t da, vector<Particle_v<FTYPE_t>>& particles)
{
    const unsigned Np = sim.box_opt.par_num;
    #pragma omp parallel for
	for (unsigned i = 0; i < Np; i++)
	{
        particles[i].position += particles[i].velocity*da;
    }
}

void stream_kick_stream(const Sim_Param &sim, const FTYPE_t da, vector<Particle_v<FTYPE_t>>& particles, function<void()> kick_step)
{// general Leapfrog method: Stream-Kick-Stream & ensure periodicity
    stream_step(sim, da/2, particles);
    kick_step();
    stream_step(sim, da/2, particles);
    get_per(particles, sim.box_opt.par_num, sim.box_opt.mesh_num);
}

/**************
 * KICK STEPS *
 **************/

void kick_step_no_momentum(const Sim_Param &sim, const FTYPE_t a, vector<Particle_v<FTYPE_t>>& particles, const vector< Mesh> &vel_field)
{
    // no memory of previus velocity, 1st order ODE
    const unsigned Np = sim.box_opt.par_num;
    Vec_3D<FTYPE_t> vel;
    const FTYPE_t dDda = growth_change(a, sim.cosmo); // dD / da
    
    #pragma omp parallel for private(vel)
    for (unsigned i = 0; i < Np; i++)
	{
        vel.fill(0.);
        assign_from(vel_field, particles[i].position, vel);
        particles[i].velocity = vel*dDda;
    }
}

void kick_step_w_momentum(const Sim_Param &sim, const FTYPE_t a, const FTYPE_t da, vector<Particle_v<FTYPE_t>>& particles, const vector< Mesh> &force_field)
{
    // classical 2nd order ODE
    const unsigned Np = sim.box_opt.par_num;
    Vec_3D<FTYPE_t> force;
    const FTYPE_t D = growth_factor(a, sim.cosmo);
    const FTYPE_t OL = sim.cosmo.Omega_L()*pow(a,3);
    const FTYPE_t Om = sim.cosmo.Omega_m;
    // -3/2a represents usual EOM, the rest are LCDM corrections
    const FTYPE_t f1 = 3/(2*a)*(Om+2*OL)/(Om+OL);
    const FTYPE_t f2 = 3/(2*a)*Om/(Om+OL)*D/a;
    
    #pragma omp parallel for private(force)
    for (unsigned i = 0; i < Np; i++)
	{
        force.fill(0.);
        assign_from(force_field, particles[i].position, force);
        force = force*f2 - particles[i].velocity*f1;		
        particles[i].velocity += force*da;
    }
}

void kick_step_w_pp(const Sim_Param &sim, const FTYPE_t a, const FTYPE_t da, vector<Particle_v<FTYPE_t>>& particles, const vector< Mesh> &force_field,
                    LinkedList& linked_list, Interp_obj& fs_interp)
{    // 2nd order ODE with long & short range potential
    const unsigned Np = sim.box_opt.par_num;
    Vec_3D<FTYPE_t> force;
    const FTYPE_t D = growth_factor(a, sim.cosmo);
    const FTYPE_t OL = sim.cosmo.Omega_L()*pow(a,3);
    const FTYPE_t Om = sim.cosmo.Omega_m;
    // -3/2a represents usual EOM, the rest are LCDM corrections
    const FTYPE_t f1 = 3/(2*a)*(Om+2*OL)/(Om+OL);
    const FTYPE_t f2 = 3/(2*a)*Om/(Om+OL)*D/a;
    
    printf("Creating linked list...\n");
	linked_list.get_linked_list(particles);

    cout << "Computing short and long range parts of the potential...\n";
    #pragma omp parallel for private(force)
    for (unsigned i = 0; i < Np; i++)
	{
        force.fill(0.);
        assign_from(force_field, particles[i].position, force); // long-range force
        force_short(sim, D, linked_list, particles, particles[i].position, force, fs_interp); // short range force

        force = force*f2 - particles[i].velocity*f1;		
        particles[i].velocity += force*da;
    }
}

/***********************************
* ADHESION APPROXIMATION FUNCTIONS *
***********************************/

const FTYPE_t ACC = 1e-10;
const FTYPE_t log_acc = log(ACC);

void gen_init_expot(const Mesh& potential, Mesh& expotential, FTYPE_t nu)
{
	printf("Storing initial expotenital in q-space...\n");
    // store exponent only
    // *expotential = potential; !!! <- do not use this in case potential and expotential are meshes of different size
    #pragma omp parallel for
    for (unsigned i = 0; i < expotential.length; i++) expotential[i] = -potential[i] / (2*nu);
}

FTYPE_t get_summation(const vector<FTYPE_t>& exp_aux)
{
    FTYPE_t max_exp = *max_element(exp_aux.begin(), exp_aux.end());
    FTYPE_t sum = 0;
    for(auto const& a_exp: exp_aux) {
        if ((a_exp - max_exp) > log_acc) sum+= exp(a_exp - max_exp);
    }
    return max_exp + log(sum);
}

void convolution_y1(Mesh& potential, const vector<FTYPE_t>& gaussian, const Mesh& expotential_0){
	// multi-thread index is y3
    // compute f1 (x1, y2, y3)

    const int N = potential.N;
    vector<FTYPE_t> exp_aux;
    
	#pragma omp parallel for private(exp_aux)
	for (int x1 = 0; x1 < N; x1++){
		for (int y2 = 0; y2 < N; y2++){
			for (int y3 = 0; y3 < N; y3++){
                exp_aux.reserve(N);
                // fill in exponents
                for (int y1 = 0; y1 < N; y1++){
                    exp_aux.push_back(expotential_0(y1, y2, y3)+gaussian[abs(x1-y1)]);
				}
				potential(x1, y2, y3) = get_summation(exp_aux); // potential is now f1
                exp_aux.clear();
			}
		}
	}
}

void convolution_y2(Mesh& potential, const vector<FTYPE_t>& gaussian){
    // compute f2 (x1, x2, y3)

    const int N = potential.N;
	vector<FTYPE_t> sum_aux;
	vector<FTYPE_t> exp_aux;

	#pragma omp parallel for private(sum_aux, exp_aux)
	for (int x1 = 0; x1 < N; x1++){
		for (int y3 = 0; y3 < N; y3++){
            sum_aux.reserve(N);
			for (int x2 = 0; x2 < N; x2++){
                exp_aux.reserve(N);
				// fill in exponents
                for (int y2 = 0; y2 < N; y2++){
                    exp_aux.push_back(potential(x1, y2, y3) + gaussian[abs(x2-y2)]);
				}
				sum_aux.push_back(get_summation(exp_aux));
                exp_aux.clear();
			}

			for (int x2 = 0; x2 < N; x2++){
				potential(x1, x2, y3) = sum_aux[x2]; // potential is now f2
			}
			sum_aux.clear();
		}
	}
}

void convolution_y3(Mesh& potential, const vector<FTYPE_t>& gaussian){
    // compute f3 (x1, x2, x3) == expotential(x, b)

    const int N = potential.N;
	vector<FTYPE_t> sum_aux;
    vector<FTYPE_t> exp_aux;

	#pragma omp parallel for private(sum_aux, exp_aux)
	for (int x1 = 0; x1 < N; x1++){
		for (int x2 = 0; x2 < N; x2++){
            sum_aux.reserve(N);
			for (int x3 = 0; x3 < N; x3++){
                exp_aux.reserve(N);
				// fill in exponents
                for (int y3 = 0; y3 < N; y3++){
                    exp_aux.push_back(potential(x1, x2, y3) + gaussian[abs(x3-y3)]);
				}
				sum_aux.push_back(get_summation(exp_aux));
                exp_aux.clear();
			}
			for (int x3 = 0; x3 < N; x3++){
				potential(x1, x2, x3) = sum_aux[x3]; // potential is now f3
			}
			sum_aux.clear();
		}
	}
}

void gen_expot(Mesh& potential,  const Mesh& expotential_0, FTYPE_t nu, FTYPE_t b)
{
	/* Computing convolution using direct sum */
	printf("Computing expotential in q-space...\n");
	/*
	f(x1, x2, x3) = \int dy^3 { g(y1, y2, y3) * h(x1 - y1) * h(x2 - y2) * h(x3 - y3)}
	..
	f1 (x1, y2, y3) = \int dy1 { g  (y1, y2, y3) * h(x1 - y1)}	:: N^3 sums of length N
	f2 (x1, x2, y3) = \int dy2 { f1 (x1, y2, y3) * h(x2 - y2)}	:: N^3 sums of length N
	f3 (x1, x2, x3) = \int dy3 { f2 (x1, x2, y3) * h(x3 - y3)}	:: N^3 sums of length N
	*/

	// store values of exponential - every convolution uses the same exp(-r^2/4bv)
	vector<FTYPE_t> gaussian(expotential_0.N);

	#pragma omp parallel for
	for (unsigned i = 0; i < expotential_0.N; i++){
		gaussian[i]=-i*i/(4*b*nu);
	}

	convolution_y1(potential, gaussian, expotential_0);
	convolution_y2(potential, gaussian);
	convolution_y3(potential, gaussian);
}


void aa_convolution(App_Var_AA& APP)
{
    printf("Computing potential...\n");	
    gen_expot(APP.app_field[0], APP.expotential, APP.sim.app_opt.nu, APP.b_half());
	// gen_expot(APP.app_field[0], APP.expotential, sim.app_opt.nu, APP.b);
    APP.app_field[0] *= -2*APP.sim.app_opt.nu;
				
	printf("Computing velocity field via FFT...\n");
	fftw_execute_dft_r2c(APP.p_F, APP.app_field[0]);
	gen_displ_k(APP.app_field, APP.app_field[0]);
	fftw_execute_dft_c2r_triple(APP.p_B, APP.app_field);
}

template<class T>
void init_adhesion(T& APP)
{
    /* Computing initial expotential */
	fftw_execute_dft_c2r(APP.p_B, APP.power_aux[0]);
    gen_init_expot(APP.power_aux[0], APP.expotential, APP.sim.app_opt.nu);
}

/*****************
* APPROXIMATIONS *
*****************/

template<class T>
void general_sim_flow(const Sim_Param &sim, const string& app_short, const string& app_long,
                             function<void(T&)> upd_pos, function<void(T&)> pot_corr = [](T&){;})
{
    // print simulation name
    string app_long_upper(app_long);
    transform(app_long_upper.begin(), app_long_upper.end(), app_long_upper.begin(), ::toupper);
    cout << "\n" << std::string(app_long_upper.length(), '*') << "\n"
         << app_long_upper
         << "\n" << std::string(app_long_upper.length(), '*') << "\n";

    // initialize approximation & print memory
    T APP(sim, app_short);
    APP.print_mem();

    // standard preparation & init conditions & possible corrections
    standard_preparation(APP);
    init_cond(APP); //< with or without velocities
    pot_corr(APP);

    // integration
    integration(APP, upd_pos);

    cout << app_long << " ended successfully.\n";
}
}// end of anonymous namespace

void zel_app(const Sim_Param &sim)
{
    typedef App_Var<Particle_v<FTYPE_t>> APP_t;
    auto upd_pos = [](APP_t& APP){
        set_pert_pos_w_vel(APP.sim, APP.b, APP.particles, APP.app_field); //< ZA with velocitites
    };
    general_sim_flow<APP_t>(sim, "ZA", "Zel`dovich approximation", upd_pos);
}

void frozen_flow(const Sim_Param &sim)
{
    // frozen-flow approximation type
    typedef App_Var<Particle_v<FTYPE_t>> APP_t;

    // Leapfrog method for frozen-flow
    auto upd_pos = [](APP_t& APP){
        auto kick_step = [&](){ kick_step_no_momentum(APP.sim, APP.b-APP.db/2, APP.particles, APP.app_field); };
        stream_kick_stream(APP.sim, APP.db, APP.particles, kick_step);
    };

    // force interpolation corrections
    auto pot_corr = [](APP_t& APP){ init_pot_w_cic(APP); };

    // start general work-flow of simulations
    general_sim_flow<APP_t>(sim, "FF", "Frozen-flow approximation", upd_pos, pot_corr);
}

void frozen_potential(const Sim_Param &sim)
{
    // frozen-potential approximation type
    typedef App_Var<Particle_v<FTYPE_t>> APP_t;

    // Leapfrog method for frozen-potential
    auto upd_pos = [](APP_t& APP){
        auto kick_step = [&](){ kick_step_w_momentum(APP.sim, APP.b-APP.db/2, APP.db, APP.particles, APP.app_field); };
        stream_kick_stream(APP.sim, APP.db, APP.particles, kick_step);
    };

    // force interpolation corrections
    auto pot_corr = [](APP_t& APP){ init_pot_w_cic(APP); };

    // start general work-flow of simulations
    general_sim_flow<APP_t>(sim, "FP", "Frozen-potential approximation", upd_pos, pot_corr);
}

void mod_frozen_potential(const Sim_Param &sim)
{
    // modified frozen-potential approximation type
    typedef App_Var_FP_mod APP_t;

    // Leapfrog method for modified frozen-potential
    auto upd_pos = [](APP_t& APP){
        auto kick_step = [&](){ kick_step_w_pp(APP.sim, APP.b-APP.db/2, APP.db, APP.particles, APP.app_field, APP.linked_list, APP.fs_interp); };
        stream_kick_stream(APP.sim, APP.db, APP.particles, kick_step);
    };

    // force interpolation corrections, long range potential for S2-shaped particles
    auto pot_corr = [](APP_t& APP){ init_pot_w_s2(APP); };

    // start general work-flow of simulations
    general_sim_flow<APP_t>(sim, "FP_pp", "Modified Frozen-potential approximation", upd_pos, pot_corr);
}

void adhesion_approximation(const Sim_Param &sim)
{
    // adhesion approximation type
    typedef App_Var_AA APP_t;

    // Leapfrog method for adhesion
    auto upd_pos = [](APP_t& APP){
        aa_convolution(APP);
        auto kick_step = [&](){ kick_step_w_momentum(APP.sim, APP.b-APP.db/2, APP.db, APP.particles, APP.app_field); };
        stream_kick_stream(APP.sim, APP.db, APP.particles, kick_step);
    };

    // initialize potentials for adhesion
    auto pot_corr = [](APP_t& APP){ init_adhesion(APP); };

    // start general work-flow of simulations
    general_sim_flow<APP_t>(sim, "AA", "Adhesion approximation", upd_pos, pot_corr);
}

/**********************
 * MODIFIED GRAVITIES *
 *********************/

void chameleon_gravity(const Sim_Param &sim)
{
    // chameleon gravity approximation type
    typedef App_Var_chi APP_t;

    // Leapfrog method for chameleon gravity (frozen-potential)
    auto upd_pos = [](APP_t& APP){
        auto kick_step = [&]()
        {
            APP.solve(APP.b-APP.db/2);   
            kick_step_w_momentum(APP.sim, APP.b-APP.db/2, APP.db, APP.particles, APP.app_field);
        };
        stream_kick_stream(APP.sim, APP.db, APP.particles, kick_step);
    };

    // force interpolation corrections
    auto pot_corr = [](APP_t& APP){ init_pot_w_cic(APP); };

    // start general work-flow of simulations
    general_sim_flow<APP_t>(sim, "CHI", "Chameleon gravity", upd_pos, pot_corr);
}