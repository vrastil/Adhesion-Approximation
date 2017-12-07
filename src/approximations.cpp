
#include "stdafx.h"
#include "core.h"

#include "core_app.h"
#include "core_mesh.h"
#include "core_out.h"

using namespace std;

static void gen_init_expot(const Mesh& potential, Mesh* expotential, double nu);
static void gen_expot(Mesh* potential,  const Mesh& expotential, double nu, double b);
static void aa_convolution(App_Var_AA* APP);

/*********************
* INITIAL CONDITIONS *
*********************/

template<class T>
void init_cond_no_vel(T& APP)
{
    printf("\nSetting initial positions of particles...\n");
    set_pert_pos(APP.sim, APP.integ_opt.b_in, APP.particles, APP.app_field);
}

template<class T>
void init_cond_w_vel(T& APP)
{
    printf("\nSetting initial positions and velocitis of particles...\n");
	set_pert_pos_w_vel(APP.sim, APP.sim.integ_opt.b_in, APP.particles, APP.app_field);
}

template<class T>
void init_pot_w_s2(T& APP)
{
    /* Computing displacement in k-space with S2 shaped particles */
	gen_displ_k_S2(&APP.app_field, APP.power_aux[0], APP.sim.app_opt.a);
    
    /* Computing force in q-space */
    printf("Computing force in q-space...\n");
    fftw_execute_dft_c2r_triple(APP.p_B, APP.app_field);
}

template<class T>
void init_pot_w_cic(T& APP)
{
    /* Computing displacement in k-space with CIC opt */
    gen_displ_k_cic(&APP.app_field, APP.power_aux[0]);

    /* Computing force in q-space */
    printf("Computing force in q-space...\n");
    fftw_execute_dft_c2r_triple(APP.p_B, APP.app_field);
}

template<class T>
void init_adhesion(T& APP)
{
    /* Computing initial expotential */
	fftw_execute_dft_c2r(APP.p_B, APP.power_aux[0]);
    gen_init_expot(APP.power_aux[0], &APP.expotential, APP.sim.app_opt.nu);
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
    pwr_spec_k_init(APP.app_field[0], &APP.power_aux[0]);
    gen_pow_spec_binned_init(APP.sim, APP.power_aux[0], APP.app_field[0].length/2, &APP.pwr_spec_binned_0);
    APP.pwr_spec_input.init(APP.pwr_spec_binned_0);
    print_pow_spec(APP.pwr_spec_binned_0, APP.out_dir_app, APP.z_suffix_const + "init");
}

/***************************************
* STANDARD PREPARATION FOR INTEGRATIOM *
***************************************/

template<class T>
void standard_preparation(T& APP)
{
    /* Generating the right density distribution in k-space */	
    gen_rho_dist_k(APP.sim, &APP.app_field[0], APP.p_F);

    /* Print input power spectrum (one realisation), before Zel`dovich push */
    print_input_realisation(APP);
    
	/* Computing initial potential in k-space */
	gen_pot_k(APP.app_field[0], &APP.power_aux[0]);
	
	/* Computing displacement in k-space */
	gen_displ_k(&APP.app_field, APP.power_aux[0]);
    
    /* Computing displacement in q-space */
    printf("Computing displacement in q-space...\n");
    fftw_execute_dft_c2r_triple(APP.p_B, APP.app_field);
}

/**************
* INTEGRATION *
**************/

template<class T>
void integration(T& APP, function<void()> upd_pos)
{
    print_init(APP); // WARNING: power_aux[0] is modified
	while(APP.integrate())
	{
		printf("\nStarting computing step with z = %.2f (b = %.3f)\n", APP.z(), APP.b);
		upd_pos();
        APP.track.update_track_par(APP.particles);
		if (APP.printing()) APP.print_output();
		APP.upd_time();
    }
    APP.print_info();
}

/*****************
* APPROXIMATIONS *
*****************/

void zel_app(const Sim_Param &sim)
{
    cout << "\n"
	"************************\n"
	"ZEL`DOVICH APPROXIMATION\n"
    "************************\n";
    App_Var<Particle_v> APP(sim, "ZA");
    APP.print_mem();
    standard_preparation(APP);
    init_cond_w_vel(APP); //< with velocities
    auto upd_pos = [&](){
        set_pert_pos_w_vel(APP.sim, APP.b, APP.particles, APP.app_field); //< ZA with velocitites
    };
    integration(APP, upd_pos);
	printf("Zel`dovich approximation ended successfully.\n");
}

void frozen_flow(const Sim_Param &sim)
{
	cout << "\n"
	"*************************\n"
	"FROZEN-FLOW APPROXIMATION\n"
	"*************************\n";
    App_Var<Particle_v> APP(sim, "FF");
    APP.print_mem();
    standard_preparation(APP);
    init_cond_w_vel(APP); //< with velocities
    init_pot_w_cic(APP); //< force interpolation corrections
    auto upd_pos = [&](){
        upd_pos_first_order(APP.sim, APP.db, APP.b, APP.particles, APP.app_field); //< FF specific
    };
    integration(APP, upd_pos);
    printf("Frozen-flow approximation ended successfully.\n");
}

void frozen_potential(const Sim_Param &sim)
{
	cout << "\n"
	"******************************\n"
	"FROZEN-POTENTIAL APPROXIMATION\n"
    "******************************\n";
    App_Var<Particle_v> APP(sim, "FP");
    APP.print_mem();
    standard_preparation(APP);
    init_cond_w_vel(APP); //< with velocities
    init_pot_w_cic(APP); //< force interpolation corrections
    auto upd_pos = [&](){
        upd_pos_second_order(APP.sim, APP.db, APP.b, APP.particles, APP.app_field); //< FP specific
    };
    integration(APP, upd_pos);
    printf("Frozen-potential approximation ended successfully.\n");
}

void mod_frozen_potential(const Sim_Param &sim)
{
	cout << "\n"
	"***************************************\n"
	"MODIFIED FROZEN-POTENTIAL APPROXIMATION\n"
	"***************************************\n";
    
    App_Var_FP_mod APP(sim, "FP_pp");
    APP.print_mem();
    standard_preparation(APP);
    init_cond_w_vel(APP); //< with velocities
    init_pot_w_s2(APP); //< force interpolation corrections, long range potential for S2-shaped particles
    auto upd_pos = [&](){
        upd_pos_second_order_w_pp(APP.sim, APP.db, APP.b, APP.particles, APP.app_field, &APP.linked_list, &APP.fs_interp); //< FP_pp specific
    };
    integration(APP, upd_pos);
    printf("Modified Frozen-potential approximation ended successfully.\n");
}

void adhesion_approximation(const Sim_Param &sim)
{
	cout << "\n"
	"**********************\n"
	"ADHESION APPROXIMATION\n"
    "**********************\n";

    App_Var_AA APP(sim, "AA");
    APP.print_mem();
    standard_preparation(APP);
    init_cond_w_vel(APP); //< with velocities
    init_adhesion(APP); //< AA specific
    auto upd_pos = [&](){
        aa_convolution(&APP); //< AA specific
        upd_pos_first_order(APP.sim, APP.db, APP.b, APP.particles, APP.app_field); //< AA specific
    };
    integration(APP, upd_pos);
	printf("Adhesion approximation ended successfully.\n");
}

/***********************************
* ADHESION APPROXIMATION FUNCTIONS *
***********************************/

const double ACC = 1e-10;
const double log_acc = log(ACC);

static void aa_convolution(App_Var_AA* APP)
{
    printf("Computing potential...\n");	
    gen_expot(&APP->app_field[0], APP->expotential, APP->sim.app_opt.nu, APP->b_half());
	// gen_expot(&APP->app_field[0], APP->expotential, sim.app_opt.nu, APP->b);
    APP->app_field[0] *= -2*APP->sim.app_opt.nu;
				
	printf("Computing velocity field via FFT...\n");
	fftw_execute_dft_r2c(APP->p_F, APP->app_field[0]);
	gen_displ_k(&APP->app_field, APP->app_field[0]);
	fftw_execute_dft_c2r_triple(APP->p_B, APP->app_field);
}

static void gen_init_expot(const Mesh& potential, Mesh* expotential, double nu)
{
	printf("Storing initial expotenital in q-space...\n");
    // store exponent only
    // *expotential = potential; !!! <- do not use this in case potential and expotential are meshes of different size
    #pragma omp parallel for
    for (unsigned i = 0; i < expotential->length; i++) (*expotential)[i] = -potential[i] / (2*nu);
}

static double get_summation(const vector<double>& exp_aux)
{
    double max_exp = *max_element(exp_aux.begin(), exp_aux.end());
    double sum = 0;
    for(auto const& a_exp: exp_aux) {
        if ((a_exp - max_exp) > log_acc) sum+= exp(a_exp - max_exp);
    }
    return max_exp + log(sum);
}

static void convolution_y1(Mesh* potential, const vector<double>& gaussian, const Mesh& expotential_0){
	// multi-thread index is y3
    // compute f1 (x1, y2, y3)

    const int N = potential->N;
    vector<double> exp_aux;
    
	#pragma omp parallel for private(exp_aux)
	for (int x1 = 0; x1 < N; x1++){
		for (int y2 = 0; y2 < N; y2++){
			for (int y3 = 0; y3 < N; y3++){
                exp_aux.reserve(N);
                // fill in exponents
                for (int y1 = 0; y1 < N; y1++){
                    exp_aux.push_back(expotential_0(y1, y2, y3)+gaussian[abs(x1-y1)]);
				}
				(*potential)(x1, y2, y3) = get_summation(exp_aux); // potential is now f1
                exp_aux.clear();
			}
		}
	}
}

static void convolution_y2(Mesh* potential, const vector<double>& gaussian){
    // compute f2 (x1, x2, y3)

    const int N = potential->N;
	vector<double> sum_aux;
	vector<double> exp_aux;

	#pragma omp parallel for private(sum_aux, exp_aux)
	for (int x1 = 0; x1 < N; x1++){
		for (int y3 = 0; y3 < N; y3++){
            sum_aux.reserve(N);
			for (int x2 = 0; x2 < N; x2++){
                exp_aux.reserve(N);
				// fill in exponents
                for (int y2 = 0; y2 < N; y2++){
                    exp_aux.push_back((*potential)(x1, y2, y3) + gaussian[abs(x2-y2)]);
				}
				sum_aux.push_back(get_summation(exp_aux));
                exp_aux.clear();
			}

			for (int x2 = 0; x2 < N; x2++){
				(*potential)(x1, x2, y3) = sum_aux[x2]; // potential is now f2
			}
			sum_aux.clear();
		}
	}
}

static void convolution_y3(Mesh* potential, const vector<double>& gaussian){
    // compute f3 (x1, x2, x3) == expotential(x, b)

    const int N = potential->N;
	vector<double> sum_aux;
    vector<double> exp_aux;

	#pragma omp parallel for private(sum_aux, exp_aux)
	for (int x1 = 0; x1 < N; x1++){
		for (int x2 = 0; x2 < N; x2++){
            sum_aux.reserve(N);
			for (int x3 = 0; x3 < N; x3++){
                exp_aux.reserve(N);
				// fill in exponents
                for (int y3 = 0; y3 < N; y3++){
                    exp_aux.push_back((*potential)(x1, x2, y3) + gaussian[abs(x3-y3)]);
				}
				sum_aux.push_back(get_summation(exp_aux));
                exp_aux.clear();
			}
			for (int x3 = 0; x3 < N; x3++){
				(*potential)(x1, x2, x3) = sum_aux[x3]; // potential is now f3
			}
			sum_aux.clear();
		}
	}
}

static void gen_expot(Mesh* potential,  const Mesh& expotential_0, double nu, double b)
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
	vector<double> gaussian(expotential_0.N);

	#pragma omp parallel for
	for (unsigned i = 0; i < expotential_0.N; i++){
		gaussian[i]=-i*i/(4.*b*nu);
	}

	convolution_y1(potential, gaussian, expotential_0);
	convolution_y2(potential, gaussian);
	convolution_y3(potential, gaussian);
}
