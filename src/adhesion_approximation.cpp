
#include "stdafx.h"
#include "core.h"
#include "core_out.h"
#include "core_app.h"
#include "core_mesh.h"

using namespace std;

const double ACC = 1e-10;
const double log_acc = log(ACC);

static void gen_init_expot(const Mesh& potential, Mesh* expotential, double nu);
static void gen_expot(Mesh* potential,  const Mesh& expotential, double nu, double b);
static void aa_convolution(App_Var_AA* APP, const Sim_Param &sim);

int adhesion_approximation(const Sim_Param &sim)
{
	cout << "\n"
	"**********************\n"
	"ADHESION APPROXIMATION\n"
	"**********************\n";
	string out_dir_app = std_out_dir("AA_run/", sim);
	work_dir_over(out_dir_app);
    
	/******************************************
    * ALLOCATION OF MEMORY + FFTW PREPARATION *
    ******************************************/

	App_Var_AA APP(sim, "_AA_");
	printf("Initialization completed...\n");
	
	/***************************************
    * STANDARD PREPARATION FOR INTEGRATIOM *
    ***************************************/
	
	/* Generating the right density distribution in k-space */	
	gen_rho_dist_k(sim, &APP.app_field[0], APP.p_F);
	
	/* Computing initial potential in k-space */
	gen_pot_k(&APP.app_field[0]);

	/* Computing initial expotential */
	fftw_execute_dft_c2r(APP.p_B, APP.app_field[0]);
	gen_init_expot(APP.app_field[0], &APP.expotential, sim.nu);

	/*********************
    * INITIAL CONDITIONS *
    *********************/

    printf("\nSetting initial positions of particles...\n");
	set_unpert_pos(sim, APP.particles);
	aa_convolution(&APP, sim);
    upd_pos_first_order(sim, sim.b_in, APP.particles, APP.app_field);

    /* Setting initial (binned) power spectrum, WARNING: power_aux is modified */
    APP.track.update_track_par(APP.particles);
    APP.print(sim, out_dir_app);
	APP.upd_time();

	/**************
    * INTEGRATION *
    **************/
	
	while(APP.integrate())
	{
		printf("\nStarting computing step with z = %.2f (b = %.3f)\n", APP.z(), APP.b);
		/* Computing convolution */
        aa_convolution(&APP, sim);
		
		/* Updating positions of particles... */
		printf("Updating positions of particles...\n");
		upd_pos_first_order(sim, APP.db, APP.particles, APP.app_field);
        
        APP.track.update_track_par(APP.particles);
		if (APP.printing()) APP.print(sim, out_dir_app);
		APP.upd_time();
	}
    
    sim.print_info(out_dir_app, "AA");
	printf("Adhesion approximation ended successfully.\n");
	return APP.err;
}

static void aa_convolution(App_Var_AA* APP, const Sim_Param &sim)
{
    printf("Computing potential...\n");	
    //	gen_expot(&APP->app_field[0], APP->expotential, sim.nu, APP->b_half());
	gen_expot(&APP->app_field[0], APP->expotential, sim.nu, APP->b);
    APP->app_field[0] *= -2*sim.nu;
				
	printf("Computing velocity field via FFT...\n");
	fftw_execute_dft_r2c(APP->p_F, APP->app_field[0]);
	gen_displ_k(&APP->app_field, APP->app_field[0]);
	fftw_execute_dft_c2r_triple(APP->p_B, APP->app_field);
}

static void gen_init_expot(const Mesh& potential, Mesh* expotential, double nu)
{
	printf("Storing initial expotenital in q-space...\n");
    // store exponent only
	*expotential = potential;
    *expotential /= -2*nu;
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

    const unsigned N = potential->N;
    vector<double> exp_aux(N);
    
	#pragma omp parallel for private(exp_aux)
	for (unsigned x1 = 0; x1 < N; x1++){
		for (unsigned y2 = 0; y2 < N; y2++){
			for (unsigned y3 = 0; y3 < N; y3++){
                // fill in exponents
                for (unsigned y1 = 0; y1 < N; y1++){
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
    const unsigned N = potential->N;
	vector<double> sum_aux(N);
	vector<double> exp_aux(N);

	#pragma omp parallel for private(sum_aux, exp_aux)
	for (unsigned x1 = 0; x1 < N; x1++){
		for (unsigned y3 = 0; y3 < N; y3++){
			for (unsigned x2 = 0; x2 < N; x2++){
				// fill in exponents
                for (unsigned y2 = 0; y2 < N; y2++){
                    exp_aux.push_back((*potential)(x1, y2, y3) + gaussian[abs(x2-y2)]);
				}
				sum_aux.push_back(get_summation(exp_aux));
                exp_aux.clear();
			}

			for (unsigned x2 = 0; x2 < N; x2++){
				(*potential)(x1, x2, y3) = sum_aux[x2]; // potential is now f2
			}
			sum_aux.clear();
		}
	}
}

static void convolution_y3(Mesh* potential, const vector<double>& gaussian){
    // compute f3 (x1, x2, x3) == expotential(x, b)
    const unsigned N = potential->N;
	vector<double> sum_aux(N);
    vector<double> exp_aux(N);

	#pragma omp parallel for private(sum_aux, exp_aux)
	for (unsigned x1 = 0; x1 < N; x1++){
		for (unsigned x2 = 0; x2 < N; x2++){
			for (unsigned x3 = 0; x3 < N; x3++){
				// fill in exponents
                for (unsigned y3 = 0; y3 < N; y3++){
                    exp_aux.push_back((*potential)(x1, x2, y3) + gaussian[abs(x3-y3)]);
				}
				sum_aux.push_back(get_summation(exp_aux));
                exp_aux.clear();
			}
			for (unsigned x3 = 0; x3 < N; x3++){
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
