
#include "stdafx.h"
#include "core.h"
#include "core_out.h"
#include "core_app.h"
#include "core_mesh.h"

using namespace std;
const double PI = acos(-1.);

static void gen_init_expot(const Mesh& potential, Mesh* expotential, double nu);
static void pot2exp(const Mesh& potential,  Mesh* expotential, double nu);
static void exp2pot(Mesh* potential,  const Mesh& expotential, double nu);
static void gen_expot(Mesh* potential,  const Mesh& expotential, double nu, double b);
static void aa_convolution(App_Var_AA* APP, const Sim_Param &sim);


/**
 * @class:	Exp_sum
 * @brief:	class for handling summation of expoential with huge exponents
 */

class Exp_sum
{
public:
	// CONSTRUCTORS
	Exp_sum(): exponent(0), amplitude(1){};
    Exp_sum(double exponent): exponent(exponent), amplitude(1){};
	
	// VARIABLES
	double exponent, amplitude;
	
	// METHODS
    inline void reset() { exponent = 0; amplitude = 1; }
    inline double get_pure_exp(){ return exponent + log(amplitude); }

	// OPERATORS
	Exp_sum& operator+=(double exponent_rhs);
    inline Exp_sum& operator*=(double exponent_rhs){ exponent+= exponent_rhs}
};

int adhesion_approximation(const Sim_Param &sim)
{
	cout << "\n"
	"**********************\n"
	"ADHESION APPROXIMATION\n"
	"**********************\n";
	
	string out_dir_app = sim.out_dir + "AA_run/" + currentDateTime() + "/";
	work_dir_over(out_dir_app);
	sim.print_info(out_dir_app);
    
	/** ALLOCATION OF MEMORY + FFTW PREPARATION **/
	App_Var_AA APP(sim, "_AA_");
	printf("Initialization completed...\n");
	
	/** STANDARD PREPARATION FOR INTEGRATIOM **/
	
	/* Generating the right density distribution in k-space */	
	gen_rho_dist_k(sim, &APP.app_field[0], APP.p_F);
	
	/* Computing initial power spectrum in k-space */
	pwr_spec_k(sim, APP.app_field[0], &APP.power_aux);
	gen_pow_spec_binned(sim, APP.power_aux, &APP.pwr_spec_binned_0);
	print_pow_spec(APP.pwr_spec_binned_0, out_dir_app,  APP.z_suffix_const + "init");
	
	/* Computing initial potential in k-space */
	gen_pot_k(&APP.app_field[0]);

	/* Computing initial expotential */
	fftw_execute_dft_c2r(APP.p_B, APP.app_field[0]);
	gen_init_expot(APP.app_field[0], &APP.expotential, sim.nu, APP.p_F);

	/* Setting initial positions of particles */
    printf("Setting initial positions of particles...\n");
	set_unpert_pos(sim, APP.particles);
	aa_convolution(&APP, sim);
    upd_pos_first_order(sim, sim.b_in, APP.particles, APP.app_field);
    APP.print(sim, out_dir_app);
	APP.upd_time();

	/** INTEGRATION **/
	
	while(APP.integrate())
	{
		printf("\nStarting computing step with z = %.2f (b = %.3f)\n", APP.z(), APP.b);
		/* Computing convolution */
        aa_convolution(&APP, sim);
		
		/* Updating positions of particles... */
		printf("Updating positions of particles...\n");
		upd_pos_first_order(sim, APP.db, APP.particles, APP.app_field);
		
		if (APP.printing()) APP.print(sim, out_dir_app);
		APP.upd_time();
	}
	print_suppression(APP.supp, sim, out_dir_app);
		
	printf("Adhesion approximation ended successfully.\n");
	return APP.err;
}

static void aa_convolution(App_Var_AA* APP, const Sim_Param &sim)
{
    //	gen_expot(&APP->app_field[0], APP->expotential, sim.nu, APP->b_half());
	gen_expot(&APP->app_field[0], APP->expotential, sim.nu, APP->b;
				
	printf("Computing potential...\n");	
	exp2pot(&APP->app_field[0], APP->app_field[0], sim.nu);
				
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

static void convolution_y1(Mesh* potential, const vector<double>& gaussian, const Mesh& expotential_0){
	// multi-thread index is y3
	// compute f1 (x1, y2, y3)
	double exponent;
	double sum;
	
	for (int x1 = 0; x1 < potential->N; x1++){
		for (int y2 = 0; y2 < potential->N; y2++){
			#pragma omp parallel for private(exponent, sum)
			for (int y3 = 0; y3 < potential->N; y3++){
				// summation over y1
				sum = 0;
				for (int y1 = 0; y1 < potential->N; y1++){
					exponent=expotential_0(y1, y2, y3)+pow(x1-y1, 2.)/(2.*b);
					exponent = -exponent/(2.*nu);
					if (exponent < 500) sum +=exp(exponent);
					else printf("WARNING! Huge exponent = %f\n", exponent);
					if (!isfinite(sum)) {printf("Error while chcecking for NAN or INF. r = %i, potential = %f, b = %f, exponent = %f, +sum = %f\n",
						abs(x1-y1), expotential_0(y1, y2, y3), b, exponent, exp(exponent)); }
				}
				(*potential)(x1, y2, y3) = sum; // potential is now f1
			}
		}
	}
}

static void convolution_y2(Mesh* potential, const vector<double>& gaussian){
	// multi-thread index is x1
	// compute f2 (x1, x2, y3)
	double sum;
	vector<double> sum_aux(potential->N);
	
	#pragma omp parallel for private(sum_aux, sum)
	for (int x1 = 0; x1 < potential->N; x1++){
		for (int y3 = 0; y3 < potential->N; y3++){
			for (int x2 = 0; x2 < potential->N; x2++){
				// summation over y2
				sum = 0;
				for (int y2 = 0; y2 < potential->N; y2++){
					sum += (*potential)(x1, y2, y3)*gaussian[abs(x2-y2)];
					if (!isfinite(sum)) {printf("Error while chcecking for NAN or INF. r = %i, potential = %f, gaussian = %f, +sum = %f\n",
						abs(x2-y2), (*potential)(x1, y2, y3), gaussian[abs(x2-y2)], (*potential)(x1, y2, y3)*gaussian[abs(x2-y2)]); }
				}
				sum_aux.push_back(sum);
			}

			for (int x2 = 0; x2 < potential->N; x2++){
				(*potential)(x1, x2, y3) = sum_aux[x2]; // potential is now f2
			}
			
			sum_aux.clear();
		}
	}
}

static void convolution_y3(Mesh* potential, const vector<double>& gaussian){
	// multi-thread index is x1
	// compute f3 (x1, x2, x3) == expotential(x, b)
	double sum;
	vector<double> sum_aux(potential->N);

	#pragma omp parallel for private(sum_aux, sum)
	for (int x1 = 0; x1 < potential->N; x1++){
		for (int x2 = 0; x2 < potential->N; x2++){
			for (int x3 = 0; x3 < potential->N; x3++){
				// summation over y3
				sum = 0;
				for (int y3 = 0; y3 < potential->N; y3++){
					sum += (*potential)(x1, x2, y3)*gaussian[abs(x3-y3)];
					if (!isfinite(sum)) {printf("Error while chcecking for NAN or INF. r = %i, potential = %f, gaussian = %f, +sum = %f\n",
						abs(x3-y3), (*potential)(x1, x2, y3), gaussian[abs(x3-y3)], (*potential)(x1, x2, y3)*gaussian[abs(x3-y3)]); }
				}
				sum_aux.push_back(sum);
			}
			for (int x3 = 0; x3 < potential->N; x3++){
				(*potential)(x1, x2, x3) = sum_aux[x3]; // potential is now f3
				if (sum_aux[x3] < 0) printf("WARNING! Negative value of convolution (%f)", sum_aux[x3]);
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
	vector<double> gaussian(expotential.N);

	#pragma omp parallel for
	for (int i = 0; i < expotential_0.N; i++){
		gaussian[i]=-i*i/(4.*b*nu);
	}

	convolution_y1(potential, gaussian, expotential_0);
	convolution_y2(potential, gaussian);
	convolution_y3(potential, gaussian);
}

Exp_sum& Exp_sum::operator+=(double exponent_rhs)
{
    if (exponent >= exponent_rhs)
    {
        amplitude += exp(exponent_rhs - exponent);
    } else {
        amplitude = 1 + exp(exponent - exponent_rhs);
        exponent = exponent_rhs
    }
	return *this;
}