
#include "stdafx.h"
#include "core.h"
#include "core_out.h"
#include "core_app.h"
#include "core_mesh.h"

// #define CONV_MOD 0 // convolution using FFT
#define CONV_MOD 1 // convolution using direct sum
// #define CONV_MOD 2 // TODO: convolution using better FFT

using namespace std;
const double PI = acos(-1.);

static void gen_init_expot(const Mesh& potential, Mesh* expotential, double nu, const fftw_plan &p_F);
static void pot2exp(const Mesh& potential,  Mesh* expotential, double nu);
static void exp2pot(Mesh* potential,  const Mesh& expotential, double nu);
static void gen_expot(Mesh* potential,  const Mesh& expotential, double nu, double b, const fftw_plan &p_B);

static int check_field(const Mesh& field){
	// check field of length max_i*max_i*(max_i+2); ignore last 2 element in last dim
	double check;
	for (int i = 0; i < field.N; i++){
		for (int j = 0; j < field.N; j++){
			for (int k = 0; k < field.N; k++){
				check = field(i, j, k);
				if (isfinite(check) == 0){
					printf("Error while performing check for NAN or INF! field posititon = (%i, %i, %i), value = %f\n", 
					i, j, k, check);
					return 1;
				}
			}
		}
	}
	return 0;
}

static int check_field_pos(const Mesh& field){
	// check field of length max_i*max_i*(max_i+2); ignore last 2 element in last dim
	double check;
	for (int i = 0; i < field.N; i++){
		for (int j = 0; j < field.N; j++){
			for (int k = 0; k < field.N; k++){
				check = field(i, j, k);
				if (check < 0){
					printf("Error while performing check for negative values! field posititon = (%i, %i, %i), value = %f\n", 
					i, j, k, check);
					return 1;
				}
			}
		}
	}
	return 0;
}

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
	if (check_field(APP.app_field[0])) return 1;

	/* Setting initial positions of particles */
    printf("Setting initial positions of particles...\n");
	set_unpert_pos(sim, APP.particles);
	
	/** INTEGRATION **/
	
	while(APP.integrate())
	{
		printf("\nStarting computing step with z = %.2f (b = %.3f)\n", APP.z(), APP.b);
		
		/* Computing convolution */
	//	gen_expot(&APP.app_field[0], APP.expotential, sim.nu, APP.b_half(), APP.p_B);
		gen_expot(&APP.app_field[0], APP.expotential, sim.nu, APP.b, APP.p_B);
		if (check_field(APP.app_field[0])) return 1;
		if (check_field_pos(APP.app_field[0])) return 1;
				
		printf("Computing potential...\n");	
		exp2pot(&APP.app_field[0], APP.app_field[0], sim.nu);
		if (check_field(APP.app_field[0])) return 1;
				
		printf("Computing velocity field via FFT...\n");
		fftw_execute_dft_r2c(APP.p_F, APP.app_field[0]);
		gen_displ_k(&APP.app_field, APP.app_field[0]);
		fftw_execute_dft_c2r_triple(APP.p_B, APP.app_field);
		
		/* Updating positions of particles... */
		printf("Updating positions of particles...\n");
		upd_pos_first_order(sim, APP.db, APP.particles, APP.app_field);
		
		if (APP.printing())
	//	if (false)
		{
			/* Printing positions */
			print_par_pos_cut_small(APP.particles, sim, out_dir_app, APP.z_suffix());
			APP.track.update_track_par(APP.particles);
			print_track_par(APP.track, sim, out_dir_app, APP.z_suffix());

			/* Printing density */
			get_rho_from_par(APP.particles, &APP.power_aux, sim);
			gen_dens_binned(APP.power_aux, APP.dens_binned, sim);
			print_rho_map(APP.power_aux, sim, out_dir_app, APP.z_suffix());
			print_dens_bin(APP.dens_binned, sim.mesh_num, out_dir_app, APP.z_suffix());
			
			/* Printing power spectrum */
			fftw_execute_dft_r2c(APP.p_F, APP.power_aux);
			pwr_spec_k(sim, APP.power_aux, &APP.power_aux);
			gen_pow_spec_binned(sim, APP.power_aux, &APP.pwr_spec_binned);
			print_pow_spec(APP.pwr_spec_binned, out_dir_app, APP.z_suffix());
			print_pow_spec_diff(APP.pwr_spec_binned, APP.pwr_spec_binned_0, APP.b, out_dir_app, APP.z_suffix());
			
			APP.upd_supp();
		}
		APP.upd_time();
	}
	print_suppression(APP.supp, sim, out_dir_app);
		
	printf("Adhesion approximation ended successfully.\n");
	return APP.err;
}

static void gen_init_expot(const Mesh& potential, Mesh* expotential, double nu, const fftw_plan &p_F)
{
	switch(CONV_MOD){
	case 0:
		printf("Computing initial expotential in k-space...\n");
		pot2exp(potential, expotential, nu);	
		fftw_execute_dft_r2c(p_F, *expotential);
		break;
	case 1:
		printf("Storing initial potenital in q-space...\n");
		*expotential = potential;
		break;
	}
}

static void pot2exp(const Mesh& potential,  Mesh* expotential, double nu)
{
	#pragma omp parallel for
	for (int i = 0; i < potential.N*potential.N; i++)
	{
		for (int j = 0; j < potential.N; j++)
		{
			(*expotential)(i, j) = exp(-potential(i, j)/(2.*nu));
		}
	}
}

static void exp2pot(Mesh* potential,  const Mesh& expotential, double nu)
{
	#pragma omp parallel for
	for (int i = 0; i < expotential.N*expotential.N; i++)
	{
		for (int j = 0; j < expotential.N; j++)
		{
			(*potential)(i, j) = -2.*nu*log(expotential(i, j));
			if (!isfinite((*potential)(i, j))) {printf(
			"Error while chcecking for NAN or INF. expotential = %f, -2.*nu*log(U) = %f, potential = %f\n",
						expotential(i, j),-2.*nu*log(expotential(i, j)), (*potential)(i, j)); }
		}
	}
}

static void convolution_y1(Mesh* potential,  const Mesh& expotential_0, double nu, double b){
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
					else printf("WARNING! Huge exponent = %f", exponent);
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
	// compute f1 (x1, y2, y3)
	double sum;
	vector<double> sum_aux(potential->N);
	
	// compute f2 (x1, x2, y3)
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
	// compute f1 (x1, y2, y3)
	double sum;
	vector<double> sum_aux(potential->N);

	// compute f3 (x1, x2, x3) == expotential(x, b)
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

static void gen_expot(Mesh* potential,  const Mesh& expotential, double nu, double b, const fftw_plan &p_B)
{
	switch(CONV_MOD){
	case 0:
		/* Computing convolution using FFT */
		printf("Computing expotential in k-space...\n");
		double k2;
		#pragma omp parallel for private(k2)
		for(int i=0; i < expotential.length / 2;i++){
			k2 = get_k_sq(expotential.N, i);
			(*potential)[2*i] = expotential[2*i]*exp(-4.*PI*PI*b*nu*k2/(expotential.N*expotential.N)); // DEFINITION of nu_phys = (L/N)^2 * v_comp
			(*potential)[2*i+1] = expotential[2*i+1]*exp(-4.*PI*PI*b*nu*k2/(expotential.N*expotential.N));
		}
		fftw_execute_dft_c2r(p_B, *potential);
		break;
	case 1:
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
		for (int i = 0; i < expotential.N; i++){
			gaussian[i]=exp(-i*i/(4.*b*nu));
		}

		convolution_y1(potential, expotential, nu, b);
		convolution_y2(potential, gaussian);
		convolution_y3(potential, gaussian);
		break;
	}
}