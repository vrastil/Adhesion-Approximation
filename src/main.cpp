	
#include <ctime>
#include <fftw3.h>
#include "stdafx.h"
#include "cmd_line.h"
#include "grid_fce.h"
#include "output.h"
#include "approximations.h"

// #include "mesh.h"
// #include "mod_frozen_pot.h"

using namespace std;
const double PI = acos(-1.);

int par_num, mesh_num, Ng, box_size, nt;
double A, ns, k2_G, s8, z_in, z_out, nu;

string out_dir;

struct timespec start, finish;
double CPU_time, REAL_time;

void print_info(){
	printf("\n");	
	if (!isPowerOfTwo(par_num)) printf("WARNING: Given number of particles per dimension IS NOT a power of 2!\n");
	if (!isPowerOfTwo(mesh_num)) printf("WARNING: Given number of mesh cells per dimension IS NOT a power of 2!\n");
	printf("Num_par:\t%i^3\n", par_num);
	printf("Num_mesh:\t%i^3\n", mesh_num);
	printf("Box size:\t%i Mpc/h\n", box_size);
	printf("Starting redshift:\t%G\n", z_in);
	printf("The primordial power spectrum 'P(k)=A*k^ns' has amplitude A = %G and spectral index ns = %G.\n", A, ns);
//	printf("The transfer function 'T(k)=(1 + [ak + (bk)^1.5 + (ck)^2]^v)^(-1/v)'\nhas coefficients a = %G, b = %G, c = %G and v = %G.\n", a, b, c, v);
	if (k2_G == 0) printf("Smoothing length was not set.\n");
	else printf("Smoothing wavenumber is %G h/Mpc.\n", k2_G);
	k2_G *= k2_G;
	printf("'viscozity' for adhesion approximation is %G px^2.\n", nu);
	printf("The program will try to use %i threads.\n", nt);
	cout << "Output will be written to folder '"<< out_dir << "'\n";
	printf("\n");	
}

int main(int argc, char* argv[]){
	const clock_t START = clock();
	clock_gettime(CLOCK_MONOTONIC, &start);
	
	/* HANDLE COMMAND LINE OPTIONS */
	
	if (int err = handle_cmd_line(argc, argv)) return err;
	par_num = mesh_num / Ng;
	double parameters[]={A, ns, (double)box_size};
	power_spectrum_norm(parameters);
	A  = parameters[0];	
	print_info();
	
	/*
	c_Sim_Param sim; // simulation parameters
	int err = sim.init(argc, argv);
	if (err) return err; // read command line options / config file
	sim.power.norm_pwr(); // compute power spectrum normalization
	sim.print_info(); // print simulation parameters
	*/
	
	try{
		
		/* ZEL`DOVICH APPROXIMATION */
//		work_dir_over(out_dir + "ZA_run/");
//		zel_app(mesh_num, Ng, box_size, power_spectrum_T, parameters, z_in, z_out, nt, out_dir + "ZA_run/");

		/* FROZEN-FLOW APPROXIMATION */
		work_dir_over(out_dir + "FF_run/");
		frozen_flow(mesh_num, Ng, box_size, power_spectrum_T, parameters, z_in, z_out, nt, out_dir + "FF_run/");
		
		/* TRUNCATED ZEL`DOVICH APPROXIMATION */
//		work_dir_over(out_dir + "ZAt_run/");
//		truncated_zeldovich(mesh_num, Ng, box_size, power_spectrum_T, parameters, z_in, z_out, nt, out_dir + "ZAt_run/", k2_G);
		
		/* ADEHSION APPROXIMATION */
		work_dir_over(out_dir + "AA_run/");
		adhesiom_approximation(mesh_num, Ng, box_size, power_spectrum_T, parameters, nu, z_in, z_out, nt, out_dir + "AA_run/");
		
		/* MODIFIED FROZEN-FLOW APPROXIMATION */
//		work_dir_over(out_dir + "FFm_run/");
//		mod_frozen_flow(mesh_num, Ng, box_size, power_spectrum_T, parameters, z_in, z_out, nt, out_dir + "FFm_run/");	
	
		/* FROZEN-POTENTIAL APPROXIMATION */
		work_dir_over(out_dir + "FPA_run/");
		frozen_potential(mesh_num, Ng, box_size, power_spectrum_T, parameters, z_in, z_out, nt, out_dir + "FPA_run/");

		/* MODIFIED FROZEN-POTENTIAL APPROXIMATION */
//		sim.out_dir_app = sim.out_dir + "FPA_mod_run/";
//		err = mod_frozen_potential(sim);
//		printf("Modified frozen-potential approximation exited with status %i", err);
	}
	catch(...){
		printf("ERROR!\n");
	}

	clock_gettime(CLOCK_MONOTONIC, &finish);
	CPU_time = (double)(clock() - START) / CLOCKS_PER_SEC;
	REAL_time = finish.tv_sec - start.tv_sec + (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

	printf("\nProgram ran for %f s and used %f s of CPU time.\n\n", REAL_time, CPU_time);
	return 0;
}
