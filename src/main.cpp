/*
	Main part of the program. Functionalities:
		- handles command line options
		- Zel`dovich approximation
		- computes power spectrum
*/
#include <ctime>
#include <fftw3.h>
#include "stdafx.h"
#include "cmd_line.h"
#include "grid_fce.h"
#include "output.h"
#include "approximations.h"

using namespace std;
const double PI = acos(-1.);

int par_num, mesh_num, box_size, nt;
double A, ns, k2_G, s8, z;

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
	printf("The simulation is starting at redshift:\t%G\n", z);
	printf("The primordial power spectrum 'P(k)=A*k^ns' has amplitude A = %G and spectral index ns = %G.\n", A, ns);
//	printf("The transfer function 'T(k)=(1 + [ak + (bk)^1.5 + (ck)^2]^v)^(-1/v)'\nhas coefficients a = %G, b = %G, c = %G and v = %G.\n", a, b, c, v);
	if (k2_G == 0) printf("Smoothing length was not set.\n");
	else printf("Smoothing wavenumber is %G h/Mpc.\n", k2_G);
	k2_G *= k2_G;
	printf("The program will try to use %i threads.\n", nt);
//	out_dir = work_dir(out_dir);
	out_dir = "output/run/";
	work_dir_over(out_dir);
	cout << "Output will be written to folder '"<< out_dir << "'\n";
	printf("\n");	
}

int main(int argc, char* argv[]){
	const clock_t START = clock();
	clock_gettime(CLOCK_MONOTONIC, &start);
	
	/* HANDLE COMMAND LINE OPTIONS */
	if (int err = handle_cmd_line(argc, argv)) return err;	
	double parameters[]={A, ns, (double)box_size};
	power_spectrum_norm(parameters);
	A  = parameters[0];	
	print_info();
	
	try{		
		/* ZEL`DOVICH APPROXIMATION */
//		zel_app(mesh_num, 1, box_size, power_spectrum_T, parameters, z, 10., nt, out_dir);

		/* FROZEN-FLOW APPROXIMATION */
	//	parameters[2]=1e3;
		frozen_flow(mesh_num, 1, box_size, single_power_spectrum, parameters, z, 10., nt, out_dir);
		
		/* TRUNCATED ZEL`DOVICH APPROXIMATION */
//		truncated_zeldovich(mesh_num, 1, box_size, power_spectrum_T, parameters, z, 0., nt, out_dir, k2_G);
		
		/* ADEHSION APPROXIMATION */
		
//		const double nu = 1.0 * pow((double)mesh_num / box_size, 2.);
//		adhesiom_approximation(mesh_num, 1, box_size, power_spectrum_T, parameters, nu, z, 0., nt, out_dir);
		
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
