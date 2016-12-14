
#include "stdafx.h"
#include "core.h"

using namespace std;

struct timespec start, finish;
double CPU_time, REAL_time;

int main(int argc, char* argv[]){
	const clock_t START = clock();
	clock_gettime(CLOCK_MONOTONIC, &start);
	
	/* HANDLE COMMAND LINE OPTIONS */
	
	/*
	c_Sim_Param sim; // simulation parameters
	int err = sim.init(argc, argv);
	if (err) return err; // read command line options / config file
	sim.power.norm_pwr(); // compute power spectrum normalization
	sim.print_info(); // print simulation parameters
	*/
	
	#include "examples.cpp"
	
	try{

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
