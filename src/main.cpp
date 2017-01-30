
#include "stdafx.h"
#include "core.h"
#include "core_power.h"
#include "zeldovich.h"
#include "frozen_flow.h"
#include "frozen_potential.h"

using namespace std;

int main(int argc, char* argv[]){
	struct timespec start, finish;
	double CPU_time, REAL_time;
	const clock_t START = clock();
	clock_gettime(CLOCK_MONOTONIC, &start);
	
	/* HANDLE COMMAND LINE OPTIONS */	
	
	Sim_Param sim; // simulation parameters	
	int err = sim.init(argc, argv);
	if (err) return err; // read command line options / config file
	norm_pwr(&sim.power); // compute power spectrum normalization
	sim.print_info(); // print simulation parameters
	
	// #include "examples.cpp"
	
	try{
		/* ZEL`DOVICH APPROXIMATION */
	//	err = zel_app(sim);
		
		/* FROZEN-FLOW APPROXIMATION */
	//	err = frozen_flow(sim);
	
		/* FROZEN-POTENTIL APPROXIMATION */
		err = frozen_potential(sim);
	}
	catch(int error){
		printf("ERROR %i!\n", error);
	}
	catch(...){
		printf("UNKNOWN ERROR!\n");
	}

	clock_gettime(CLOCK_MONOTONIC, &finish);
	CPU_time = (double)(clock() - START) / CLOCKS_PER_SEC;
	REAL_time = finish.tv_sec - start.tv_sec + (finish.tv_nsec - start.tv_nsec) / 1000000000.0;

	printf("\nProgram ran for %f s and used %f s of CPU time.\n\n", REAL_time, CPU_time);
	return err;
}
