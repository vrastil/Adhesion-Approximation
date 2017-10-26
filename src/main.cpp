
#include "stdafx.h"
#include "core.h"
#include "core_power.h"
#include "approximations.hpp"

using namespace std;

int main(int argc, char* argv[]){
	struct timespec start, finish;
	double CPU_time, REAL_time;
	const clock_t START = clock();
	clock_gettime(CLOCK_MONOTONIC, &start);
	
	int err = 0;
	try{
        /* SIMULATION PARAMETERS
            - read and handle command line options / config file
            - compute power spectrum normalization
        */
        Sim_Param sim(argc, argv);
        sim.print_info();

        // double D, f, Oma, dDda;
        // printf("a\t\tD\t\tf\t\tdD/da\t\tD/a*f\t\tOm\t\tOm^0.6\n");
        // for (double a =0; a <= 0.1; a += 0.01)
        // {
        //     D = growth_factor(a, sim.power);
        //     f = growth_rate(a, sim.power);
        //     dDda = growth_change(a, sim.power);
        //     Oma = 1 - 1/(1+sim.power.Omega_m()/(pow(a, 3)*sim.power.Omega_L()));
        //     printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\n", a, D, f, dDda, D/a*f, Oma, pow(Oma, 0.6));
        // }
		
		/* ZEL`DOVICH APPROXIMATION */
		if(sim.comp_ZA)	err = zel_app(sim);
		
		/* FROZEN-FLOW APPROXIMATION */
		if(sim.comp_FF)	err = frozen_flow(sim);
	
		/* FROZEN-POTENTIAL APPROXIMATION */
		if(sim.comp_FP)	err = frozen_potential(sim);
		
		/* ADHESION APPROXIMATION */
		if(sim.comp_AA)	err = adhesion_approximation(sim);
		
		/* MODIFIED FROZEN-POTENTIAL APPROXIMATION */
		if(sim.comp_FP_pp)	err = mod_frozen_potential(sim);
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

    int CPU_int, CPU_dec, REAL_int, REAL_dec;

    CPU_int = floor(log10(CPU_time)) + 1;
    CPU_dec = 3 - CPU_int;
    if (CPU_int < 0) CPU_int = 0;
    if (CPU_dec < 0) CPU_dec = 0;

    REAL_int = floor(log10(REAL_time)) + 1;
    REAL_dec = 3 - REAL_int;
    if (REAL_int < 0) REAL_int = 0;
    if (REAL_dec < 0) REAL_dec = 0;

	printf("\nProgram ran for %.*fs and used %.*fs of CPU time.\n\n", \
           REAL_dec, REAL_time, CPU_dec, CPU_time);
	return err;
}
