#include "core.h"
#include "core_power.h"
#include "approximations.hpp"

using namespace std;

int main(int argc, char* argv[]){
	try{
    	struct timespec start, finish;
	    double CPU_time, REAL_time;
	    const clock_t START = clock();
	    clock_gettime(CLOCK_MONOTONIC, &start);

        /* SIMULATION PARAMETERS
            - read and handle command line options / config file
            - compute power spectrum normalization
        */
        Sim_Param sim(argc, argv);
        sim.print_info();
        
        do{
            /* ZEL`DOVICH APPROXIMATION */
            if(sim.comp_app.ZA)	zel_app(sim);
            
            /* FROZEN-FLOW APPROXIMATION */
            if(sim.comp_app.FF)	frozen_flow(sim);
        
            /* FROZEN-POTENTIAL APPROXIMATION */
            if(sim.comp_app.FP)	frozen_potential(sim);
            
            /* ADHESION APPROXIMATION */
            if(sim.comp_app.AA)	adhesion_approximation(sim);
            
            /* MODIFIED FROZEN-POTENTIAL APPROXIMATION */
            if(sim.comp_app.FP_pp)	mod_frozen_potential(sim);
        } while (sim.simulate());

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

        return 0;
	}
    catch(const string& e){
        if (e == "help") return 0;
        cerr << "Error: " << e << "\n";
        return 1;
    }
	catch(const exception& e){
		cerr << "Error: " << e.what() << "\n";
        return 1;
	}
	catch(...){
		cerr << "Exception of unknown type!\n";
        return 1;
	}
}
