
#include <ctime>
#include "stdafx.h"

#include "params.hpp"
#include "ApproximationsSchemes/adhesion.hpp"
#include "ApproximationsSchemes/chameleon.hpp"
#include "ApproximationsSchemes/frozen_flow.hpp"
#include "ApproximationsSchemes/frozen_potential.hpp"
#include "ApproximationsSchemes/mod_frozen_potential.hpp"
#include "ApproximationsSchemes/zeldovich.hpp"

#ifdef MPI_ENABLED
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
namespace mpi = boost::mpi;
#endif

template<class T>
static void init_and_run_app(const Sim_Param& sim)
{
    T APP(sim);
    APP.run_simulation();  
}

int main(int argc, char* argv[]){
    #ifdef MPI_ENABLED
    mpi::environment env;
    mpi::communicator world;
    std::cout << "I am process " << world.rank() << " of " << world.size()
            << "." << std::endl;

    if (world.rank() == 0)
    {
    #endif
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
            if(sim.comp_app.ZA)	init_and_run_app<App_Var_ZA>(sim);
            
            /* FROZEN-FLOW APPROXIMATION */
            if(sim.comp_app.FF)	init_and_run_app<App_Var_FF>(sim);
        
            /* FROZEN-POTENTIAL APPROXIMATION */
            if(sim.comp_app.FP)	init_and_run_app<App_Var_FP>(sim);
            
            /* ADHESION APPROXIMATION */
            if(sim.comp_app.AA)	init_and_run_app<App_Var_AA>(sim);
            
            /* MODIFIED FROZEN-POTENTIAL APPROXIMATION */
            if(sim.comp_app.FP_pp)	init_and_run_app<App_Var_FP_mod>(sim);

            /* CHAMELEON GRAVITY (FROZEN-POTENTIAL APPROXIMATION) */
            if(sim.comp_app.chi) init_and_run_app<App_Var_Chi>(sim);

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
    catch(const std::string& e){
        if (e == "help") return 0;
        std::cerr << "Error: " << e << "\n";
        return 1;
    }
	catch(const std::exception& e){
		std::cerr << "Error: " << e.what() << "\n";
        return 1;
	}
	catch(...){
		std::cerr << "Exception of unknown type!\n";
        return 1;
	}
     #ifdef MPI_ENABLED
    }
    #endif
}
