/**
 * @brief main file
 * 
 * @file main.cpp
 * @author Michal Vrastil
 * @date 2018-07-11
 */

#include "stdafx.h"

#include "params.hpp"
#include "adhesion.hpp"
#include "chameleon.hpp"
#include "frozen_flow.hpp"
#include "frozen_potential.hpp"
#include "mod_frozen_potential.hpp"
#include "zeldovich.hpp"

#include <boost/log/core.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/timer/timer.hpp>

namespace logging = boost::log;
namespace sinks = boost::log::sinks;

namespace
{
void init_logging()
{
    logging::core::get()->set_filter
    (
        logging::trivial::severity >= logging::trivial::info
    );
}

class Timer{
public:
   Timer() = default;
   ~Timer(){
       t.stop();
       BOOST_LOG_TRIVIAL(info) << "Time elapsed:" << t.format();
   }
private:
    boost::timer::cpu_timer t;
};

template<class T>
void init_and_run_app(Sim_Param& sim)
{
    Timer t;
    T APP(sim);
    APP.update_cosmo(sim.cosmo);
    APP.run_simulation();   
}

}
/**
 * @brief initialize program and run all simulations
 * 
 * @param argc number of command line arguments
 * @param argv command line arguments
 * @return int exit status
 */
int main(int argc, char* argv[]){
    // initialize logging
    init_logging();

    // timer --automatically displays timing information
    Timer t;

	try{
        /* SIMULATION PARAMETERS
            - read and handle command line options / config file
            - compute power spectrum normalization
        */
        Sim_Param sim(argc, argv);
        sim.print_info();
        
        do{
            /* ZEL`DOVICH APPROXIMATION */
            if(sim.comp_app.ZA)	init_and_run_app<App_Var_ZA>(sim);

            /* TRUNCATED ZEL`DOVICH APPROXIMATION */
            if(sim.comp_app.TZA) init_and_run_app<App_Var_TZA>(sim);
            
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
        BOOST_LOG_TRIVIAL(info) << "All simulations completed.";
        return 0;
	}
    catch(const std::string& e){
        if (e == "help") return 0;
        BOOST_LOG_TRIVIAL(error) << "Error: " << e;
        return 1;
    }
	catch(const std::exception& e){
		BOOST_LOG_TRIVIAL(error) << "Error: " << e.what();
        return 1;
	}
	catch(...){
		BOOST_LOG_TRIVIAL(error) << "Exception of unknown type!";
        return 2;
	}
}
