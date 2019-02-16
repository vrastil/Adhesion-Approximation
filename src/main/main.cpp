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
#include <boost/log/utility/setup/common_attributes.hpp>
#include <boost/log/utility/setup/console.hpp>
#include <boost/log/utility/setup/file.hpp>
#include <boost/timer/timer.hpp>
#include <boost/date_time/posix_time/posix_time_types.hpp>

namespace logging = boost::log;
namespace sinks = boost::log::sinks;
namespace keywords = boost::log::keywords;

typedef sinks::synchronous_sink< sinks::text_file_backend > sink_t;
typedef sinks::synchronous_sink< sinks::text_ostream_backend > sink_os_t;

namespace
{
void init_logging()
{
    // set logging core to log everything
    logging::core::get()->set_filter
    (
        logging::trivial::severity >= logging::trivial::trace
    );

    // add atributes
    boost::log::add_common_attributes();
    boost::log::register_simple_formatter_factory< boost::log::trivial::severity_level, char >("Severity");

    // register 'std:cout' as sink
    boost::shared_ptr< sink_os_t > g_file_sink = logging::add_console_log(std::cout);
    g_file_sink->set_filter(logging::trivial::severity >= logging::trivial::info);
}

class Logger
{
public:
    Logger(const std::string& filename):
        g_file_sink(logging::add_file_log(
            keywords::file_name = filename,
            keywords::time_based_rotation = sinks::file::rotation_at_time_interval(boost::posix_time::hours(1)),
            keywords::format = "[%TimeStamp%] <%Severity%>: %Message%"
        ))
    {
        g_file_sink->set_filter(logging::trivial::severity >= logging::trivial::debug);
    }
    ~Logger()
    {
        // close the log file
        logging::core::get()->remove_sink(g_file_sink);
        g_file_sink.reset();
    }

private:
    boost::shared_ptr< sink_t > g_file_sink;
};

class Timer{
public:
   Timer() = default;
   ~Timer(){
       stop();
   }
   void stop(bool print = true){
       if (!_t.is_stopped()){
            _t.stop();
            if (print) BOOST_LOG_TRIVIAL(info) << "Time elapsed:" << _t.format();
       }
   }
private:
    boost::timer::cpu_timer _t;
};

template<class T>
void init_and_run_app(Sim_Param& sim)
{
    // timer for wall time, CPU time
    Timer t;

    // create approximation class and whether we have truncation in initial power spectrum
    T APP(sim);
    APP.update_cosmo(sim.cosmo);

    // register log file
    Logger l(APP.get_out_dir() + "log/log" + "_%N.log");
    
    // run the simulation
    APP.run_simulation();

    // stop timer manually so it is in the log file
    t.stop();
}

} ///< end of anonymous namespace

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

    /* SIMULATION PARAMETERS
        - read and handle command line options / config file
        - compute power spectrum normalization
    */
    Sim_Param sim(argc, argv);
    if (sim.is_ready())
    {
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
    }
    else
    {
        t.stop(false);
    }
    return 0;
}
