#include <catch.hpp>
#include "test.hpp"
#include "core_power.cpp" ///< implementation testing

TEST_CASE( "UNIT TEST: growth functions {growth_factor, growth_rate, growth_change}", "[core_power]" )
{
    print_unit_msg("growth functions {growth_factor, growth_rate, growth_change}");

    int argc = 1;
    const char* const argv[1] = {"test"};
    try{
        Sim_Param sim(argc, argv);
        FTYPE_t D, D_to_a, f, Oma, OLa, dDda, factor;
        BOOST_LOG_TRIVIAL(debug) << "a\t\tD\t\tf\t\tdD/da\t\tD/a*f\t\tOm\t\tOm^0.6\t\tfactor";
        for (FTYPE_t a =0; a <= 1.0; a += 0.1)
        {
            D = growth_factor(a, sim.cosmo);
            f = growth_rate(a, sim.cosmo);
            dDda = growth_change(a, sim.cosmo);
            D_to_a = a ? D/a : dDda;
            OLa = 1/(1+sim.cosmo.Omega_m/(pow(a, 3)*sim.cosmo.Omega_L()));
            Oma = 1 - OLa;
            factor = sim.cosmo.Omega_m/(sim.cosmo.Omega_m+2*pow(a, 3)*sim.cosmo.Omega_L())*D_to_a;
            BOOST_LOG_TRIVIAL(debug) << a << "\t" << D << "\t" << f << "\t" << dDda << "\t" << D_to_a*f << "\t" << Oma << "\t" << pow(Oma, 0.6) << "\t" << factor;
        }
    }
    catch(const std::exception& e){
		BOOST_LOG_TRIVIAL(error) << "Error: " << e.what();
    }
}