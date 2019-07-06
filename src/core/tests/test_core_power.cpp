#include <catch.hpp>

#define BOOST_LOG_DYN_LINK 1
#include <boost/log/trivial.hpp>

#include "../core_power.cpp"

TEST_CASE( "UNIT TEST: growth functions {growth_factor, growth_rate, growth_change}", "[core_power]" )
{
    BOOST_LOG_TRIVIAL(info) << "growth functions {growth_factor, growth_rate, growth_change}";

    int argc = 1;
    const char* const argv[1] = {"test"};
    Sim_Param sim(argc, argv);
    FTYPE_t D, D_to_a, f, Oma, OLa, dDda, factor;
    for (FTYPE_t a =0; a <= 1.0; a += 0.1)
    {
        D = growth_factor(a, sim.cosmo);
        f = growth_rate(a, sim.cosmo);
        dDda = growth_change(a, sim.cosmo);
        D_to_a = a ? D/a : dDda;
        CHECK( dDda == Approx(D_to_a*f) );
    }

    CHECK(0.0 == growth_factor(0, sim.cosmo));
    CHECK(1.0 == growth_factor(1, sim.cosmo));
}