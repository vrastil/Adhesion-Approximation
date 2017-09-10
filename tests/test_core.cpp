#include <catch.hpp>

TEST_CASE( "sign function is working properly", "[sign]" ){

    REQUIRE( sgn<int>(0) == 0);
    REQUIRE( sgn<int>(10) == 1);
    REQUIRE( sgn<int>(-5) == -1);
    REQUIRE( sgn<double>(0.) == 0);
    REQUIRE( sgn<double>(0.34534535) == 1);
    REQUIRE( sgn<double>(-1.34534E16) == -1);

}
