#include <catch.hpp>

TEST_CASE( "UNIT TEST: periodicity functions {get_per}", "[core_mesh]" )
{
    CHECK( get_per(9.4, 10) == Approx(9.4) );
    CHECK( get_per(31.4, 10) == Approx(1.4) );
    CHECK( get_per(-7.4, 10) == Approx(2.6) );
    CHECK( get_per(10.0, 10) == Approx(0.0) );
    CHECK( get_per(0.0, 10) == Approx(0.0) );
    CHECK( get_per(-0.0, 10) == Approx(0.0) );
    CHECK( get_per(-10.0, 10) == Approx(0.0) );

    CHECK( get_per(9, 10) == 9 );
    CHECK( get_per(31, 10) == 1 );
    CHECK( get_per(-7, 10) == 3 );
    CHECK( get_per(10, 10) == 0 );
    CHECK( get_per(0, 10) == 0 );
    CHECK( get_per(-0, 10) == 0 );
    CHECK( get_per(-10, 10) == 0 );
}