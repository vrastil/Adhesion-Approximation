#include <catch.hpp>
#include "test.hpp"

TEST_CASE( "UNIT TEST: sign function {sgn<T>}", "[core]" )
{
    print_unit_msg("sign function {sgn<T>}");

    CHECK( sgn(0) == 0 );
    CHECK( sgn(10) == 1 );
    CHECK( sgn(-5) == -1 );
    CHECK( sgn(0.) == 0 );
    CHECK( sgn(0.34534f) == 1 );
    CHECK( sgn(-1.34534E16) == -1 );
}

TEST_CASE( "UNIT TEST: periodicity functions {get_per}", "[core_mesh]" )
{
    print_unit_msg("periodicity functions {get_per}");

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

    Vec_3D<double> pos(0., -10., 10.);
    get_per(pos, 10);
    Vec_3D<double> pos2(4.3, -7.8, 18.4);
    get_per(pos2, 10);

    CHECK( pos[0] == Approx(0.) );
    CHECK( pos[1] == Approx(0.) );
    CHECK( pos[2] == Approx(0.) );
    CHECK( pos2[0] == Approx(4.3) );
    CHECK( pos2[1] == Approx(2.2) );
    CHECK( pos2[2] == Approx(8.4) );

    Vec_3D<int> posi(0, -10, 10);
    get_per(posi, 10);
    Vec_3D<int> pos2i(4, -7, 18);
    get_per(pos2i, 10);

    CHECK( posi == Vec_3D<int>(0, 0, 0) );
    CHECK( pos2i == Vec_3D<int>(4, 3, 8) );    

    Vec_3D<int> pos4i(0, -10, 10);
    get_per(pos4i, 5, 10, 4);
    CHECK( pos4i == Vec_3D<int>(0, 0, 2) );
}

TEST_CASE( "UNIT TEST: assign functions iterator {IT}", "[core_mesh]" )
{
    print_unit_msg("assign functions iterator {IT}");

    IT<1> it1(Vec_3D<FTYPE>(3.2, 7.8, 4.0));
    REQUIRE( it1.vec == Vec_3D<int>(3, 8, 4) );
    CHECK( it1.counter == 0);
    CHECK_FALSE( it1.iter() );

    IT<2> it2(Vec_3D<FTYPE>(3.2, 7.8, 4.0));
    REQUIRE( it2.vec == Vec_3D<int>(3, 7, 4) );
    CHECK( it2.iter() );

    IT<3> it3(Vec_3D<FTYPE>(3.2, 7.8, 4.0), 2);
    REQUIRE( it3.vec == Vec_3D<int>(0, 2, 1) );

    CHECK( it2.vec == Vec_3D<int>(3, 7, 5) );
    CHECK( it2.iter() );
    CHECK( it2.vec == Vec_3D<int>(3, 8, 4) );
    
    do{} while( it2.iter() );
    CHECK( it2.vec == Vec_3D<int>(4, 8, 5) );
}