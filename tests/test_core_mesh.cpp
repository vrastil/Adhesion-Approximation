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

    Vec_3D<int> pos3i(0, -10, 10);
    Vec_3D<int> per(4, 5, 3);
    get_per(pos3i, per);
    CHECK( pos3i == Vec_3D<int>(0, 0, 1) );

    Vec_3D<int> pos4i(0, -10, 10);
    get_per(pos4i, 5, 10, 4);
    CHECK( pos4i == Vec_3D<int>(0, 0, 2) );
}

TEST_CASE( "UNIT TEST: assign functions iterator {IT}", "[core_mesh]" )
{
    IT it0(Vec_3D<double>(3.2, 7.8, 4.0), 0);
    CHECK( it0.counter == 0);
    CHECK( it0.points == 1);
    CHECK( it0.max_counter == 1);
    REQUIRE( it0.vec == Vec_3D<int>(3, 8, 4) );
    CHECK_FALSE( it0.iter() );

    IT it1(Vec_3D<double>(3.2, 7.8, 4.0), 1);
    CHECK( it1.counter == 0);
    CHECK( it1.points == 2);
    CHECK( it1.max_counter == 8);
    REQUIRE( it1.vec == Vec_3D<int>(3, 7, 4) );

    CHECK( it1.iter() );
    CHECK( it1.vec == Vec_3D<int>(3, 7, 5) );
    CHECK( it1.iter() );
    CHECK( it1.vec == Vec_3D<int>(3, 8, 4) );
    
    do{} while( it1.iter() );
    CHECK( it1.vec == Vec_3D<int>(4, 8, 5) );
    // CHECK( it1.vec[0] == 4 );
    // CHECK( it1.vec[1] == 8 );
    // CHECK( it1.vec[2] == 5 );
}