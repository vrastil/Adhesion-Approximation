#include <catch.hpp>

TEST_CASE( "UNIT TEST: vector class {Vec_3D<T>}", "[core]" )
{
    print_unit_msg("vector class {Vec_3D<T>}");

    Vec_3D<double> vec_d(sqrt(2.), -sqrt(2.), sqrt(5.));
    Vec_3D<float> vec_f(sqrt(2.f), -sqrt(2.f), sqrt(5.f));
    Vec_3D<int> vec_i(3,0,-4);
    CHECK( vec_d.norm() == Approx(3.) );
    CHECK( vec_f.norm() == Approx(3.) );
    CHECK( vec_i.norm() == Approx(5) );

    REQUIRE( vec_f[0] == sqrt(2.f) );
    REQUIRE( vec_f[2] == Approx(sqrt(5.)) );

    vec_d.fill(-1.345E1);
    REQUIRE( vec_d[0] == -1.345E1 );
    REQUIRE( vec_d[2] == -13.45 );

    vec_i.fill(0);
    vec_i+=Vec_3D<int>(2, 3, -4);
    REQUIRE ( (vec_i[0] + vec_i[1] + vec_i[2]) == Approx(1) );
    Vec_3D<double> vec_d2 = (Vec_3D<double>)vec_i + Vec_3D<double>(1., 1.5, -3.5)*2.;
    REQUIRE( vec_d2[0] == Approx(4.) );
    REQUIRE( vec_d2[1] == Approx(6.) );
    REQUIRE( vec_d2[2] == Approx(-11.) );

    double sumd = 0;
    double sumi = 0;

    for (double val : vec_d2) sumd += val;
    for (int val : vec_i) sumi += val;

    CHECK( sumd == Approx(-1.) );
    CHECK( sumi == 1 );

    CHECK( Vec_3D<int>(4, -3, 8) ==  Vec_3D<int>(4, -3, 8) );
    CHECK_FALSE( Vec_3D<int>(4, -3, 8) !=  Vec_3D<int>(4, -3, 8) );
    CHECK( Vec_3D<int>(4, -3, 8) <=  Vec_3D<int>(6, 0, 8) );
    CHECK_FALSE( Vec_3D<int>(4, -3, 8) <  Vec_3D<int>(2, 0, 10) );
    CHECK_FALSE( Vec_3D<int>(4, -3, 8) >  Vec_3D<int>(4, 0, 8) );
}