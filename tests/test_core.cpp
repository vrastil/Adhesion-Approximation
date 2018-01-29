#include <catch.hpp>

TEST_CASE( "UNIT TEST: sign function {sgn<T>}", "[core]" )
{
    CHECK( sgn(0) == 0 );
    CHECK( sgn(10) == 1 );
    CHECK( sgn(-5) == -1 );
    CHECK( sgn(0.) == 0 );
    CHECK( sgn(0.34534f) == 1 );
    CHECK( sgn(-1.34534E16) == -1 );
}

TEST_CASE( "UNIT TEST: vector class {Vec_3D<T>}", "[core]" )
{
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

TEST_CASE( "UNIT TEST: mesh class {Mesh_base<T>}", "[core]" )
{
    // dimension
    Mesh_base<double> mesh(8, 16, 20);
    CHECK( mesh.N1 == 8 );
    CHECK( mesh.N2 == 16 );
    CHECK( mesh.N3 == 20 );
    CHECK( mesh.length == 2560 );
    mesh.assign(-0.5);
    CHECK( mesh.length == 2560 );
    CHECK( mesh.N2 == 16 );
    CHECK( mesh[5] == -0.5 );

    // writing, reading
    mesh[5] = 3.14;
    REQUIRE( mesh[5] == 3.14 );
    CHECK( mesh(0,0,5) == 3.14 );
    CHECK( *(mesh.real()+5) == 3.14 );

    mesh(3,4,1) = 2.71;
    REQUIRE( mesh[1041] == 2.71 );
    CHECK( mesh(3,4,1) == 2.71 );
    CHECK( mesh(52,1) == 2.71 );
    Vec_3D<int> pos(3,4,1);
    CHECK( mesh(pos) == 2.71 );
    
    // copy constructor
    Mesh_base<double> mesh2(mesh);
    CHECK( mesh2[0] == -0.5 );
    CHECK( mesh2[5] == 3.14 );
    CHECK( mesh2[1041] == 2.71 );

    // assign operator
    Mesh_base<double> mesh3(2,5,7);
    mesh3 = mesh;
    CHECK( mesh3.N1 == 8 );
    CHECK( mesh3.N2 == 16 );
    CHECK( mesh3.N3 == 20 );
    CHECK( mesh3.length == 2560 );
    CHECK( mesh3[0] == -0.5 );
    CHECK( mesh3[5] == 3.14 );
    CHECK( mesh3[1041] == 2.71 );

    // operations
    mesh*=2.;
    CHECK( mesh[0] == Approx(-1.) );
    CHECK( mesh[5] == Approx(6.28) );
    mesh+=1.;
    CHECK( mesh[0] == Approx(0.) );
    CHECK( mesh[5] == Approx(7.28) );
    mesh/=7.28;
    CHECK( mesh[0] == Approx(0.) );
    CHECK( mesh[5] == Approx(1.) );
}

TEST_CASE( "UNIT TEST: mesh class {Mesh}", "[core]" )
{
    // dimension
    Mesh mesh_c(8);
    mesh_c.assign(0.);
    CHECK( mesh_c.N == 8 );
    CHECK( mesh_c.N1 == 8 );
    CHECK( mesh_c.N2 == 8 );
    CHECK( mesh_c.N3 == 10 );

    // writing, reading
    mesh_c[90] = 3.14;
    REQUIRE( mesh_c[90] == (FTYPE)3.14 );
    CHECK( mesh_c(1,1,0) == (FTYPE)3.14 );

    Vec_3D<int> pos(1,1,0);
    CHECK( mesh_c(pos) == (FTYPE)3.14 );
    pos = Vec_3D<int>(9,-7,8);
    CHECK( mesh_c(pos) == (FTYPE)3.14 );
    pos = Vec_3D<int>(-6,10,0);
    mesh_c(pos) = (FTYPE)2.5;
    CHECK( mesh_c(2,2,0) == (FTYPE)2.5 );
    CHECK( mesh_c[180] == (FTYPE)2.5 );

    // copy constructor
    Mesh mesh2_c(mesh_c);
    CHECK( mesh2_c[90] == (FTYPE)3.14 );
    CHECK( mesh2_c[180] == (FTYPE)2.5 );

    // assign operator
    Mesh mesh3_c(14);
    mesh3_c = mesh_c;
    CHECK( mesh3_c.N1 == 8 );
    CHECK( mesh3_c.N2 == 8 );
    CHECK( mesh3_c.N3 == 10 );
    CHECK( mesh3_c.N == 8 );
    CHECK( mesh3_c[90] == (FTYPE)3.14 );
    CHECK( mesh3_c[180] == (FTYPE)2.5 );

    // operations
    mesh3_c*=2.;
    CHECK( mesh3_c[90] == Approx(6.28) );
    CHECK( mesh3_c[180] == Approx(5) );
    mesh3_c-=5.;
    CHECK( mesh3_c[90] == Approx(1.28) );
    CHECK( mesh3_c[180] == Approx(0) );
    mesh3_c/=1.28;
    CHECK( mesh3_c[90] == Approx(1.) );
    CHECK( mesh3_c[180] == Approx(0) );
}
    
TEST_CASE( "UNIT TEST: particle class {Particle_x}", "[core]" )
{
    Vec_3D<FTYPE> position(0., -3.14, 4E5);
    Particle_x par1(position);

    CHECK( par1.position[0] == 0. );
    CHECK( par1.position[1] == (FTYPE)-3.14 );
    CHECK( par1.position[2] == (FTYPE)4E5 );
    CHECK( par1[1] == (FTYPE)-3.14 );

    par1.position[0] = par1.position[1]*2+6.28;
    CHECK( par1[0] == Approx(0) );
    
    // copy constructor
    Particle_x par3(par1);
    CHECK( par3[0] == Approx(0) );
    CHECK( par3[1] == (FTYPE)-3.14 );
    CHECK( par3[2] == (FTYPE)4E5 );

    // assign operator
    par1 = par3;
    CHECK( par1[0] == Approx(0) );
    CHECK( par1[1] == (FTYPE)-3.14 );
    CHECK( par1[2] == (FTYPE)4E5 );
}

TEST_CASE( "UNIT TEST: particle class {Particle_v}", "[core]" )
{    
    Vec_3D<FTYPE> position(0., -3.14, 4E5);
    Vec_3D<FTYPE> velocity(2.3E-6, -4.56E-7, 6.87903E-6);
    Particle_v par1(position, velocity);

    CHECK( par1[0] == (FTYPE)0. );
    CHECK( par1[1] == (FTYPE)-3.14 );
    CHECK( par1[2] == (FTYPE)4E5 );
    CHECK( par1(0) == (FTYPE)2.3E-6 );
    CHECK( par1(1) == (FTYPE)-4.56E-7 );
    CHECK( par1(2) == (FTYPE)6.87903E-6 );
    
    // copy constructor
    Particle_v par3(par1);
    CHECK( par3[1] == (FTYPE)-3.14 );
    CHECK( par3(1) == (FTYPE)-4.56E-7 );

    // assign operator
    par1 = par3;
    CHECK( par1[1] == (FTYPE)-3.14 );
    CHECK( par1(1) == (FTYPE)-4.56E-7 );
}

TEST_CASE( "UNIT TEST: tracking class {Tracking}", "[core]" )
{
    // 2, 12, 2
    Tracking track(2, 24);

    CHECK( track.num_track_par == 4 );
    CHECK( track.par_ids.size() == 4 );
    CHECK( track.par_ids[0] == 1442 );
}
