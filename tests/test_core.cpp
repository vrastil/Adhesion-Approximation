#include <catch.hpp>

TEST_CASE( "UNIT TEST: sign function {sgn<T>}", "[core]" )
{
    CHECK( sgn<int>(0) == 0 );
    CHECK( sgn<int>(10) == 1 );
    CHECK( sgn<int>(-5) == -1 );
    CHECK( sgn<double>(0.) == 0 );
    CHECK( sgn<double>(0.34534535) == 1 );
    CHECK( sgn<double>(-1.34534E16) == -1 );

}

TEST_CASE( "UNIT TEST: vector class {Vec_3D<T>}", "[core]" )
{
    Vec_3D<double> vec_d(sqrt(2.), -sqrt(2.), sqrt(5.));
    Vec_3D<int> vec_i(3,0,-4);
    CHECK( vec_d.norm() == Approx(3.) );
    CHECK( vec_i.norm() == Approx(5) );

    vec_d.assign(3., 5.56, -1.345E1);
    REQUIRE( vec_d.x == 3. );
    REQUIRE( vec_d[2] == -13.45 );
    REQUIRE( Vec_3D<int>(vec_d).y == 5 );

    vec_i.assign(0, 0, 0);
    vec_i+=Vec_3D<int>(2, 3, -4);
    REQUIRE ( (vec_i.x + vec_i[1] + vec_i.z) == Approx(1) );

    Vec_3D<double> vec_d2 = (Vec_3D<double>)vec_i + Vec_3D<double>(1., 1.5, -3.5)*2.;
    CHECK( (vec_d2[0] + vec_d2[1] + vec_d2[2]) == Approx(-1.) );

    //Invalid acces
    CHECK( vec_d[-6] == 3. );
    CHECK( vec_d[7] == -13.45 );
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
    REQUIRE( mesh_c[90] == 3.14 );
    CHECK( mesh_c(1,1,0) == 3.14 );

    Vec_3D<int> pos(1,1,0);
    CHECK( mesh_c(pos) == 3.14 );
    pos = Vec_3D<int>(9,-7,8);
    CHECK( mesh_c(pos) == 3.14 );
    pos = Vec_3D<int>(-6,10,0);
    mesh_c(pos) = 2.5;
    CHECK( mesh_c(2,2,0) == 2.5 );
    CHECK( mesh_c[180] == 2.5 );

    // copy constructor
    Mesh mesh2_c(mesh_c);
    CHECK( mesh2_c[90] == 3.14 );
    CHECK( mesh2_c[180] == 2.5 );

    // assign operator
    Mesh mesh3_c(14);
    mesh3_c = mesh_c;
    CHECK( mesh3_c.N1 == 8 );
    CHECK( mesh3_c.N2 == 8 );
    CHECK( mesh3_c.N3 == 10 );
    CHECK( mesh3_c.N == 8 );
    CHECK( mesh3_c[90] == 3.14 );
    CHECK( mesh3_c[180] == 2.5 );

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
    Particle_x par1(0., -3.14, 4E5);
    Vec_3D<double> position(0., -3.14, 4E5);
    Particle_x par2(position);

    CHECK( par1.position.x == 0. );
    REQUIRE( par1.position.y == -3.14 );
    CHECK( par1.position.z == 4E5 );
    CHECK( par1[1] == -3.14 );
    CHECK( par1.position.x == par2[0] );
    CHECK( par1.position.y == par2[1] );
    CHECK( par1.position.z == par2[2] );

    par1.position.x = par1.position.y*2+6.28;
    CHECK( par1[0] == Approx(0) );
    
    // copy constructor
    Particle_x par3(par1);
    CHECK( par3[0] == Approx(0) );
    CHECK( par3[1] == -3.14 );
    CHECK( par3[2] == 4E5 );

    // assign operator
    par2 = par3;
    CHECK( par2[0] == Approx(0) );
    CHECK( par2[1] == -3.14 );
    CHECK( par2[2] == 4E5 );
}

TEST_CASE( "UNIT TEST: particle class {Particle_v}", "[core]" )
{    
    Particle_v par1(0., -3.14, 4E5, 2.3E-6, -4.56E-7, 6.87903E-6);
    Vec_3D<double> position(0., -3.14, 4E5);
    Vec_3D<double> velocity(2.3E-6, -4.56E-7, 6.87903E-6);
    Particle_v par2(position, velocity);

    CHECK( par1[0] == 0. );
    CHECK( par1[1] == -3.14 );
    CHECK( par1[2] == 4E5 );
    CHECK( par1(0) == 2.3E-6 );
    CHECK( par1(1) == -4.56E-7 );
    CHECK( par1(2) == 6.87903E-6 );

    CHECK( par1.position.x == par2[0] );
    CHECK( par1.position.y == par2[1] );
    CHECK( par1.position.z == par2[2] );
    CHECK( par1.velocity.x == par2(0) );
    CHECK( par1.velocity.y == par2(1) );
    CHECK( par1.velocity.z == par2(2) );
    
    // copy constructor
    Particle_v par3(par1);
    CHECK( par3[1] == -3.14 );
    CHECK( par3(1) == -4.56E-7 );

    // assign operator
    par2 = par3;
    CHECK( par2[1] == -3.14 );
    CHECK( par2(1) == -4.56E-7 );
}

TEST_CASE( "UNIT TEST: power spectrum class {e_power_spec}", "[core]" ){
    e_power_spec ps0 = power_law_T;
    e_power_spec ps1 = power_law;
    e_power_spec ps2 = flat;
    e_power_spec ps3 = single;

    CHECK( ps0 == 0 );
    CHECK( ps1 == 1 );
    CHECK( ps2 == 2 );
    CHECK( ps3 == 3 );
}

TEST_CASE( "UNIT TEST: tracking class {Tracking}", "[core]" )
{
    // 2, 12, 2
    Tracking track(2, 24);

    CHECK( track.num_track_par == 4 );
    CHECK( track.par_ids.size() == 4 );
    CHECK( track.par_ids[0] == 1442 );
}
