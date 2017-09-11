#include <catch.hpp>

TEST_CASE( "UNIT TEST: sign function {sgn<T>}", "[core]" ){

    CHECK( sgn<int>(0) == 0 );
    CHECK( sgn<int>(10) == 1 );
    CHECK( sgn<int>(-5) == -1 );
    CHECK( sgn<double>(0.) == 0 );
    CHECK( sgn<double>(0.34534535) == 1 );
    CHECK( sgn<double>(-1.34534E16) == -1 );

}

TEST_CASE( "UNIT TEST: vector class {Vec_3D<T>}", "[core]" ){

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
}

TEST_CASE( "UNIT TEST: mesh class {Mesh_base<T>}", "[core]" ){

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
    CHECK( mesh[0] == -1. );
    CHECK( mesh[5] == 6.28 );
    mesh+=1.;
    CHECK( mesh[0] == 0. );
    CHECK( mesh[5] == 7.28 );
    mesh/=7.28;
    CHECK( mesh[0] == 0. );
    CHECK( mesh[5] == 1. );
}