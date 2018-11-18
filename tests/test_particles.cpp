#include <catch.hpp>
#include "test.hpp"
#include "class_particles.hpp"
#include "precision.hpp"
    
TEST_CASE( "UNIT TEST: particle class {Particle_x}", "[core]" )
{
    print_unit_msg("particle class {Particle_x}");

    Vec_3D<FTYPE_t> position(0., -3.14, 4E5);
    Particle_x<FTYPE_t> par1(position);

    CHECK( par1.position[0] == 0. );
    CHECK( par1.position[1] == (FTYPE_t)-3.14 );
    CHECK( par1.position[2] == (FTYPE_t)4E5 );
    CHECK( par1[1] == (FTYPE_t)-3.14 );

    par1.position[0] = par1.position[1]*2+6.28;
    CHECK( par1[0] == Approx(0) );
    
    // copy constructor
    Particle_x<FTYPE_t> par3(par1);
    CHECK( par3[0] == Approx(0) );
    CHECK( par3[1] == (FTYPE_t)-3.14 );
    CHECK( par3[2] == (FTYPE_t)4E5 );

    // assign operator
    par1 = par3;
    CHECK( par1[0] == Approx(0) );
    CHECK( par1[1] == (FTYPE_t)-3.14 );
    CHECK( par1[2] == (FTYPE_t)4E5 );
}

TEST_CASE( "UNIT TEST: particle class {Particle_v}", "[core]" )
{
    print_unit_msg("particle class {Particle_v}");

    Vec_3D<FTYPE_t> position(0., -3.14, 4E5);
    Vec_3D<FTYPE_t> velocity(2.3E-6, -4.56E-7, 6.87903E-6);
    Particle_v<FTYPE_t> par1(position, velocity);

    CHECK( par1[0] == (FTYPE_t)0. );
    CHECK( par1[1] == (FTYPE_t)-3.14 );
    CHECK( par1[2] == (FTYPE_t)4E5 );
    CHECK( par1(0) == (FTYPE_t)2.3E-6 );
    CHECK( par1(1) == (FTYPE_t)-4.56E-7 );
    CHECK( par1(2) == (FTYPE_t)6.87903E-6 );
    
    // copy constructor
    Particle_v<FTYPE_t> par3(par1);
    CHECK( par3[1] == (FTYPE_t)-3.14 );
    CHECK( par3(1) == (FTYPE_t)-4.56E-7 );

    // assign operator
    par1 = par3;
    CHECK( par1[1] == (FTYPE_t)-3.14 );
    CHECK( par1(1) == (FTYPE_t)-4.56E-7 );
}