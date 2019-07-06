#include <catch.hpp>
#include "test.hpp"
#include "app_var.cpp"

TEST_CASE( "UNIT TEST: tracking class {Tracking}", "[core]" )
{
    print_unit_msg("tracking class {Tracking}");
    const size_t par_per_dim = 8;
    Tracking track(2, par_per_dim);
    std::vector<Particle_x<FTYPE_t>> particles;
    particles.resize(par_per_dim*par_per_dim*par_per_dim);

    CHECK( track.get_num_track_par() == 4 );
    CHECK( track.get_num_steps() == 0 );

    track.update_track_par(particles);

    CHECK( track.get_num_steps() == 1 );
}