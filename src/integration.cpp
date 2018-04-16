#include "core_mesh.h"
#include "core_power.h"
#include "integration.hpp"
#include "params.hpp"

void stream_step(const FTYPE_t da, std::vector<Particle_v<FTYPE_t>>& particles)
{
    const unsigned Np = particles.size();
    #pragma omp parallel for
	for (unsigned i = 0; i < Np; i++)
	{
        particles[i].position += particles[i].velocity*da;
    }
}

void stream_kick_stream(const FTYPE_t da, std::vector<Particle_v<FTYPE_t>>& particles, std::function<void()> kick_step, unsigned per)
{// general Leapfrog method: Stream-Kick-Stream & ensure periodicity
    stream_step(da/2, particles);
    kick_step();
    stream_step(da/2, particles);
    get_per(particles, per);
}

void kick_step_no_momentum(const Cosmo_Param &cosmo, const FTYPE_t a, std::vector<Particle_v<FTYPE_t>>& particles, const std::vector< Mesh> &vel_field)
{
    // no memory of previus velocity, 1st order ODE
    const unsigned Np = particles.size();
    Vec_3D<FTYPE_t> vel;
    const FTYPE_t dDda = growth_change(a, cosmo); // dD / da
    
    #pragma omp parallel for private(vel)
    for (unsigned i = 0; i < Np; i++)
	{
        vel.fill(0.);
        assign_from(vel_field, particles[i].position, vel);
        particles[i].velocity = vel*dDda;
    }
}

void kick_step_w_momentum(const Cosmo_Param &cosmo, const FTYPE_t a, const FTYPE_t da, std::vector<Particle_v<FTYPE_t>>& particles, const std::vector< Mesh> &force_field)
{
    // classical 2nd order ODE
    const unsigned Np = particles.size();
    Vec_3D<FTYPE_t> force;
    const FTYPE_t D = growth_factor(a, cosmo);
    const FTYPE_t OL = cosmo.Omega_L()*pow(a,3);
    const FTYPE_t Om = cosmo.Omega_m;
    // -3/2a represents usual EOM, the rest are LCDM corrections
    const FTYPE_t f1 = 3/(2*a)*(Om+2*OL)/(Om+OL);
    const FTYPE_t f2 = 3/(2*a)*Om/(Om+OL)*D/a;
    
    #pragma omp parallel for private(force)
    for (unsigned i = 0; i < Np; i++)
	{
        force.fill(0.);
        assign_from(force_field, particles[i].position, force);
        force = force*f2 - particles[i].velocity*f1;		
        particles[i].velocity += force*da;
    }
}