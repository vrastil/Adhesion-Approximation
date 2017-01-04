
#include "stdafx.h"
#include "core.h"


void set_unpert_pos(const Sim_Param &sim, Particle_x* particles)
{
	int x, y, z;
	int par_per_dim = sim.mesh_num / sim.Ng;
	
	#pragma omp parallel for private(x,y,z)
	for(int i=0; i< sim.par_num; i++)
	{
		x = i / (par_per_dim * par_per_dim);
		y = (i / par_per_dim) % par_per_dim;
		z = i % par_per_dim;
		
		particles[i] = Particle_x(x, y, z);
	}
}

void set_unpert_pos_w_vel(const Sim_Param &sim, Particle_v* particles, const std::vector<Mesh> &vel_field)
{
	int x, y, z;
	double vx, vy, vz;
	int par_per_dim = sim.mesh_num / sim.Ng;
	
	#pragma omp parallel for private(x,y,z,vx,vy,vz)
	for(int i=0; i< sim.par_num; i++)
	{
		x = i / (par_per_dim * par_per_dim);
		y = (i / par_per_dim) % par_per_dim;
		z = i % par_per_dim;
		vx = vel_field[0](x, y, z);
		vy = vel_field[1](x, y, z);
		vz = vel_field[2](x, y, z);
		particles[i] = Particle_v(x, y, z, vx, vy, vz);
	}
}