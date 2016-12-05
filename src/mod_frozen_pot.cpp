
#include "stdafx.h"
#include "mesh.h"

using namespace std;
const double PI = acos(-1.);

static void leapfrog_int(vector<c_Part_v>* particles, const vector<c_Mesh<double>> &force_field,
	double db, double b_half, int order, c_Pool* pool);

int mod_frozen_potential(const c_Sim_Param &sim){
	cout << "\n"
	"***************************************\n"
	"MODIFIED FROZEN-POTENTIAL APPROXIMATION\n"
	"***************************************\n";

	/** ALLOCATION OF MEMORY + FFTW PREPARATION **/
	c_App_Var_v APP(sim, "_FP_mod_");
	
	if (APP.err){
			printf("Errors during initializing!\n");
			return APP.err;
	}
	
	/** STANDARD PREPARATION FOR INTEGRATIOM **/
	
	/* Generating the right density distribution in k-space */
	gen_rho_dist_k(sim, &APP.app_field[0], APP.p_F, &APP.pool);
	
	/* Computing initial power spectrum in k-space */
	pwr_spec_k(sim, APP.app_field[0], &APP.power_aux, &APP.pool);
	gen_pow_spec_binned(sim, APP.power_aux, &APP.pwr_spec_binned_0);

	/* Computing initial potential in k-space */
	gen_pot_k(&APP.app_field[0], &APP.pool);
	
	/* Computing displacement in k-space */
	gen_displ_k(&APP.app_field, APP.app_field[0], &APP.pool);
	
	/* Computing displacement in q-space */
	printf("Computing displacement in q-space...\n");
	fftw_execute_dft_r2c_triple(APP.p_B, &APP.app_field, &APP.pool);
	
	/* Setting initial positions and velocitites of particles */
	APP.set_unpert_pos_w_vel(sim, APP.app_field);
	
	/** INTEGRATION **/
	
	while(APP.integrate())
	{
		printf("Starting computing step with z = %.2f (b = %.2f)\n", APP.z(), APP.b);
		leapfrog_int(&APP.particles, APP.app_field, APP.db, APP.b - APP.db/2., sim.order, &APP.pool);
		
		if (APP.printing())
		{
			/* Printing positions */
			
			/* Printing rho map */
			get_rho_from_par(APP.particles, &APP.power_aux, sim.order, &APP.pool);
			
			/* Printing power spectrum */
			pwr_spec(sim, &APP.power_aux, &APP.power_aux, APP.p_F, &APP.pool);
			gen_pow_spec_binned(sim, APP.power_aux, &APP.pwr_spec_binned);
		}
		
		APP.upd_time();
	}
	
	return APP.err;
}


static void leapfrog_int_th(int i_min, int i_max, vector<c_Part_v>* particles,
	const vector<c_Mesh<double>> &force_field, double db, double b_half, int order)
{
	double v, f;
	for(int i=i_min; i< i_max; i++)
	{
		(*particles)[i].drift((*particles)[i].velocity, db / 2.);
		for(int k=0; k<3;k++)
		{
			f = 0;
			force_field[k].assign_from((*particles)[i].position, &f, order);
			v = (*particles)[i].velocity[k];
			f = -3/(2.*b_half)*(v-f);
			
			v += f*db; // <- KICK at middle-step
			(*particles)[i].position[k]+=v*db/2.; // <- DRIFT
			(*particles)[i].velocity[k] = v;
		}
	}
}

static void leapfrog_int(vector<c_Part_v>* particles, const vector<c_Mesh<double>> &force_field,
	double db, double b_half, int order, c_Pool* pool)
{
	printf("Updating positions and velocitites of particles...\n");
	auto tmp_func = bind(leapfrog_int_th, placeholders::_1, placeholders::_2, particles, force_field, db, b_half, order);
	pool->add_task(0, particles->size(), tmp_func);
}