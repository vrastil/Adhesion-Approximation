
#include "stdafx.h"
#include "core.h"
#include "core_out.h"
#include "core_app.h"
#include "core_mesh.h"

using namespace std;
const double PI = acos(-1.);

int frozen_flow(const Sim_Param &sim)
{
	cout << "\n"
	"*************************\n"
	"FROZEN-FLOW APPROXIMATION\n"
	"*************************\n";
	
	string out_dir_app = sim.out_dir + "FF_run/" + currentDateTime() + "/";
	work_dir_over(out_dir_app);
	sim.print_info(out_dir_app);
    
	/** ALLOCATION OF MEMORY + FFTW PREPARATION **/
	App_Var APP(sim, "_FF_");
	printf("Initialization completed...\n");
	
	/** STANDARD PREPARATION FOR INTEGRATIOM **/
	
	/* Generating the right density distribution in k-space */	
	gen_rho_dist_k(sim, &APP.app_field[0], APP.p_F);
	
	/* Computing initial power spectrum in k-space */
	pwr_spec_k(sim, APP.app_field[0], &APP.power_aux);
	gen_pow_spec_binned(sim, APP.power_aux, &APP.pwr_spec_binned_0);
	print_pow_spec(APP.pwr_spec_binned_0, out_dir_app,  APP.z_suffix_const + "init");
	
	/* Computing initial potential in k-space */
	gen_pot_k(&APP.app_field[0]);
	
	/* Computing displacement in k-space */
	gen_displ_k(&APP.app_field, APP.app_field[0]);
	
	/* Computing displacement in q-space */
	printf("Computing displacement in q-space...\n");
	fftw_execute_dft_c2r_triple(APP.p_B, APP.app_field);
	
	/* Setting initial positions of particles */
    printf("Setting initial positions of particles...\n");
	set_unpert_pos(sim, APP.particles);
	
	/** INTEGRATION **/
	
	while(APP.integrate())
	{
		printf("\nStarting computing step with z = %.2f (b = %.3f)\n", APP.z(), APP.b);
		
		/* Updating positions of particles... */
		printf("Updating positions of particles...\n");
		upd_pos_first_order(sim, APP.db, APP.particles, APP.app_field);
		
		if (APP.printing()) APP.print(sim, out_dir_app);
		APP.upd_time();
	}
	print_suppression(APP.supp, sim, out_dir_app);
		
	printf("Frozen-flow approximation ended successfully.\n");
	return APP.err;
}