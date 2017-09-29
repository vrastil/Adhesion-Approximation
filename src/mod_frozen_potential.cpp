
#include "stdafx.h"
#include "core.h"
#include "core_out.h"
#include "core_app.h"
#include "core_mesh.h"

using namespace std;
const double PI = acos(-1.);

int mod_frozen_potential(const Sim_Param &sim)
{
	cout << "\n"
	"***************************************\n"
	"MODIFIED FROZEN-POTENTIAL APPROXIMATION\n"
	"***************************************\n";
	
    string out_dir_app = std_out_dir("FP_pp_run/", sim);
	work_dir_over(out_dir_app);
    sim.print_info(out_dir_app, "FP_pp");
	
	/******************************************
    * ALLOCATION OF MEMORY + FFTW PREPARATION *
    ******************************************/

	App_Var_FP_mod APP(sim, "_FP_pp_");
	printf("Initialization completed...\n");

	/***************************************
    * STANDARD PREPARATION FOR INTEGRATIOM *
    ***************************************/

	/* Generating the right density distribution in k-space */	
	gen_rho_dist_k(sim, &APP.app_field[0], APP.p_F);
	
	/* Computing initial potential in k-space */
	gen_pot_k(APP.app_field[0], &APP.power_aux);
	
	/* Computing displacement in k-space */
	gen_displ_k(&APP.app_field, APP.power_aux);
	
	/* Computing displacement in q-space */
	printf("Computing displacement in q-space...\n");
	fftw_execute_dft_c2r_triple(APP.p_B, APP.app_field);
	
	/*********************
    * INITIAL CONDITIONS *
    *********************/

    printf("\nSetting initial positions and velocitis of particles...\n");
	set_pert_pos_w_vel(sim, sim.b_in, APP.particles, APP.app_field);

    /*******************************************************
    * PREPARATION FOR INTEGRATIOM WITH S2 SHAPED PARTICLES *
    *******************************************************/

	/* Computing displacement in k-space with S2 shaped particles */
	gen_displ_k_S2(&APP.app_field, APP.power_aux, sim.a);
	
	/* Computing force in q-space */
	printf("Computing force for S2 shaped particles in q-space...\n");
    fftw_execute_dft_c2r_triple(APP.p_B, APP.app_field);
    
    /* Setting initial (binned) power spectrum, WARNING: power_aux is modified */
    APP.track.update_track_par(APP.particles);
	APP.print(sim, out_dir_app);
	APP.upd_time();
	
	/**************
    * INTEGRATION *
    **************/
	
	while(APP.integrate())
	{
		printf("\nStarting computing step with z = %.2f (b = %.3f)\n", APP.z(), APP.b);
		
		/* Updating positions of particles... */
		printf("Updating positions of particles...\n");
		upd_pos_second_order_w_short_force(sim, &APP.linked_list, APP.db, APP.b, APP.particles, APP.app_field);
        
        APP.track.update_track_par(APP.particles);
		if (APP.printing()) APP.print(sim, out_dir_app);
		APP.upd_time();
	}
		
	printf("Modified Frozen-potential approximation ended successfully.\n");
	return APP.err;
}