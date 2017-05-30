
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
	
	string out_dir_app = sim.out_dir + "FP_pp_run/" + currentDateTime() + "/";
	work_dir_over(out_dir_app);
    sim.print_info(out_dir_app);
	
	/** ALLOCATION OF MEMORY + FFTW PREPARATION **/
	App_Var_FP_mod APP(sim, "_FP_pp_");
	printf("Initialization completed...\n");

	/** STANDARD PREPARATION FOR INTEGRATIOM **/
	/* Generating the right density distribution in k-space */	
	gen_rho_dist_k(sim, &APP.app_field[0], APP.p_F);
	
	/* Computing initial power spectrum in k-space */
	pwr_spec_k(sim, APP.app_field[0], &APP.power_aux);
	gen_pow_spec_binned(sim, APP.power_aux, &APP.pwr_spec_binned_0);
	print_pow_spec(APP.pwr_spec_binned_0, out_dir_app,  APP.z_suffix_const + "init");
	
	/* Computing initial potential in k-space */
	gen_pot_k(APP.app_field[0], &APP.power_aux);
	
	/* Computing displacement in k-space */
	gen_displ_k(&APP.app_field, APP.power_aux);
	
	/* Computing displacement in q-space */
	printf("Computing displacement in q-space...\n");
	fftw_execute_dft_c2r_triple(APP.p_B, APP.app_field);
	
	/** INITIAL CONDITIONS **/
	/* Setting initial positions of particles */
    printf("Setting initial positions and velocitis of particles...\n");
	set_pert_pos_w_vel(sim, sim.b_in, APP.particles, APP.app_field);

	/** PREPARATION FOR INTEGRATIOM WITH S2 SHAPED PARTICLES **/
	/* Computing displacement in k-space with S2 shaped particles */
	gen_displ_k_S2(&APP.app_field, APP.power_aux, sim.a);
	
	/* Computing force in q-space */
	printf("Computing force for S2 shaped particles in q-space...\n");
	fftw_execute_dft_c2r_triple(APP.p_B, APP.app_field);
	
	APP.print(sim, out_dir_app);
	APP.upd_time();
	
	/** INTEGRATION **/
	
	while(APP.integrate())
	{
		printf("\nStarting computing step with z = %.2f (b = %.3f)\n", APP.z(), APP.b);
		
		/* Updating positions of particles... */
		printf("Updating positions of particles...\n");
		upd_pos_second_order_w_short_force(sim, &APP.linked_list, APP.db, APP.b, APP.particles, APP.app_field);
		
		if (APP.printing()) APP.print(sim, out_dir_app);
		APP.upd_time();
	}
	print_suppression(APP.supp, sim, out_dir_app);
		
	printf("Modified Frozen-potential approximation ended successfully.\n");
	return APP.err;
}