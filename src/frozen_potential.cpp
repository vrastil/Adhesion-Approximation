
#include "stdafx.h"
#include "core.h"
#include "core_out.h"
#include "core_app.h"
#include "core_mesh.h"

using namespace std;
const double PI = acos(-1.);

int frozen_potential(const Sim_Param &sim)
{
	cout << "\n"
	"******************************\n"
	"FROZEN-POTENTIAL APPROXIMATION\n"
	"******************************\n";
	
	string out_dir_app = sim.out_dir + "FP_run/";
	work_dir_over(out_dir_app);
	
	/** ALLOCATION OF MEMORY + FFTW PREPARATION **/
	App_Var_v APP(sim, "_FP_");
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
    printf("Setting initial positions and velocitis of particles...\n");
	set_unpert_pos_w_vel(sim, APP.particles, APP.app_field);

	/** INTEGRATION **/
	
	while(APP.integrate())
	{
		printf("\nStarting computing step with z = %.2f (b = %.3f)\n", APP.z(), APP.b);
		
		/* Updating positions of particles... */
		printf("Updating positions of particles...\n");
		upd_pos_second_order(sim, APP.db, APP.b, APP.particles, APP.app_field);
		
		if (APP.printing())
		{
			/* Printing positions */
			print_par_pos_cut_small(APP.particles, sim, out_dir_app, APP.z_suffix());
			APP.track.update_track_par(APP.particles);
			print_track_par(APP.track, sim, out_dir_app, APP.z_suffix());

			/* Printing density */
			get_rho_from_par(APP.particles, &APP.power_aux, sim);
			gen_dens_binned(APP.power_aux, APP.dens_binned, sim);
			print_rho_map(APP.power_aux, sim, out_dir_app, APP.z_suffix());
			print_dens_bin(APP.dens_binned, sim.mesh_num, out_dir_app, APP.z_suffix());
			
			/* Printing power spectrum */
			fftw_execute_dft_r2c(APP.p_F, APP.power_aux);
			pwr_spec_k(sim, APP.power_aux, &APP.power_aux);
			gen_pow_spec_binned(sim, APP.power_aux, &APP.pwr_spec_binned);
			print_pow_spec(APP.pwr_spec_binned, out_dir_app, APP.z_suffix());
			print_pow_spec_diff(APP.pwr_spec_binned, APP.pwr_spec_binned_0, APP.b, out_dir_app, APP.z_suffix());
			
			APP.upd_supp();
		}
		APP.upd_time();
	}
	print_suppression(APP.supp, sim, out_dir_app);
		
	printf("Frozen-potential approximation ended successfully.\n");
	return APP.err;
}