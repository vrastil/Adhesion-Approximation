
#include "stdafx.h"
#include "core.h"
#include "core_out.h"
#include "core_app.h"
#include "core_mesh.h"

using namespace std;

int zel_app(const Sim_Param &sim)
{
	cout << "\n"
	"************************\n"
	"ZEL`DOVICH APPROXIMATION\n"
	"************************\n";
	
    string out_dir_app = std_out_dir("ZA_run/", sim);
	work_dir_over(out_dir_app);
    
    /******************************************
    * ALLOCATION OF MEMORY + FFTW PREPARATION *
    ******************************************/

	App_Var APP(sim, "_ZA_");
	printf("Initialization completed...\n");
	
    /***************************************
    * STANDARD PREPARATION FOR INTEGRATIOM *
    ***************************************/
    
	/* Generating the right density distribution in k-space */	
    gen_rho_dist_k(sim, &APP.app_field[0], APP.p_F);
    
	/* Computing initial potential in k-space */
	gen_pot_k(&APP.app_field[0]);
	
	/* Computing displacement in k-space */
	gen_displ_k(&APP.app_field, APP.app_field[0]);
	
	/* Computing displacement in q-space */
	printf("Computing displacement in q-space...\n");
    fftw_execute_dft_c2r_triple(APP.p_B, APP.app_field);
    
    /*********************
    * INITIAL CONDITIONS *
    *********************/

    printf("\nSetting initial positions of particles...\n");
    set_pert_pos(sim, sim.b_in, APP.particles, APP.app_field);

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
		
		/* Computing displaced positions of particles */
		printf("Computing displaced positions of particles...\n");
		set_pert_pos(sim, APP.b, APP.particles, APP.app_field);
        
        APP.track.update_track_par(APP.particles);
		if (APP.printing()) APP.print(sim, out_dir_app);
		APP.upd_time();
	}
        
    sim.print_info(out_dir_app, "ZA");
	printf("Zel`dovich approximation ended successfully.\n");
	return APP.err;
}