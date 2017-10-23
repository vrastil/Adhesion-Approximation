
#include "stdafx.h"
#include "core.h"

#include "core_app.h"
#include "core_mesh.h"

using namespace std;

template<class T>
void standard_preparation(T& APP)
{
    /***************************************
    * STANDARD PREPARATION FOR INTEGRATIOM *
    ***************************************/

    /* Generating the right density distribution in k-space */	
    gen_rho_dist_k(APP.sim, &APP.app_field[0], APP.p_F);
    
    /* Computing initial potential in k-space */
    gen_pot_k(&APP.app_field[0]);
    
    /* Computing displacement in k-space */
    gen_displ_k(&APP.app_field, APP.app_field[0]);
    
    /* Computing displacement in q-space */
    printf("Computing displacement in q-space...\n");
    fftw_execute_dft_c2r_triple(APP.p_B, APP.app_field);
}

template<class T>
void init_cond_no_vel(T& APP)
{
    /*********************
    * INITIAL CONDITIONS *
    *********************/

    printf("\nSetting initial positions of particles...\n");
    set_pert_pos(APP.sim, APP.sim.b_in, APP.particles, APP.app_field);
}

template<class T>
void print_init(T& APP)
{
    /* Setting initial (binned) power spectrum, WARNING: power_aux is modified */
    APP.track.update_track_par(APP.particles);
    APP.print(APP.sim);
    APP.upd_time();
}

template<class T>
void integration(T& APP, function<void()> upd_pos)
{
    /**************
    * INTEGRATION *
    **************/
	
	while(APP.integrate())
	{
		printf("\nStarting computing step with z = %.2f (b = %.3f)\n", APP.z(), APP.b);
		printf("Updating positions of particles...\n");
		upd_pos();
        APP.track.update_track_par(APP.particles);
		if (APP.printing()) APP.print(APP.sim);
		APP.upd_time();
	}
}

int zel_app(const Sim_Param &sim)
{
    cout << "\n"
	"************************\n"
	"ZEL`DOVICH APPROXIMATION\n"
    "************************\n";

    App_Var<Particle_x> APP(sim, "ZA");
    APP.print_mem();
    standard_preparation(APP);
    init_cond_no_vel(APP);
    print_init(APP);
    auto upd_pos = [&](){set_pert_pos(sim, APP.b, APP.particles, APP.app_field);};
    integration(APP, upd_pos);
    APP.print_info();
	printf("Zel`dovich approximation ended successfully.\n");
	return APP.err;
}

int frozen_flow(const Sim_Param &sim)
{
	cout << "\n"
	"*************************\n"
	"FROZEN-FLOW APPROXIMATION\n"
	"*************************\n";

    App_Var<Particle_x> APP(sim, "FF");
    APP.print_mem();
    standard_preparation(APP);
    init_cond_no_vel(APP);
    print_init(APP);
    auto upd_pos = [&](){upd_pos_first_order(sim, APP.db, APP.particles, APP.app_field);};
    integration(APP, upd_pos);
    APP.print_info();
    printf("Frozen-flow approximation ended successfully.\n");
    return APP.err;
}