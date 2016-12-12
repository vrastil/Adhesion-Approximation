
#include "stdafx.h"
#include <fftw3.h>
#include "output.h"
#include "grid_fce.h"

using namespace std;
const double PI = acos(-1.);

void mod_frozen_flow(int mesh_num, int Ng, int L, t_power power_spectrum, double* parameters, double z_in, double z_out, int nt, string out_dir){
	cout << "\n"
	"***********************************\n"
	" MODIFIED FROZEN-FLOW APPROXIMATION\n"
	"***********************************\n";

	/* ALLOCATION OF MEMORY + FFTW PREPARATION */
	int par_num = mesh_num / Ng;
	int bin_num = 100;
	void** arrays = alloc_mod_frozen_flow(mesh_num, par_num, bin_num, nt);
	if (arrays == NULL) {printf("Exiting...\n"); return;}
	double** vel_field = (double**)arrays[0];
	double** par_pos = (double**)arrays[1];
	double* power_aux = (double*)arrays[2];
	fftw_complex* pwr_spec_binned = (fftw_complex*)arrays[3];
	fftw_complex* pwr_spec_binned_0 = (fftw_complex*)arrays[4];
	fftw_plan p_F = (fftw_plan)arrays[5];
	fftw_plan p_B = (fftw_plan)arrays[6];
	
	const double b_in = 1./(z_in + 1);
	const double b_out = 1./(z_out + 1);
	double b = b_in;
	double db = b_in;
	double z;
	int i = 0;

	string z_suffix_const = "_mFF_";
	string z_suffix;
	stringstream z_suffix_num;	
	/* END OF ALLOCATION */
	
	/* Generating the right density distribution in k-space */
	gen_rho_dist_k(mesh_num, L, vel_field[0], p_F, power_spectrum, parameters, 0., nt); // displ_vec[0] used as rho(k)
	
	/* Computing initial power spectrum in k-space */
	pwr_spec_0_nt(mesh_num, L, vel_field[0], reinterpret_cast<fftw_complex*>(power_aux), nt, 1);
	gen_pow_spec_binned(mesh_num, reinterpret_cast<fftw_complex*>(power_aux), pwr_spec_binned_0,  2.*PI/L, 2.*PI*mesh_num/L, bin_num);
	
	/* Computing initial potential in k-space */
	printf("Computing initial potential in k-space...\n");
	gen_pot_k_nt(mesh_num, vel_field[0], nt);
	
	/* Computing displacement in k-space */
	printf("Computing displacement in k-space...\n");	
	gen_displ_k_w_nt(mesh_num, vel_field, vel_field[0], nt);
	
	/* Computing displacement in q-space */
	printf("Computing displacement in q-space...\n");
	
	fftw_execute_dft_c2r_triple(p_B, reinterpret_cast<fftw_complex**>(vel_field), vel_field);
	normalize_FFT_BACKWARD_nt(mesh_num, vel_field, nt);
	
	/* Setting initial positions of particles */
    printf("Setting initial positions of particles...\n");
	gen_mff_ic_nt(par_num, Ng,vel_field, par_pos, nt);
	
	/* Tracking initialization */
	
	 int track_num = 4;
	 int num_step = 25;
	 int* id = new int[track_num*track_num];
	 double*** track_pos = new double**[track_num*track_num];
	 for (int j=0; j<track_num*track_num; j++){
		 track_pos[j]=new double*[num_step];
		 for (int k=0; k<num_step; k++) track_pos[j][k]=new double[2];
	 }
	 int step = 0;
	 printf("Initializing IDs of tracked particles...\n");
	 get_track_par_id(par_num, id, track_num);
	
	
	/****************
	* INTEGRATION...*
	****************/
	while((b <= b_out) && (db > 0)){
		z = 1./b - 1.;
		printf("Starting computing next step with z = %.2f (b = %.2f)\n", z, b);
		
		/* Updating positions of particles... */
		printf("Updating positions of particles...\n");
		upd_pos_mff_nt(par_num, mesh_num, par_pos, vel_field, db, 1, nt);

		if (((i % 1) == 0) or (b == b_out)){
			z_suffix_num.str("");
			z_suffix_num << fixed << setprecision(2) << z;
			z_suffix = z_suffix_const + "z" + z_suffix_num.str();
			
			/* Printing positions */
			print_par_pos_cut_small(par_num, mesh_num, L, par_pos, out_dir, z_suffix);
			if (step < num_step){
			update_track_par(track_pos, par_pos, step, id, track_num);
			step++;
			print_track_par(track_pos, step, track_num, mesh_num, L, out_dir, z_suffix);
			} else printf("Maximal number of tracking steps reched!\n");
			
			/* Printing power spectrum */
			get_rho_par(mesh_num, par_num, par_pos, power_aux, nt, 1);
			print_rho_map(mesh_num, L,  power_aux, out_dir, z_suffix);
			
			pwr_spec(mesh_num, par_num, L, power_aux, p_F, nt, 1);
			gen_pow_spec_binned(mesh_num, reinterpret_cast<fftw_complex*>(power_aux), pwr_spec_binned,  2.*PI/L, 2.*PI*mesh_num/L, bin_num);			
			print_pow_spec(mesh_num, pwr_spec_binned, out_dir, z_suffix, 2.*PI/L, 2.*PI*mesh_num/L, bin_num);
			print_pow_spec_diff(pwr_spec_binned, pwr_spec_binned_0, out_dir, "_diff" + z_suffix, bin_num, b);
		}
		i++;
		if ((b_out - b) < db) db = b_out - b;
	//	else if (db > 0.1) db = (b_out - b)/100.;
		else db = 0.01;
		b += db;
	}
		
	/* Cleanup */
	dealloc_zeldovich(arrays);
}


void frozen_potential(int mesh_num, int Ng, int L, t_power power_spectrum, double* parameters, double z_in, double z_out, int nt, string out_dir){
	cout << "\n"
	"******************************\n"
	"FROZEN-POTENTIAL APPROXIMATION\n"
	"******************************\n";

	/* ALLOCATION OF MEMORY + FFTW PREPARATION */
	int par_num = mesh_num / Ng;
	int bin_num = 100;
	void** arrays = alloc_frozen_pot(mesh_num, par_num, bin_num, nt);
	if (arrays == NULL) {printf("Exiting...\n"); return;}
	double** force_field = (double**)arrays[0];
	double*** par_pos = (double***)arrays[1];
	double* power_aux = (double*)arrays[2];
	fftw_complex* pwr_spec_binned = (fftw_complex*)arrays[3];
	fftw_complex* pwr_spec_binned_0 = (fftw_complex*)arrays[4];
	fftw_plan p_F = (fftw_plan)arrays[5];
	fftw_plan p_B = (fftw_plan)arrays[6];
	
	const double b_in = 1./(z_in + 1);
	const double b_out = 1./(z_out + 1);
	double b = b_in;
	double db = b_in;
	double z;
	int i = 0;

	string z_suffix_const = "_FP_";
	string z_suffix;
	stringstream z_suffix_num;
	/* END OF ALLOCATION */
	
	/* Generating the right density distribution in k-space */
	gen_rho_dist_k(mesh_num, L, force_field[0], p_F, power_spectrum, parameters, 0., nt); // displ_vec[0] used as rho(k)
	
	/* Computing initial power spectrum in k-space */
	pwr_spec_0_nt(mesh_num, L, force_field[0], reinterpret_cast<fftw_complex*>(power_aux), nt, 1);
	gen_pow_spec_binned(mesh_num, reinterpret_cast<fftw_complex*>(power_aux), pwr_spec_binned_0,  2.*PI/L, 2.*PI*mesh_num/L, bin_num);
	
	/* Computing initial potential in k-space */
	printf("Computing initial potential in k-space...\n");
	gen_pot_k_nt(mesh_num, force_field[0], nt);
	
	/* Computing displacement in k-space */
	printf("Computing displacement in k-space...\n");	
	gen_displ_k_w_nt(mesh_num, force_field, force_field[0], nt);
	
	/* Computing displacement in q-space */
	printf("Computing displacement in q-space...\n");
	
	fftw_execute_dft_c2r_triple(p_B, reinterpret_cast<fftw_complex**>(force_field), force_field);
	normalize_FFT_BACKWARD_nt(mesh_num, force_field, nt);
	
	/* Setting initial positions and velocities of particles */
    printf("Setting initial positions and velocities of particles...\n");
	gen_init_con_nt(par_num, Ng, force_field, par_pos, nt);
	
	/* Tracking initialization */
	
	 int track_num = 4;
	 int num_step = 25;
	 int* id = new int[track_num*track_num];
	 double*** track_pos = new double**[track_num*track_num];
	 for (int j=0; j<track_num*track_num; j++){
		 track_pos[j]=new double*[num_step];
		 for (int k=0; k<num_step; k++) track_pos[j][k]=new double[2];
	 }
	 int step = 0;
	 printf("Initializing IDs of tracked particles...\n");
	 get_track_par_id(par_num, id, track_num);
	
	/* Suppresion initialization */
	typedef double t_supp[2];
	t_supp* supp = new t_supp[num_step];
	
	/* Statistical properties initialization */
	int* dens_binned = new int[500]; //5 bin per 1, [-1, 99]
	
	/****************
	* INTEGRATION...*
	****************/
	while((b <= b_out) && (db > 0)){
		z = 1./b - 1.;
		printf("Starting computing next step with z = %.2f (b = %.2f)\n", z, b);
		
		/* Updating positions of particles... */
		printf("Updating positions of particles...\n");
		upd_pos_leapfrog_nt(par_num, mesh_num, par_pos, force_field, b - db/2., db, 1, nt);
	
		
		if (((i % 2) == 0) or (b == b_out)){
			z_suffix_num.str("");
			z_suffix_num << fixed << setprecision(2) << z;
			z_suffix = z_suffix_const + "z" + z_suffix_num.str();
			
			/* Printing positions */
			print_par_pos_cut_small(par_num, mesh_num, L, par_pos[0], out_dir, z_suffix);
			
			/* Printing power spectrum */
			get_rho_par(mesh_num, par_num, par_pos[0], power_aux, nt, 1);
			gen_dens_binned(mesh_num, Ng, power_aux, dens_binned, 500);
			
			print_rho_map(mesh_num, L,  power_aux, out_dir, z_suffix);
			print_dens_bin(dens_binned, out_dir, z_suffix, 500, mesh_num);
			
			pwr_spec(mesh_num, par_num, L, power_aux, p_F, nt, 1);
			gen_pow_spec_binned(mesh_num, reinterpret_cast<fftw_complex*>(power_aux), pwr_spec_binned,  2.*PI/L, 2.*PI*mesh_num/L, bin_num);
			print_pow_spec(mesh_num, pwr_spec_binned, out_dir, z_suffix, 2.*PI/L, 2.*PI*mesh_num/L, bin_num);
			print_pow_spec_diff(pwr_spec_binned, pwr_spec_binned_0, out_dir, "_diff" + z_suffix, bin_num, b);
			
			/* Tracking and suppresion */
			if (step < num_step){
			update_track_par(track_pos, par_pos[0], step, id, track_num);
			upd_suppresion(supp, pwr_spec_binned, pwr_spec_binned_0, b, step);
			step++;
			print_track_par(track_pos, step, track_num, mesh_num, L, out_dir, z_suffix);
			
			} else printf("Maximal number of tracking steps reched!\n");
			
		}
		i++;
		if ((b_out - b) < db) db = b_out - b;
		else db = 0.05;
		b += db;
	}
	double resolution = (double)mesh_num / L;
	
	z_suffix_num.str("");
	z_suffix_num << fixed << setprecision(2) << resolution;
	z_suffix = "_res" + z_suffix_num.str();
	z_suffix_num.str("");
	z_suffix_num << fixed << setprecision(0) << Ng;
	z_suffix = z_suffix + "_R" + z_suffix_num.str();
	
	print_suppression(supp, out_dir, step, z_suffix);
	
	/* Cleanup */
	printf("Cleaning...\n");
	dealloc_frozen_pot(arrays);
	delete[] id;	
	for (int j=0; j<track_num*track_num; j++){
		for (int k=0; k<num_step; k++) delete[] track_pos[j][k];
		delete[] track_pos[j];
	}
	delete[] track_pos;
	delete[] supp;
	delete[] dens_binned;
}

void adhesiom_approximation(int mesh_num, int Ng, int L, t_power power_spectrum, double* parameters, double nu, double z_in, double z_out, int nt, string out_dir){
	cout << "\n"
	"**********************\n"
	"ADHESION APPROXIMATION\n"
	"**********************\n";
	
	/* ALLOCATION OF MEMORY + FFTW PREPARATION */
	int par_num = mesh_num / Ng;
	int bin_num = 100;
	void** arrays = alloc_adhesion(mesh_num, par_num, bin_num, nt);
	if (arrays == NULL) {printf("Exiting...\n"); return;}
	double** vel_field = (double**)arrays[0];
	double* potential = vel_field[0];
	double** par_pos = (double**)arrays[1];
	double* power_aux = (double*)arrays[2];
	double* expotential_0 = (double*)arrays[3];
	fftw_complex* pwr_spec_binned = (fftw_complex*)arrays[4];
	fftw_complex* pwr_spec_binned_0 = (fftw_complex*)arrays[5];
	fftw_plan p_F = (fftw_plan)arrays[6];
	fftw_plan p_B = (fftw_plan)arrays[7];
	
	const double b_in = 1./(z_in + 1);
	const double b_out = 1./(z_out + 1);
	double b = b_in;
	double db = b_in;
	double z;
	int i = 0;
	bool conv_mod = false; // false == direct sum, true == FFT
		
	string z_suffix_const = "_AA_";
	string z_suffix;
	stringstream z_suffix_num;
	
	/* END OF ALLOCATION */
	
	/* Generating the right density distribution in k-space */
	gen_rho_dist_k(mesh_num, L, potential, p_F, power_spectrum, parameters, 0., nt);
	
	pwr_spec_0_nt(mesh_num, L, potential, reinterpret_cast<fftw_complex*>(power_aux), nt, 1);
	gen_pow_spec_binned(mesh_num, reinterpret_cast<fftw_complex*>(power_aux), pwr_spec_binned_0,  2.*PI/L, 2.*PI*mesh_num/L, bin_num);
	
	/* Computing initial potential in k-space */
	printf("Computing initial potential in k-space...\n");
	gen_pot_k_nt(mesh_num, potential, nt);
	
	/* Setting initial positions of particles */
    printf("Setting initial positions of particles...\n");
	gen_no_displ_pos_nt(par_num, Ng, par_pos, nt);
	
	/* Computing initial expotential */
	fftw_execute_dft_c2r(p_B, reinterpret_cast<fftw_complex*>(potential), potential);
	normalize_FFT_BACKWARD_nt(mesh_num, potential, nt);	
	gen_init_expot(mesh_num, potential, expotential_0, nu, conv_mod, p_F, nt);
	
	if (check_field(expotential_0, mesh_num)) return;
	// expotential_0 now stores initial values which we will be using later
	
	/* Tracking initialization */
	
	 int track_num = 4;
	 int num_step = 25;
	 int* id = new int[track_num*track_num];
	 double*** track_pos = new double**[track_num*track_num];
	 for (int j=0; j<track_num*track_num; j++){
		 track_pos[j]=new double*[num_step];
		 for (int k=0; k<num_step; k++) track_pos[j][k]=new double[2];
	 }
	 int step = 0;
	 printf("Initializing IDs of tracked particles...\n");
	 get_track_par_id(par_num, id, track_num);
	
	/* Suppresion initialization */
	typedef double t_supp[2];
	t_supp* supp = new t_supp[num_step];
	
	/* Statistical properties initialization */
	int* dens_binned = new int[500]; //5 bin per 1, [-1, 99]
	 
	/****************
	* INTEGRATION...*
	****************/
	while((b <= b_out) && (db > 0)){
		z = 1./b - 1.;
		printf("Starting computing next step with z = %.2f (b = %.2f)\n",z, b);
		/* Computing convolution */
		gen_expot(mesh_num, potential, expotential_0, nu, b - db/2., conv_mod, p_B, nt);
		if (check_field(potential, mesh_num)) return;
		
		/* Computing potential */
		printf("Computing potential...\n");	
		exp2pot_nt(mesh_num, potential, nu, nt);
		if (check_field(potential, mesh_num)) return;
		
		/* Computing velocity field via finite difference */
	//	printf("Computing velocity field via finite difference...\n");
	//	fin_diff(mesh_num, vel_field, potential);
		
		/* Computing velocity field via FFT (+CIC corrections) */
		printf("Computing velocity field via FFT...\n");
		
		fftw_execute_dft_r2c(p_F, potential, reinterpret_cast<fftw_complex*>(potential));
		normalize_FFT_FORWARD_nt(mesh_num, reinterpret_cast<fftw_complex*>(potential), nt);
		
		gen_displ_k_w_nt(mesh_num, vel_field, potential, nt);
		
		fftw_execute_dft_c2r_triple(p_B, reinterpret_cast<fftw_complex**>(vel_field), vel_field);
		normalize_FFT_BACKWARD_nt(mesh_num, vel_field, nt);

		/* Updating positions of particles... */
		printf("Updating positions of particles...\n");
		upd_pos_nt(par_num, mesh_num, par_pos, vel_field, db, 1, nt);
	//	upd_pos_vel0_nt(par_num, mesh_num, par_pos, vel_field, db, 1, nt); // using initial particle position
		
		if (((i % 2) == 0) or (b == b_out)){
			z_suffix_num.str("");
			z_suffix_num << fixed << setprecision(2) << z;
			z_suffix = z_suffix_const + "z" + z_suffix_num.str();
			
			/* Printing positions */
			print_par_pos_cut_small(par_num, mesh_num, L, par_pos, out_dir, z_suffix);
			
			/* Printing power spectrum */
			get_rho_par(mesh_num, par_num, par_pos, power_aux, nt, 1);
			gen_dens_binned(mesh_num, Ng, power_aux, dens_binned, 500);
			
			print_rho_map(mesh_num, L,  power_aux, out_dir, z_suffix);
			print_dens_bin(dens_binned, out_dir, z_suffix, 500, mesh_num);
			
			pwr_spec(mesh_num, par_num, L, power_aux, p_F, nt, 1);
			gen_pow_spec_binned(mesh_num, reinterpret_cast<fftw_complex*>(power_aux), pwr_spec_binned,  2.*PI/L, 2.*PI*mesh_num/L, bin_num);
			print_pow_spec(mesh_num, pwr_spec_binned, out_dir, z_suffix, 2.*PI/L, 2.*PI*mesh_num/L, bin_num);
			print_pow_spec_diff(pwr_spec_binned, pwr_spec_binned_0, out_dir, "_diff" + z_suffix, bin_num, b);
			
			/* Tracking and suppresion */
			if (step < num_step){
			update_track_par(track_pos, par_pos, step, id, track_num);
			upd_suppresion(supp, pwr_spec_binned, pwr_spec_binned_0, b, step);
			step++;
			print_track_par(track_pos, step, track_num, mesh_num, L, out_dir, z_suffix);
			
			} else printf("Maximal number of tracking steps reched!\n");
		}
		
		i++;
		if ((b_out - b) < db) db = b_out - b;
	//	else if (db > 0.1) db = (b_out - b)/100.;
		else db = 0.05;
		b += db;
	}
	
	double resolution = (double)mesh_num / L;
	
	z_suffix_num.str("");
	z_suffix_num << fixed << setprecision(2) << resolution;
	z_suffix = "_res" + z_suffix_num.str();
	z_suffix_num.str("");
	z_suffix_num << fixed << setprecision(0) << Ng;
	z_suffix = z_suffix + "_R" + z_suffix_num.str();
	
	print_suppression(supp, out_dir, step, z_suffix);
	
	/* Cleanup */
	printf("Cleaning...\n");
	dealloc_adhesion(arrays);
	delete[] id;	
	for (int j=0; j<track_num*track_num; j++){
		for (int k=0; k<num_step; k++) delete[] track_pos[j][k];
		delete[] track_pos[j];
	}
	delete[] track_pos;
	delete[] supp;
	delete[] dens_binned;
}

void frozen_flow(int mesh_num, int Ng, int L, t_power power_spectrum, double* parameters, double z_in, double z_out, int nt, string out_dir){
	cout << "\n"
	"*************************\n"
	"FROZEN-FLOW APPROXIMATION\n"
	"*************************\n";

	/* ALLOCATION OF MEMORY + FFTW PREPARATION */
	int par_num = mesh_num / Ng;
	int bin_num = 100;
	void** arrays = alloc_zeldovich(mesh_num, par_num, bin_num, nt);
	if (arrays == NULL) {printf("Exiting...\n"); return;}
	double** vel_field = (double**)arrays[0];
	double** par_pos = (double**)arrays[1];
	double* power_aux = (double*)arrays[2];
	fftw_complex* pwr_spec_binned = (fftw_complex*)arrays[3];
	fftw_complex* pwr_spec_binned_0 = (fftw_complex*)arrays[4];
	fftw_plan p_F = (fftw_plan)arrays[5];
	fftw_plan p_B = (fftw_plan)arrays[6];
	
	const double b_in = 1./(z_in + 1);
	const double b_out = 1./(z_out + 1);
	double b = b_in;
	double db = b_in;
	double z;
	int i = 0;

	string z_suffix_const = "_FF_";
	string z_suffix;
	stringstream z_suffix_num;
	/* END OF ALLOCATION */
	
	/* Generating the right density distribution in k-space */
	gen_rho_dist_k(mesh_num, L, vel_field[0], p_F, power_spectrum, parameters, 0., nt); // displ_vec[0] used as rho(k)
	
	/* Computing initial power spectrum in k-space */
	pwr_spec_0_nt(mesh_num, L, vel_field[0], reinterpret_cast<fftw_complex*>(power_aux), nt, 1);
	gen_pow_spec_binned(mesh_num, reinterpret_cast<fftw_complex*>(power_aux), pwr_spec_binned_0,  2.*PI/L, 2.*PI*mesh_num/L, bin_num);
	
	/* Computing initial potential in k-space */
	printf("Computing initial potential in k-space...\n");
	gen_pot_k_nt(mesh_num, vel_field[0], nt);
	
	/* Computing displacement in k-space */
	printf("Computing displacement in k-space...\n");	
	gen_displ_k_w_nt(mesh_num, vel_field, vel_field[0], nt);
	
	/* Computing displacement in q-space */
	printf("Computing displacement in q-space...\n");
	
	fftw_execute_dft_c2r_triple(p_B, reinterpret_cast<fftw_complex**>(vel_field), vel_field);
	normalize_FFT_BACKWARD_nt(mesh_num, vel_field, nt);
	
	/* Setting initial positions of particles */
    printf("Setting initial positions of particles...\n");
	gen_no_displ_pos_nt(par_num, Ng, par_pos, nt);
	
	/* Tracking initialization */
	
	 int track_num = 4;
	 int num_step = 25;
	 int* id = new int[track_num*track_num];
	 double*** track_pos = new double**[track_num*track_num];
	 for (int j=0; j<track_num*track_num; j++){
		 track_pos[j]=new double*[num_step];
		 for (int k=0; k<num_step; k++) track_pos[j][k]=new double[2];
	 }
	 int step = 0;
	 printf("Initializing IDs of tracked particles...\n");
	 get_track_par_id(par_num, id, track_num);
	
	/* Suppresion initialization */
	typedef double t_supp[2];
	t_supp* supp = new t_supp[num_step];
	
	/* Statistical properties initialization */
	int* dens_binned = new int[500]; //5 bin per 1, [-1, 99]
	
	
	/****************
	* INTEGRATION...*
	****************/
	while((b <= b_out) && (db > 0)){
		z = 1./b - 1.;
		printf("Starting computing next step with z = %.2f (b = %.2f)\n", z, b);
		
		/* Updating positions of particles... */
		printf("Updating positions of particles...\n");
		upd_pos_nt(par_num, mesh_num, par_pos, vel_field, db, 1, nt);
	//	upd_pos_vel0_nt(par_num, mesh_num, par_pos, vel_field, db, 1, nt); // using initial particle positions
			
		if (((i % 2) == 0) or (b == b_out)){
			z_suffix_num.str("");
			z_suffix_num << fixed << setprecision(2) << z;
			z_suffix = z_suffix_const + "z" + z_suffix_num.str();
			
			/* Printing positions */
			print_par_pos_cut_small(par_num, mesh_num, L, par_pos, out_dir, z_suffix);
			
			/* Printing power spectrum */
			get_rho_par(mesh_num, par_num, par_pos, power_aux, nt, 1);
			gen_dens_binned(mesh_num, Ng, power_aux, dens_binned, 500);
			
			print_rho_map(mesh_num, L,  power_aux, out_dir, z_suffix);
			print_dens_bin(dens_binned, out_dir, z_suffix, 500, mesh_num);
			
			pwr_spec(mesh_num, par_num, L, power_aux, p_F, nt, 1);
			gen_pow_spec_binned(mesh_num, reinterpret_cast<fftw_complex*>(power_aux), pwr_spec_binned,  2.*PI/L, 2.*PI*mesh_num/L, bin_num);
			print_pow_spec(mesh_num, pwr_spec_binned, out_dir, z_suffix, 2.*PI/L, 2.*PI*mesh_num/L, bin_num);
			print_pow_spec_diff(pwr_spec_binned, pwr_spec_binned_0, out_dir, "_diff" + z_suffix, bin_num, b);
			
			/* Tracking and suppresion */
			if (step < num_step){
			update_track_par(track_pos, par_pos, step, id, track_num);
			upd_suppresion(supp, pwr_spec_binned, pwr_spec_binned_0, b, step);
			step++;
			print_track_par(track_pos, step, track_num, mesh_num, L, out_dir, z_suffix);
			
			} else printf("Maximal number of tracking steps reched!\n");
		}
		i++;
		if ((b_out - b) < db) db = b_out - b;
	//	else if (db > 0.1) db = (b_out - b)/100.;
		else db = 0.05;
		b += db;
	}
	double resolution = (double)mesh_num / L;
	
	z_suffix_num.str("");
	z_suffix_num << fixed << setprecision(2) << resolution;
	z_suffix = "_res" + z_suffix_num.str();
	z_suffix_num.str("");
	z_suffix_num << fixed << setprecision(0) << Ng;
	z_suffix = z_suffix + "_R" + z_suffix_num.str();
	
	print_suppression(supp, out_dir, step, z_suffix);
	
	/* Cleanup */
	printf("Cleaning...\n");
	dealloc_zeldovich(arrays);
	delete[] id;	
	for (int j=0; j<track_num*track_num; j++){
		for (int k=0; k<num_step; k++) delete[] track_pos[j][k];
		delete[] track_pos[j];
	}
	delete[] track_pos;
	delete[] supp;
	delete[] dens_binned;
}

void truncated_zeldovich(int mesh_num, int Ng, int L, t_power power_spectrum, double* parameters, double z_in, double z_out, int nt, string out_dir, double k2_G){
	if (k2_G == 0) cout << "\n"
	"************************\n"
	"ZEL`DOVICH APPROXIMATION\n"
	"************************\n";
	else cout << "\n"
	"**********************************\n"
	"TRUNCATED ZEL`DOVICH APPROXIMATION\n"
	"**********************************\n";

	/* ALLOCATION OF MEMORY + FFTW PREPARATION */
	int par_num = mesh_num / Ng;
	int bin_num = 100;
	void** arrays = alloc_zeldovich(mesh_num, par_num, bin_num, nt);
	if (arrays == NULL) {printf("Exiting...\n"); return;}
	double** vel_field = (double**)arrays[0];
	double** par_pos = (double**)arrays[1];
	double* power_aux = (double*)arrays[2];
	fftw_complex* pwr_spec_binned = (fftw_complex*)arrays[3];
	fftw_complex* pwr_spec_binned_0 = (fftw_complex*)arrays[4];
	fftw_plan p_F = (fftw_plan)arrays[5];
	fftw_plan p_B = (fftw_plan)arrays[6];
	
	const double b_in = 1./(z_in + 1);
	const double b_out = 1./(z_out + 1);
	double b = b_in;
	double db = b_in;
	double z;
	
	string z_suffix_const = "_ZA_";
	string z_suffix;
	stringstream z_suffix_num;
	/* END OF ALLOCATION */
	
	/* Generating the right density distribution in k-space */
	gen_rho_dist_k(mesh_num, L, vel_field[0], p_F, power_spectrum, parameters, k2_G, nt); // displ_vec[0] used as rho(k)
	
	pwr_spec_0_nt(mesh_num, L, vel_field[0], reinterpret_cast<fftw_complex*>(power_aux), nt, 1);
	gen_pow_spec_binned(mesh_num, reinterpret_cast<fftw_complex*>(power_aux), pwr_spec_binned_0,  2.*PI/L, 2.*PI*mesh_num/L, bin_num);
	
	/* Computing initial potential in k-space */
	printf("Computing initial potential in k-space...\n");
	gen_pot_k_nt(mesh_num, vel_field[0], nt);
	
	/* Computing displacement in k-space */
	printf("Computing displacement in k-space...\n");	
	gen_displ_k_w_nt(mesh_num, vel_field, vel_field[0], nt); // CIC optimalization
//	gen_displ_k_nt(mesh_num, vel_field, vel_field[0], nt); // no corrections
	
	/* Computing displacement in q-space */
	printf("Computing displacement in q-space...\n");
	
	fftw_execute_dft_c2r_triple(p_B, reinterpret_cast<fftw_complex**>(vel_field), vel_field);
	normalize_FFT_BACKWARD_nt(mesh_num, vel_field, nt);
		
	/* Tracking initialization */
	 int track_num = 4;
	 int num_step = 25;
	 int* id = new int[track_num*track_num];
	 double*** track_pos = new double**[track_num*track_num];
	 for (int j=0; j<track_num*track_num; j++){
		 track_pos[j]=new double*[num_step];
		 for (int k=0; k<num_step; k++) track_pos[j][k]=new double[2];
	 }
	 int step = 0;
	 printf("Initializing IDs of tracked particles...\n");
	 get_track_par_id(par_num, id, track_num);
	 
	/****************
	* INTEGRATION...*
	****************/
		
	while((b <= b_out) && (db > 0)){
		z = 1./b - 1.;
		printf("Starting computing next step with z = %.2f (b = %.2f)\n", z, b);
		
		/* Computing displaced positions of particles */
		printf("Computing displaced positions of particles...\n");
		gen_displ_pos_ng_nt(par_num, Ng, vel_field, par_pos, b, nt);
		
		z_suffix_num.str("");
		z_suffix_num << fixed << setprecision(2) << z;
		z_suffix = z_suffix_const + "z" + z_suffix_num.str();
		
		/* Printing positions */
		print_par_pos_cut_small(par_num, mesh_num, L, par_pos, out_dir, z_suffix);
	//	print_par_pos_cut(par_num, mesh_num, L, par_pos, out_dir, z_suffix);	
		if (step < num_step){
		update_track_par(track_pos, par_pos, step, id, track_num);
		step++;
		print_track_par(track_pos, step, track_num, mesh_num, L, out_dir, z_suffix);
		} else printf("Maximal number of tracking steps reched!\n");
		
		/* Printing power spectrum */
		get_rho_par(mesh_num, par_num, par_pos, power_aux, nt, 1);
		print_rho_map(mesh_num, L,  power_aux, out_dir, z_suffix);
		
		pwr_spec(mesh_num, par_num, L, power_aux, p_F, nt, 1);
		gen_pow_spec_binned(mesh_num, reinterpret_cast<fftw_complex*>(power_aux), pwr_spec_binned,  2.*PI/L, 2.*PI*mesh_num/L, bin_num);			
		print_pow_spec(mesh_num, pwr_spec_binned, out_dir, z_suffix, 2.*PI/L, 2.*PI*mesh_num/L, bin_num);
		print_pow_spec_diff(pwr_spec_binned, pwr_spec_binned_0, out_dir, "_diff" + z_suffix, bin_num, b);
		if ((b_out - b) < db) db = b_out - b;
		//	else if (db > 0.1) db = (b_out - b)/100.;
		else db = 0.01;
		b += db;
	}
	/* Cleanup */
	dealloc_zeldovich(arrays);
}


void zel_app(int mesh_num, int Ng, int L, t_power power_spectrum, double* parameters, double z_in, double z_out, int nt, string out_dir){
	truncated_zeldovich (mesh_num, Ng, L, power_spectrum, parameters, z_in, z_out, nt, out_dir, 0.);
}
