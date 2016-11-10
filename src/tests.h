	/*
	
	printf("Printing some values:\n");
	for (int i = 0; i<10; i++) printf("Potential[%i] = %.15f\n", i, potential[i]);
	fftw_execute(p_F); // we do not care about normalization because of logarithmic derivative
	normalize_FFT_FORWARD(par_num, reinterpret_cast<fftw_complex*>(potential));
	fftw_execute(p_B);
	normalize_FFT_BACKWORD(par_num, potential);
	printf("Printing some values:\n");
	for (int i = 0; i<10; i++) printf("Potential[%i] = %.15f\n", i, potential[i]);

	
	
	printf("Printing some values of potential in q-space:\n");
	for (int i = 1; i<par_num*par_num*(par_num+2); i*=par_num/16) printf("Potential[%i] = %f\n", i-1, potential[i-1]);
	
	
	
	
	// GAUSSIAN RANDOM FIELD //
		
//		pwr_spec(mesh_num, par_num, box_size, displ_vec, power, nt, 1, true, k2_G);
//		print_pow_spec(mesh_num, reinterpret_cast<fftw_complex*>(power), out_dir, "_Pk", 2.*PI/box_size, 2.*PI*mesh_num/box_size, 100);



void pwr_spec(int mesh_num, int par_num, int L, double** par_pos, double* delta, fftw_plan p_F, int nt, const int order, bool GAUSS, t_power power_spectrum, double* parameters, double k2_G){
	try{
		fftw_complex* delta_k = reinterpret_cast<fftw_complex*>(delta);
		
		// Computing the density field
		if (GAUSS){
			cout << "\n"
			"*********************\n"
			"GAUSSIAN RANDOM FIELD\n"
			"*********************\n";			
			gen_rho_dist_k(mesh_num, L, delta, p_F, power_spectrum, parameters, k2_G, nt); // displ_vec[0] used as rho(k)
		} else{
			get_rho_par(mesh_num, par_num, par_pos, delta, nt, order);
			fftw_execute_dft_r2c(p_F, delta, delta_k);
			normalize_FFT_FORWARD_nt(mesh_num, delta_k, nt);
		}		
		
		// Computing the power spectrum P(k)
		printf("Computing the power spectrum P(k)...\n");
		comp_power_spec(mesh_num, delta_k, L, order, nt);
	}
	catch(...){
		printf("ERROR!\n");
	}
}


	
	// RANDOM POSITIONS //
		#ifndef ZEL
		gen_rnd_pos(mesh_num, par_num, displ_vec);
		#endif		

*/

