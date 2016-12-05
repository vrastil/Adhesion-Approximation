#include "stdafx.h"
#include <fftw3.h>
#include "output.h"
#include "grid_fce.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>

using namespace std;

typedef double(*t_power)(double, double*);
const double PI = acos(-1.);

/* POWER SPECTRUM FUNCTIONS */

const double Omega_0 = 1;
const double h = 0.67;
extern double s8;

double transfer_function_2(double k, double *parameters){
	if (k == 0) return 1.;
	double q = k / (Omega_0*h);
	double T_k =	log(1+2.34*q)/(2.34*q)*
					pow(1 + 3.89*q + pow(16.1*q, 2.) + pow(5.4*q, 3.) + pow(6.71*q, 4.)
					, -1./4.);
	return pow(T_k, 2.);
}

double power_spectrum_T(double k, double* parameters){
	/* scale-free power spectrum with trunsfer function */
	if (k == 0) return 0;
	double A = parameters[0];
	double ns = parameters[1];
	return A*pow(k, ns)*transfer_function_2(k, NULL);
}

double power_spectrum(double k, double* parameters){
	/* scale-free power spectrum */
	if (k == 0) return 0;
	double A = parameters[0];
	double ns = parameters[1];
	return A*pow(k, ns);
}

double flat_power_spectrum(double k, double* parameters){
	return parameters[2];
}

double single_power_spectrum(double k, double* parameters){
	double k_min = 0.01;
	double k_max = 0.04;
	
	if ((k > k_min) and (k < k_max)) return power_spectrum_T(k, parameters);
	return 0.;
}

double power_spectrum_integ(double k, void* parameters){
	double const R = 8; // R = 8 Mpc/h
//	double const L = ((double*)parameters)[2];
//	double k_ = 2*PI*k/L;
	return	k*k/(2.*PI*PI)* // spherical factor
			power_spectrum_T(k, (double*)parameters)* // P(k)
			pow(3.*gsl_sf_bessel_j1(k*R)/(k*R), 2.); // window function
}
void power_spectrum_norm(double* parameters){
	/* Normalize the power spectrum */
	printf("Computing normalization of the given power spectrum...\n");
	parameters[0] = 1.;
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	double result, error;
	
	gsl_function F;
	F.function = &power_spectrum_integ;
	F.params = parameters;
	
	gsl_integration_qagiu (&F, 0, 0, 1e-7, 1000, w, &result, &error); 
	
	gsl_integration_workspace_free (w);
	parameters[0] = s8*s8/result;
}

void gen_sin_rho_k(int mesh_num, int L, double* rho_x, fftw_plan p_F, double k2_G, int nt){
	for(int i=0; i < pow(mesh_num, 2); i++){
		for(int j=0; j<mesh_num; j++){			
			rho_x[i*(mesh_num+2)+j]=sin(k2_G*j*(double)L/mesh_num);
		}
	}
	fftw_execute_dft_r2c(p_F, rho_x, reinterpret_cast<fftw_complex*>(rho_x));
	normalize_FFT_FORWARD_nt(mesh_num, reinterpret_cast<fftw_complex*>(rho_x), nt);
}

void comp_power_spec_one(int mesh_num, fftw_complex* delta_k, int i_min, int i_max, int L, int order){
	int k_vec[3];
	double w_k;
	for(int i=i_min; i < i_max;i++){
		w_k = 1.;
		get_coord_fftw_complex(mesh_num, i, k_vec);
		
		for (int j = 0; j < 3; j++){
				if(k_vec[j] > mesh_num/2.) k_vec[j] -= mesh_num; // symmetry
				if (k_vec[j] != 0) w_k *= pow(sin(PI*k_vec[j]/mesh_num)/(PI*k_vec[j]/mesh_num), order + 1); // k = 0 => w(k) = 1
		}
		
		// CIC + Power Spectrum units [L/N]^3
	//	delta_k[i][1] = (delta_k[i][0]*delta_k[i][0] + delta_k[i][1]*delta_k[i][1])/(w_k*w_k)*pow((double)L/(double)mesh_num, 3);
		// CIC
		delta_k[i][1] = (delta_k[i][0]*delta_k[i][0] + delta_k[i][1]*delta_k[i][1])/(w_k*w_k);
		
		delta_k[i][0] = 0;
		for (int j = 0; j < 3; j++) delta_k[i][0] += pow(k_vec[j], 2);
		delta_k[i][0] = 2.*PI/L*sqrt(delta_k[i][0]); // physical k
	}	
}

void comp_power_spec(int mesh_num, fftw_complex* delta_k, int L, int order, int nt){
	int i_min = 0;
	int i_max = pow(mesh_num, 2)*(mesh_num/2 + 1);	
	int i_inc = i_max / nt;
	i_max = i_max % nt + i_inc;
	thread* th_dens = new thread[nt];
	for (int i = 0; i < nt; i++){
		th_dens[i] = thread(comp_power_spec_one, mesh_num, delta_k, i_min, i_max, L, order);
		i_min = i_max;
		i_max += i_inc;
	}
	for (int i = 0; i < nt; i++) th_dens[i].join();
	delete[] th_dens;
}


void pwr_spec_0_one(int mesh_num, int L, double* delta_k, fftw_complex* power, int i_min, int i_max, const int order){
	/* Computing the power spectrum P(k) */		
	int k_vec[3];
	double w_k;
	for(int i=i_min; i < i_max;i++){
		w_k = 1.;
		get_coord_fftw_complex(mesh_num, i, k_vec);
		
		for (int j = 0; j < 3; j++){
				if(k_vec[j] > mesh_num/2.) k_vec[j] -= mesh_num; // symmetry
				if (k_vec[j] != 0) w_k *= pow(sin(PI*k_vec[j]/mesh_num)/(PI*k_vec[j]/mesh_num), order + 1); // k = 0 => w(k) = 1
		}
		
		// CIC + Power Spectrum units [L/N]^3
	//	power[i][1] = (delta_k[2*i]*delta_k[2*i] + delta_k[2*i+1]*delta_k[2*i+1])/(w_k*w_k)*pow((double)L/(double)mesh_num, 3);
		// CIC
		power[i][1] = (delta_k[2*i]*delta_k[2*i] + delta_k[2*i+1]*delta_k[2*i+1])/(w_k*w_k);
		
		power[i][0] = 0;
		for (int j = 0; j < 3; j++) power[i][0] += pow(k_vec[j], 2);
		power[i][0] = 2.*PI/L*sqrt(power[i][0]); // physical k
	}
}

void pwr_spec_0_nt(int mesh_num, int L, double* delta_k, fftw_complex* power, int nt, const int order){
	/* Variant with already computed \delta(k) */
	/* Preserve values in delta_k */
	printf("Computing the initial power spectrum P(k)...\n");
	int i_min = 0;
	int i_max = pow(mesh_num, 2)*(mesh_num/2 + 1);	
	int i_inc = i_max / nt;
	i_max = i_max % nt + i_inc;
	thread* th_dens = new thread[nt];
	for (int i = 0; i < nt; i++){
		th_dens[i] = thread(pwr_spec_0_one, mesh_num, L, delta_k, power, i_min, i_max, order);
		i_min = i_max;
		i_max += i_inc;
	}
	for (int i = 0; i < nt; i++) th_dens[i].join();
	delete[] th_dens;
}

void pwr_spec(int mesh_num, int par_num, int L, double** par_pos, double* delta, fftw_plan p_F, int nt, const int order){
	/* AFTER routine:
	power[i][0] = k [h/Mpc]
	power[i][1] = P(k) [(Mpc/h)^3]
	i = 0..par_num*par_num*(par_num/2+1)
	*/
	try{
		fftw_complex* delta_k = reinterpret_cast<fftw_complex*>(delta);
		
		/* Computing the density field */
		get_rho_par(mesh_num, par_num, par_pos, delta, nt, order);
		fftw_execute_dft_r2c(p_F, delta, delta_k);
		normalize_FFT_FORWARD_nt(mesh_num, delta_k, nt);
		
		/* Computing the power spectrum P(k) */		
		printf("Computing the power spectrum P(k)...\n");
		comp_power_spec(mesh_num, delta_k, L, order, nt);
	}
	catch(...){
		printf("ERROR! Routine <void pwr_spec (int, int, int, double**, double*, fftw_plan, int, const int)> failed.\n");
	}
}

void pwr_spec(int mesh_num, int par_num, int L, double* delta, fftw_plan p_F, int nt, const int order){
	// overloaded version with already computed density field (in real-space)
	try{
		fftw_complex* delta_k = reinterpret_cast<fftw_complex*>(delta);
		
		/* Computing the density field */
		fftw_execute_dft_r2c(p_F, delta, delta_k);
		normalize_FFT_FORWARD_nt(mesh_num, delta_k, nt);
		
		/* Computing the power spectrum P(k) */		
		printf("Computing the power spectrum P(k)...\n");
		comp_power_spec(mesh_num, delta_k, L, order, nt);
	}
	catch(...){
		printf("ERROR! Routine <void pwr_spec(int, int, int, double*, fftw_plan, int, const int)> failed.\n");
	}
}

void get_rho_par(int mesh_num, int par_num, double** par_pos, double* delta, int nt, const int order){
	try{
		double m = pow((double)mesh_num/par_num, 3.);
		printf("Computing the density field...\n");
		init_dens_field(mesh_num, delta, nt);
		double x[3];
		for (int i = 0; i < pow(par_num, 2); i++){
			for (int j=0; j < par_num; j++){
				for (int k = 0; k < 3; k++) x[k] = par_pos[k][i*(par_num + 2)+j];
				assign_fc(delta, x, m, mesh_num, order, true);
			}
		}
	}
	catch(...){
		printf("ERROR!\n");
	}
}

/* FFT FUNCTIONS */

void fftw_execute_dft_c2r_triple(fftw_plan p_B, fftw_complex** in, double** out){
	for (int i=0; i<3; i++) fftw_execute_dft_c2r(p_B, in[i], out[i]);
}

void normalize_FFT_FORWARD_one(int mesh_num, fftw_complex* out, int i_min, int i_max){
	for(int i=i_min; i < i_max;i++){
		out[i][0] /= pow(mesh_num, 1.5); //symetric normalization
		out[i][1] /= pow(mesh_num, 1.5);
	}
}

void normalize_FFT_BACKWARD_one(int mesh_num, double* out, int i_min, int i_max){
	for(int i=i_min; i< i_max; i++){
		for (int j=0; j < mesh_num; j++){
			out[i*(mesh_num + 2)+j] /= pow(mesh_num, 1.5); //symetric normalization
		}
	}
}

void normalize_FFT_FORWARD_nt(int mesh_num, fftw_complex* out, int nt){
	int i_min = 0;
	int i_max = pow(mesh_num, 2)*(mesh_num/2 + 1);
	int i_inc = i_max / nt;
	i_max = i_max % nt + i_inc;
	thread* th_dens = new thread[nt];
	for (int i = 0; i < nt; i++){
		th_dens[i] = thread(normalize_FFT_FORWARD_one, mesh_num, out, i_min, i_max);
		i_min = i_max;
		i_max += i_inc;
	}
	for (int i = 0; i < nt; i++) th_dens[i].join();
	delete[] th_dens;	
}

void normalize_FFT_BACKWARD_nt(int mesh_num, double* out, int nt){
	int i_min = 0;
	int i_max = pow(mesh_num, 2.);
	int i_inc = i_max / nt;
	i_max = i_max % nt + i_inc;
	thread* th_dens = new thread[nt];
	for (int i = 0; i < nt; i++){
		th_dens[i] = thread(normalize_FFT_BACKWARD_one, mesh_num, out, i_min, i_max);
		i_min = i_max;
		i_max += i_inc;
	}
	for (int i = 0; i < nt; i++) th_dens[i].join();
	delete[] th_dens;
}

void normalize_FFT_FORWARD_nt(int mesh_num, fftw_complex** out, int nt){
	for(int i=0; i<3; i++) normalize_FFT_FORWARD_nt(mesh_num, out[i], nt);	
}

void normalize_FFT_BACKWARD_nt(int mesh_num, double** out, int nt){
	for(int i=0; i<3; i++) normalize_FFT_BACKWARD_nt(mesh_num, out[i], nt);	
}

/* GENERAL FUNCTIONS */

void copy_one(int mesh_num, double* copy, double* paste, int i_min, int i_max){
	for(int i=i_min; i< i_max; i++){
		for (int j=0; j < mesh_num; j++){
			paste[i*(mesh_num + 2)+j] = copy[i*(mesh_num + 2)+j];
		}
	}	
}

void copy_nt(int mesh_num, double* copy, double* paste, int nt){
	int i_min = 0;
	int i_max = pow(mesh_num, 2.);
	int i_inc = i_max / nt;
	i_max = i_max % nt + i_inc;
	thread* th_dens = new thread[nt];
	for (int i = 0; i < nt; i++){
		th_dens[i] = thread(copy_one, mesh_num, copy, paste, i_min, i_max);
		i_min = i_max;
		i_max += i_inc;
	}
	for (int i = 0; i < nt; i++) th_dens[i].join();
	delete[] th_dens;
}

int check_field(double* field, int max_i){
	// check field of length max_i*max_i*(max_i+2); ignore last 2 element in last dim
	double check;
	for (int i = 0; i < max_i; i++){
		check = field[i * max_i * (max_i + 2) + i * (max_i + 2) + i]; // diagonal
		if (isfinite(check) == 0){
			printf("Error while performing random check for NAN or INF! field posititon = %i, value = %f\n", i * max_i * (max_i + 2) + i * (max_i + 2) + i, check);
			return 1;
		}
	}
	return 0;
}

int isPowerOfTwo (unsigned int x){
  return ((x != 0) && ((x & (~x + 1)) == x));
}

const char *humanSize(uint64_t bytes){
	char const *suffix[] = {"B", "KB", "MB", "GB", "TB"};
	char length = sizeof(suffix) / sizeof(suffix[0]);

	int i = 0;
	double dblBytes = bytes;

	if (bytes > 1024) {
		for (i = 0; (bytes / 1024) > 0 && i<length-1; i++, bytes /= 1024)
			dblBytes = bytes / 1024.0;
	}

	static char output[200];
	sprintf(output, "%.02lf %s", dblBytes, suffix[i]);
	return output;
}

void get_coord_fftw_complex(int mesh_num, int pos, int vec[3]){
	// in vec[] are returned 3-d coordinates coresponding
	// to the position in 1-d long field of dimension [mesh_num*mesh_num*(mesh_num/2 + 1)]
	vec[0] = pos / ((mesh_num/2 + 1)*mesh_num);
	vec[1] = (pos / (mesh_num/2 + 1)) % mesh_num;
	vec[2] = pos % (mesh_num/2 + 1);
}

void get_coord_double(int mesh_num, int pos, int vec[3]){
	// in vec[] are returned 3-d coordinates coresponding
	// to the position in 1-d long field of dimension [mesh_num*mesh_num*(mesh_num + 2)]
	vec[0] = pos / ((mesh_num + 2)*mesh_num);
	vec[1] = (pos / (mesh_num + 2)) % mesh_num;
	vec[2] = pos % (mesh_num + 2);
}

int get_per(int vec, int mesh_num){
	vec = vec % mesh_num; // get into range [-N, N]
	vec = (vec + mesh_num) % mesh_num; // range [0, N]
	return vec;
}

double get_per(double vec, int mesh_num){
	vec = fmod(vec, mesh_num); // get into range [-N, N]
	vec = fmod(vec + mesh_num, mesh_num); // range [0, N]
	return vec;
}

double get_k_sq(int mesh_num, int pos){
	// returns computation k
	// includes symetry around N/2
	int k_vec[3];
	double tmp = 0;
	get_coord_fftw_complex(mesh_num, pos, k_vec);
	for (int i =0; i<3; i++) tmp += ((k_vec[i]<mesh_num/2.) ? pow(k_vec[i],2) : pow(k_vec[i] - mesh_num, 2.));
	return tmp;
}

double get_k_sq(int mesh_num, int k_vec[3]){
	// returns computation k
	// includes symetry around N/2
	double tmp = 0;
	for (int i =0; i<3; i++) tmp += ((k_vec[i]<mesh_num/2.) ? pow(k_vec[i],2) : pow(k_vec[i] - mesh_num, 2.));
	return tmp;
}

double get_k_sq_phys(int mesh_num, int L, int pos){
	// returns physical k (in units 1/[L])
	// includes symetry around N/2, box size, and 2*PI definition
	return get_k_sq(mesh_num, pos)*pow(2.*PI/L, 2.);
}

int get_pos(int mesh_num, int vec[3]){
	// return position in 1-d long field coresponding
	// to the position in 3-d field vec[]
	// periodicity involved
	int vec_per[3];
	for (int i = 0; i<3; i++) {
		vec_per[i] = get_per(vec[i], mesh_num);
	}
	return vec_per[0]*mesh_num*(mesh_num+2) + vec_per[1]*(mesh_num+2) + vec_per[2];
}

double wgh_sch(double *x, int *y, int mesh_num, const int order){
	// The weighting scheme used to assign values to the mesh points or vice versa
	// Return value of assigment function on mesh point y from particle in x
	double dx, w = 1;
	int y_per;
	double x_per;

	switch (order){
	case 0: {	// NGP: Nearest grid point
				for (int i = 0; i < 3; i++){
					x_per = get_per((x[i] + 0.5), mesh_num);
					y_per = get_per(y[i], mesh_num);
					if ((int)x_per != y_per) w *= 0;
				}
				return w;
	}
	case 1: {	// CIC: Cloud in cells
				for (int i = 0; i < 3; i++){
					x_per = get_per(x[i], mesh_num);
					y_per = get_per(y[i], mesh_num);					
					dx = fmin(fmin(abs(x_per - y_per), x_per + mesh_num - y_per), y_per + mesh_num - x_per);
					if (dx > 1) w *= 0;
					else w *= 1 - dx;
				}
				return w;
	}
	case 2: {	// TSC: Triangular shaped clouds
				for (int i = 0; i < 3; i++){
					x_per = get_per(x[i], mesh_num);
					y_per = get_per(y[i], mesh_num);
					dx = fmin(fmin(abs(x_per - y_per), x_per + mesh_num - y_per), y_per + mesh_num - x_per);
					if (dx > 1.5) w *= 0;
					else if (dx > 0.5) w *= (1.5 - dx)*(1.5 - dx) / 2.0;
					else w *= 3 / 4.0 - dx*dx;
				}
				return w;
	}
	}
	return 0;
}

void assign_fc(double *data, double *x, double &value, int N, const int order, bool asgn){
	// blur 'value' in point 'x' on surrounding mesh points 'data' (asgn == TRUE)
	// or assign values from surroundigs mesh points 'data' to 'value' in point 'x' (asgn == FALSE)
	// NGP: order = 0, Nearest grid point, 1 cell involved
	// CIC: order = 1, Cloud in cells, 1st order, 8 cells involved 
	// TSC: order = 2, Triangular shaped clouds, 2nd order, 27 cells involved
	int y[3], z[3], data_pos;
	for (int i = 0; i < 3; i++){
		z[i] = (int)(x[i] - 0.5*(order - 1));
	}
	for (y[0] = z[0]; y[0] < z[0] + 1 + order; y[0]++){
		for (y[1] = z[1]; y[1] < z[1] + 1 + order; y[1]++){
			for (y[2] = z[2]; y[2] < z[2] + 1 + order; y[2]++){
				/* y[] contains position of current mesh point */
				data_pos = get_pos(N, y);
				if (asgn) data[data_pos] += value * wgh_sch(x, y, N, order);
				else value += data[data_pos] * wgh_sch(x, y, N, order);
			}
		}
	}
}

/* RANDOM FIELD FUNCTIONS */

double genrand_real(){
	return rand() * (1.0 / RAND_MAX);
}

double generateGaussianNoise(const double& mean, const double &stdDev){
	double rn1, rn2, rn;
	do {
        rn1 = -1 + 2*genrand_real();
        rn2 = -1 + 2*genrand_real();
        rn = rn1*rn1 + rn2*rn2;
      } while (rn >= 1.0 || rn == 0);
	
	return stdDev * rn1 * sqrt(-2.0*log(rn)/rn) + mean;
}

void gen_gauss_white_noise_one(int mesh_num, double* rho, int i_min, int i_max, int seed, const double& mean, const double &stdDev){
	printf("Initializing the random field with SEED = %u\n", seed);
	mt19937 mt_rand(seed);
	normal_distribution<double> normal_dist(mean, stdDev);
	
	for(int i=i_min; i< i_max; i++){
		for (int j=0; j < mesh_num; j++){
			rho[i*(mesh_num+2) + j] = normal_dist(mt_rand);
		}
	}	
}

void gen_gauss_white_noise_nt(int mesh_num, double* rho, const double& mean, const double &stdDev, int nt){
	int i_min = 0;
	int i_max = pow(mesh_num, 2.);
	int i_inc = i_max / nt;
	i_max = i_max % nt + i_inc;
	thread* th_dens = new thread[nt];	
	
	unsigned int seed = time(NULL);
	seed = 87654321;
	
	for (int i = 0; i < nt; i++){
		th_dens[i] = thread(gen_gauss_white_noise_one, mesh_num, rho, i_min, i_max, seed+i, mean, stdDev);
		i_min = i_max;
		i_max += i_inc;
	}
	for (int i = 0; i < nt; i++) th_dens[i].join();
	delete[] th_dens;	
}

/* APPROXIMATION FUNCTIONS */

void exp2pot_one(int mesh_num, double* potential, double* expotential, double nu, int i_min, int i_max){
	for(int i=i_min; i< i_max; i++){
		for (int j=0; j < mesh_num; j++){
			potential[i*(mesh_num + 2)+j] = -2.*nu*log(expotential[i*(mesh_num + 2)+j]);
		}
	}
}

void pot2exp_one(int mesh_num, double* potential, double* expotential, double nu, int i_min, int i_max){
	for(int i=i_min; i< i_max; i++){
		for (int j=0; j < mesh_num; j++){
			expotential[i*(mesh_num + 2)+j] = exp(-potential[i*(mesh_num + 2)+j]/(2.*nu));
		}
	}	
}

void pot2exp_nt(int mesh_num, double* potential, double* expotential, double nu, int nt){
	int i_min = 0;
	int i_max = pow(mesh_num, 2.);
	int i_inc = i_max / nt;
	i_max = i_max % nt + i_inc;
	thread* th_dens = new thread[nt];
	for (int i = 0; i < nt; i++){
		th_dens[i] = thread(pot2exp_one, mesh_num, potential, expotential, nu, i_min, i_max);
		i_min = i_max;
		i_max += i_inc;
	}
	for (int i = 0; i < nt; i++) th_dens[i].join();
	delete[] th_dens;
}

void exp2pot_nt(int mesh_num, double* potential, double* expotential, double nu, int nt){
	int i_min = 0;
	int i_max = pow(mesh_num, 2.);
	int i_inc = i_max / nt;
	i_max = i_max % nt + i_inc;
	thread* th_dens = new thread[nt];
	for (int i = 0; i < nt; i++){
		th_dens[i] = thread(exp2pot_one, mesh_num, potential, expotential, nu, i_min, i_max);
		i_min = i_max;
		i_max += i_inc;
	}
	for (int i = 0; i < nt; i++) th_dens[i].join();
	delete[] th_dens;
}

void exp2pot_nt(int mesh_num, double* data, double nu, int nt){
	exp2pot_nt(mesh_num, data, data, nu, nt);
}

void pot2exp_nt(int mesh_num, double* data, double nu, int nt){
	pot2exp_nt(mesh_num, data, data, nu, nt);
}

void gen_rho_dist_k(int mesh_num, int L, double* rho_x, fftw_plan p_F, t_power power_spectrum, double* parameters, double k2_G, int nt){
	/** Generate density distributions \rho(k) in k-space.
	At first, a gaussian white noise (mean = 0, stdDev = 1) is generated,
	then it is convoluted with given power spectrum.
	**/
		
	fftw_complex* rho_k = reinterpret_cast<fftw_complex*>(rho_x);
	
	/* Gaussian white noies */
	printf("Generating gaussian white noise...\n");
	gen_gauss_white_noise_nt(mesh_num, rho_x, 0., 1., nt);
	
	printf("Generating gaussian white noise in k-sapce...\n");
	fftw_execute_dft_r2c(p_F, rho_x, rho_k);
	normalize_FFT_FORWARD_nt(mesh_num, rho_k, nt);
	
	/* Generate density distributions with given power spectrum */
	printf("Generating density distributions with given power spectrum...\n");
	gen_rho_k_nt(mesh_num, L, rho_k, power_spectrum, parameters, k2_G, nt);
}

void gen_displ_k_w_one(int mesh_num, double** displ_vec, double* potential, int i_min, int i_max){
	// optimalization for CIC
	int k_vec[3];
	double tmp;
	double potential_tmp[2];
	const int n_max = 1;
	double k_n[3];
	double U2, U_n, G_n, k2n;

	for(int i=i_min; i <i_max;i++){
		get_coord_fftw_complex(mesh_num, i, k_vec);		
		G_n = 0;
		U2 = 1;
		for (int j = 0; j < 3; j++){
			if(k_vec[j] > mesh_num/2.) k_vec[j] -= mesh_num;
			U2 *= 1./3.*(1+2*cos(PI*k_vec[j]/mesh_num));
		}
		for (int n1 = -n_max; n1 < n_max + 1; n1++){
			k_n[0] = 2 * PI / mesh_num*k_vec[0] + 2 * PI*n1;
			for (int n2 = -n_max; n2 < n_max + 1; n2++){
				k_n[1] = 2 * PI / mesh_num*k_vec[1] + 2 * PI*n2;
				for (int n3 = -n_max; n3 < n_max + 1; n3++){
					k_n[2] = 2 * PI / mesh_num*k_vec[2] + 2 * PI*n3;
					U_n = 1.;
					k2n = 0;
					for(int j=0; j<3; j++){
						if (k_n[j] != 0) U_n *= sin(k_n[j] / 2.) / (k_n[j] / 2.);
						k2n += pow(k_n[j],2);
					}
					if (k2n != 0){
						for(int j=0; j<3; j++){										
							G_n += 2 * PI / mesh_num*k_vec[j]* // D(k)
							k_n[j]/k2n* // R(k_n)
							pow(U_n, 2.); // W(k) for CIC
						}
					}
				}
			}
		}
		if ((G_n != G_n) || (U2 != U2)) printf("Gn = %f\tU2 = %f, k = (%i, %i, %i) \n", G_n, U2, k_vec[0], k_vec[1], k_vec[2]);
		potential_tmp[0] = potential[2*i]; // prevent overwriting if displ_vec[i] == potential
		potential_tmp[1] = potential[2*i+1]; // prevent overwriting if displ_vec[i] == potential
		for(int j=0; j<3;j++){
			if (k_vec[j] == mesh_num/2.) tmp = 0;
			else tmp = k_vec[j];			
	
			displ_vec[j][2*i] = tmp*potential_tmp[1]*(2.*PI/mesh_num)*G_n/U2; // 2*PI/N comes from derivative WITH RESPECT to the mesh coordinates
			displ_vec[j][2*i+1] = -tmp*potential_tmp[0]*(2.*PI/mesh_num)*G_n/U2; // 2*PI/N comes from derivative WITH RESPECT to the mesh coordinates	
		}
	}
}

void gen_displ_k_one(int mesh_num, double** displ_vec, double* potential, int i_min, int i_max){
	
	int k_vec[3];
	double tmp;
	double potential_tmp[2];

	for(int i=i_min; i <i_max;i++){
		get_coord_fftw_complex(mesh_num, i, k_vec);		
		potential_tmp[0] = potential[2*i]; // prevent overwriting if displ_vec[i] == potential
		potential_tmp[1] = potential[2*i+1]; // prevent overwriting if displ_vec[i] == potential
		for(int j=0; j<3;j++){
			if (k_vec[j] < mesh_num/2.) tmp = k_vec[j];
			else if (k_vec[j] > mesh_num/2.) tmp = k_vec[j]-mesh_num;
			else tmp = 0;			
	
			displ_vec[j][2*i] = tmp*potential_tmp[1]*(2.*PI/mesh_num); // 2*PI/N comes from derivative WITH RESPECT to the mesh coordinates
			displ_vec[j][2*i+1] = -tmp*potential_tmp[0]*(2.*PI/mesh_num); // 2*PI/N comes from derivative WITH RESPECT to the mesh coordinates	
		}
	}
}

void gen_displ_k_w_nt(int mesh_num, double** displ_vec, double* potential, int nt){
	int i_min = 0;
	int i_max = pow(mesh_num, 2)*(mesh_num/2 + 1);	
	int i_inc = i_max / nt;
	i_max = i_max % nt + i_inc;
	thread* th_dens = new thread[nt];
	for (int i = 0; i < nt; i++){
		th_dens[i] = thread(gen_displ_k_w_one, mesh_num, displ_vec, potential, i_min, i_max);
		i_min = i_max;
		i_max += i_inc;
	}
	for (int i = 0; i < nt; i++) th_dens[i].join();
	delete[] th_dens;
}

void gen_displ_k_nt(int mesh_num, double** displ_vec, double* potential, int nt){
	int i_min = 0;
	int i_max = pow(mesh_num, 2)*(mesh_num/2 + 1);	
	int i_inc = i_max / nt;
	i_max = i_max % nt + i_inc;
	thread* th_dens = new thread[nt];
	for (int i = 0; i < nt; i++){
		th_dens[i] = thread(gen_displ_k_one, mesh_num, displ_vec, potential, i_min, i_max);
		i_min = i_max;
		i_max += i_inc;
	}
	for (int i = 0; i < nt; i++) th_dens[i].join();
	delete[] th_dens;
}

/*void gen_no_displ_pos_vel_one(int par_num, int Ng, double*** field, int i_min, int i_max){
	int p_vec[3];
	for(int i=i_min; i< i_max; i++){
		p_vec[0] = (i / par_num) * Ng;
		p_vec[1] = (i % par_num) * Ng;
		for (int j=0; j < par_num; j++){
			p_vec[2] = j * Ng; // initial coordinates for i-th particle
			for(int k=0; k<3;k++){
				field[0][k][i*(par_num + 2)+j] = p_vec[k];
				field[1][k][i*(par_num + 2)+j] = 0;
			}
		}
	}
}
*/

/* void gen_no_displ_pos_vel_nt(int par_num, int Ng, double*** field, int nt){
	int i_min = 0;
	int i_max = pow(par_num, 2.);
	int i_inc = i_max / nt;
	i_max = i_max % nt + i_inc;
	thread* th_dens = new thread[nt];
	for (int i = 0; i < nt; i++){
		th_dens[i] = thread(gen_no_displ_pos_vel_one, par_num, Ng, field, i_min, i_max);
		i_min = i_max;
		i_max += i_inc;
	}
	for (int i = 0; i < nt; i++) th_dens[i].join();
	// printf(">>> All threads in routine <gen_no_displ_pos_vel_nt> joined.\n");
	delete[] th_dens;	
}
*/

void gen_no_displ_pos_one(int par_num, int Ng, double** displ_vec, int i_min, int i_max){
	int p_vec[3];
	for(int i=i_min; i< i_max; i++){
		p_vec[0] = (i / par_num) * Ng;
		p_vec[1] = (i % par_num) * Ng;
		for (int j=0; j < par_num; j++){
			p_vec[2] = j * Ng; // initial coordinates for i-th particle
			for(int k=0; k<3;k++){
				displ_vec[k][i*(par_num + 2)+j] = p_vec[k] ;
			}
		}
	}
}

void gen_no_displ_pos_nt(int par_num, int Ng, double** displ_vec, int nt){
	int i_min = 0;
	int i_max = pow(par_num, 2.);
	int i_inc = i_max / nt;
	i_max = i_max % nt + i_inc;
	thread* th_dens = new thread[nt];
	for (int i = 0; i < nt; i++){
		th_dens[i] = thread(gen_no_displ_pos_one, par_num, Ng, displ_vec, i_min, i_max);
		i_min = i_max;
		i_max += i_inc;
	}
	for (int i = 0; i < nt; i++) th_dens[i].join();
	delete[] th_dens;	
}

void gen_init_expot(int mesh_num, double* potential, double* expotential_0, double nu, bool conv_mod, fftw_plan p_F, int nt){
	if (conv_mod){ // convolution using FFT
		/* Computing initial expotential in k-space */
		printf("Computing initial expotential in k-space...\n");
		pot2exp_nt(mesh_num, potential, expotential_0, nu, nt);	
		fftw_execute_dft_r2c(p_F, expotential_0, reinterpret_cast<fftw_complex*>(expotential_0));
		// we do not care about normalization because of logarithmic derivative	
	}else{ // convolution using direct sum
		printf("Storing initial potenital in q-space...\n");
		copy_nt(mesh_num, potential, expotential_0, nt);
		/* Computing initial expotential in q-space */
	//	printf("Computing initial expotential in q-space...\n");
	//	pot2exp_nt(mesh_num, potential, expotential_0, nu, nt);
	}
}


void convolution_y1_one(int mesh_num, double* potential, double* expotential_0, double nu, double b, int i_min, int i_max){
	// multi-thread index is y3
	// compute f1 (x1, y2, y3)
	double exponent;
	int pos;
	double sum;
	
	for (int x1 = 0; x1 < mesh_num; x1++){
		for (int y2 = 0; y2 < mesh_num; y2++){
			for (int y3 = i_min; y3 < i_max; y3++){
				// summation over y1
				sum = 0;
				for (int y1 = 0; y1 < mesh_num; y1++){
					pos = y1*mesh_num*(mesh_num+2)+y2*(mesh_num+2)+y3;
					exponent=expotential_0[pos]+pow(x1-y1, 2.)/(2.*b);
					exponent = -exponent/(2.*nu);
					sum +=exp(exponent);
					if (isfinite(sum) == 0) {printf("Error while chcecking for NAN or INF. r = %i, potential = %f, b = %f, exponent = %f, +sum = %f\n",
						abs(x1-y1), expotential_0[pos], b, exponent, exp(-exponent/(2.*nu))); return;}
				}
				pos = x1*mesh_num*(mesh_num+2)+y2*(mesh_num+2)+y3;
				potential[pos] = sum; // potential is now f1
			}
		}
	}
}

void convolution_y1_nt(int mesh_num, double* potential, double* expotential_0, double nu, double b, int nt){
	int i_min = 0;
	int i_max = mesh_num;
	int i_inc = i_max / nt;
	i_max = i_max % nt + i_inc;
	thread* th_dens = new thread[nt];
	for (int i = 0; i < nt; i++){
		th_dens[i] = thread(convolution_y1_one, mesh_num, potential, expotential_0, nu, b, i_min, i_max);
		i_min = i_max;
		i_max += i_inc;
	}
	for (int i = 0; i < nt; i++) th_dens[i].join();
	delete[] th_dens;	
}

void convolution_y2_one(int mesh_num, double* potential, double* gaussian, double nu, double b, int i_min, int i_max){
	// multi-thread index is x1
	// compute f1 (x1, y2, y3)
	int pos;
	double sum;
	double* sum_aux = new double[mesh_num];
	
	// compute f2 (x1, x2, y3)
	for (int x1 = i_min; x1 < i_max; x1++){
		for (int y3 = 0; y3 < mesh_num; y3++){
			for (int x2 = 0; x2 < mesh_num; x2++){
				// summation over y2
				sum = 0;
				for (int y2 = 0; y2 < mesh_num; y2++){
					pos = x1*mesh_num*(mesh_num+2)+y2*(mesh_num+2)+y3;
					sum += potential[pos]*gaussian[abs(x2-y2)];
					if (isfinite(sum) == 0) {printf("Error while chcecking for NAN or INF. r = %i, potential = %f, gaussian = %f, +sum = %f\n",
						abs(x2-y2), potential[pos], gaussian[abs(x2-y2)], potential[pos]*gaussian[abs(x2-y2)]); return;}
				}
				sum_aux[x2] = sum;
			}

			for (int x2 = 0; x2 < mesh_num; x2++){
				pos = x1*mesh_num*(mesh_num+2)+x2*(mesh_num+2)+y3;
				potential[pos] = sum_aux[x2]; // potential is now f2
			}
		}
	}
	delete[] sum_aux;
}

void convolution_y2_nt(int mesh_num, double* potential, double* gaussian, double nu, double b, int nt){
	int i_min = 0;
	int i_max = mesh_num;
	int i_inc = i_max / nt;
	i_max = i_max % nt + i_inc;
	thread* th_dens = new thread[nt];
	for (int i = 0; i < nt; i++){
		th_dens[i] = thread(convolution_y2_one, mesh_num, potential, gaussian, nu, b, i_min, i_max);
		i_min = i_max;
		i_max += i_inc;
	}
	for (int i = 0; i < nt; i++) th_dens[i].join();
	delete[] th_dens;	
}

void convolution_y3_one(int mesh_num, double* potential, double* gaussian, double nu, double b, int i_min, int i_max){
	// multi-thread index is x1
	// compute f1 (x1, y2, y3)
	int pos;
	double sum;
	double* sum_aux = new double[mesh_num];

	// compute f3 (x1, x2, x3) == expotential(x, b)
	for (int x1 = i_min; x1 < i_max; x1++){
		for (int x2 = 0; x2 < mesh_num; x2++){
			for (int x3 = 0; x3 < mesh_num; x3++){
				// summation over y3
				sum = 0;
				for (int y3 = 0; y3 < mesh_num; y3++){
					pos = x1*mesh_num*(mesh_num+2)+x2*(mesh_num+2)+y3;
					sum += potential[pos]*gaussian[abs(x3-y3)];
					if (isfinite(sum) == 0) {printf("Error while chcecking for NAN or INF. r = %i, potential = %f, gaussian = %f, +sum = %f\n",
						abs(x3-y3), potential[pos], gaussian[abs(x3-y3)], potential[pos]*gaussian[abs(x3-y3)]); return;}
				}
				sum_aux[x3] = sum;
			}
			for (int x3 = 0; x3 < mesh_num; x3++){
				pos = x1*mesh_num*(mesh_num+2)+x2*(mesh_num+2)+x3;
				potential[pos] = sum_aux[x3]; // potential is now f3
			}
		}
	}
	delete[] sum_aux;
}

void convolution_y3_nt(int mesh_num, double* potential, double* gaussian, double nu, double b, int nt){
	int i_min = 0;
	int i_max = mesh_num;
	int i_inc = i_max / nt;
	i_max = i_max % nt + i_inc;
	thread* th_dens = new thread[nt];
	for (int i = 0; i < nt; i++){
		th_dens[i] = thread(convolution_y3_one, mesh_num, potential, gaussian, nu, b, i_min, i_max);
		i_min = i_max;
		i_max += i_inc;
	}
	for (int i = 0; i < nt; i++) th_dens[i].join();
	delete[] th_dens;	
}

void gen_expot(int mesh_num, double* potential, double* expotential_0, double nu, double b, bool conv_mod, fftw_plan p_B, int nt){
	if (conv_mod){
		/* Computing convolution using FFT */
		printf("Computing expotential in k-space...\n");
		double k2;
		for(int i=0; i < pow(mesh_num, 2)*(mesh_num/2 + 1);i++){
			k2 = get_k_sq(mesh_num, i); // get k2 ncluding symetry around N/2
			potential[2*i] = expotential_0[2*i]*exp(-4.*PI*PI*b*nu*k2/(mesh_num*mesh_num)); // DEFINITION of nu_phys = (L/N)^2 * v_comp
			potential[2*i+1] = expotential_0[2*i+1]*exp(-4.*PI*PI*b*nu*k2/(mesh_num*mesh_num));
			
		}
		fftw_execute_dft_c2r(p_B, reinterpret_cast<fftw_complex*>(potential), potential);
	}else{
		/* Computing convolution using direct sum */
		printf("Computing expotential in q-space...\n");
		/*
		f(x1, x2, x3) = \int dy^3 { g(y1, y2, y3) * h(x1 - y1) * h(x2 - y2) * h(x3 - y3)}
		..
		f1 (x1, y2, y3) = \int dy1 { g  (y1, y2, y3) * h(x1 - y1)}	:: N^3 sums of length N
		f2 (x1, x2, y3) = \int dy2 { f1 (x1, y2, y3) * h(x2 - y2)}	:: N^3 sums of length N
		f3 (x1, x2, x3) = \int dy3 { f2 (x1, x2, y3) * h(x3 - y3)}	:: N^3 sums of length N
		*/
		
		// store values of exponential - every convolution uses the same exp(-r^2/4bv)
		double* gaussian = new double[mesh_num];
		for (int i = 0; i < mesh_num; i++){
			gaussian[i]=exp(-i*i/(4.*b*nu));
		}

		convolution_y1_nt(mesh_num, potential, expotential_0, nu, b, nt);
		convolution_y2_nt(mesh_num, potential, gaussian, nu, b, nt);
		convolution_y3_nt(mesh_num, potential, gaussian, nu, b, nt);
		
		delete [] gaussian;
	}
}

void fin_diff(int mesh_num, double** vec_field, double* potential){
	int i_m, i_p, j_m, j_p, k_m, k_p;
	
	for (int i = 0; i < mesh_num; i++){
		i_m = (i - 1 + mesh_num) % mesh_num;
		i_p = (i + 1) % mesh_num;
		for (int j = 0; j < mesh_num; j++){
			j_m = (j - 1 + mesh_num) % mesh_num;
			j_p = (j + 1) % mesh_num;
			for (int k = 0; k < mesh_num; k++){
				k_m = (k - 1 + mesh_num) % mesh_num;
				k_p = (k + 1) % mesh_num;
				
				vec_field[2][i*(mesh_num +2)*mesh_num+j*(mesh_num+2) + k] = -(
					potential[i*(mesh_num+2)*mesh_num+j*(mesh_num+2) + k_p] -
					potential[i*(mesh_num+2)*mesh_num+j*(mesh_num+2) + k_m]) / 2.;
					
				vec_field[1][i*(mesh_num +2)*mesh_num+j*(mesh_num+2) + k]= -(
					potential[i*(mesh_num+2)*mesh_num+j_p*(mesh_num+2)+ k] -
					potential[i*(mesh_num+2)*mesh_num+j_m*(mesh_num+2)+ k]) / 2.;
					
				vec_field[0][i*(mesh_num +2)*mesh_num+j*(mesh_num + 2)+ k]= -(
					potential[i_p*(mesh_num+2)*mesh_num+j*(mesh_num+2)+ k] -
					potential[i_m*(mesh_num+2)*mesh_num+j*(mesh_num+2)+ k]) / 2.;
			}
		}
	}
}

/* MULTI-THREAD FUNCTIONS */

void init_dens_field_one(int mesh_num, double *delta, int i_min, int i_max){
	for(int i=i_min; i<i_max; i++){
		for (int j=0; j < mesh_num; j++){
			delta[i*(mesh_num + 2)+j] = -1.;
		}
	}	
}

void init_dens_field(int mesh_num, double *delta, int nt){
	int i_min = 0;
	int i_max = pow(mesh_num, 2.);
	int i_inc = i_max / nt;
	i_max = i_max % nt + i_inc;
	thread* th_dens = new thread[nt];
	for (int i = 0; i < nt; i++){
		th_dens[i] = thread(init_dens_field_one, mesh_num, delta, i_min, i_max);
		i_min = i_max;
		i_max += i_inc;
	}
	for (int i = 0; i < nt; i++) th_dens[i].join();
	delete[] th_dens;
}

void gen_rho_k_one(int mesh_num, int L, fftw_complex* rho, int i_min, int i_max, t_power power_spectrum, double* parameters){
	double k2;
	for(int i=i_min; i < i_max;i++){
		k2 = get_k_sq_phys(mesh_num, L, i);
	// Power Spectrum units [L/N]^3
	//	rho[i][0] *= sqrt(power_spectrum(sqrt(k2), parameters) * pow((double)mesh_num / L, 3.));
	//	rho[i][1] *= sqrt(power_spectrum(sqrt(k2), parameters) * pow((double)mesh_num / L, 3.));
	// Power Spectrum units dimensionless
		rho[i][0] *= sqrt(power_spectrum(sqrt(k2), parameters));
		rho[i][1] *= sqrt(power_spectrum(sqrt(k2), parameters));
	}	
}

void gen_rho_k_one_s(int mesh_num, int L, fftw_complex* rho, int i_min, int i_max, t_power power_spectrum, double* parameters, double k2_G){
	double k2;
	printf("Gen rho(k) with smoothing k = %.2f\n", k2_G);
	for(int i=i_min; i < i_max;i++){
		k2 = get_k_sq_phys(mesh_num, L, i);
	// Power Spectrum units [L/N]^3
	//	rho[i][0] *= sqrt(power_spectrum(sqrt(k2), parameters) * pow((double)mesh_num / L, 3.)*exp(-k2/(2.*k2_G)));
	//	rho[i][1] *= sqrt(power_spectrum(sqrt(k2), parameters) * pow((double)mesh_num / L, 3.)*exp(-k2/(2.*k2_G)));
	// Power Spectrum units dimensionless
		rho[i][0] *= sqrt(power_spectrum(sqrt(k2), parameters)*exp(-k2/(2.*k2_G)));
		rho[i][1] *= sqrt(power_spectrum(sqrt(k2), parameters)*exp(-k2/(2.*k2_G)));	
	}	
}


void gen_rho_k_nt(int mesh_num, int L, fftw_complex* rho, t_power power_spectrum, double* parameters, double k2_G, int nt){
	int i_min = 0;
	int i_max = pow(mesh_num, 2)*(mesh_num/2 + 1);	
	int i_inc = i_max / nt;
	i_max = i_max % nt + i_inc;
	thread* th_dens = new thread[nt];
	for (int i = 0; i < nt; i++){
		if (k2_G == 0) th_dens[i] = thread(gen_rho_k_one, mesh_num, L, rho, i_min, i_max, power_spectrum, parameters);
		else th_dens[i] = thread(gen_rho_k_one_s, mesh_num, L, rho, i_min, i_max, power_spectrum, parameters, k2_G);
		i_min = i_max;
		i_max += i_inc;
	}
	for (int i = 0; i < nt; i++) th_dens[i].join();
	delete[] th_dens;
}

void upd_pos_one(int par_num, int mesh_num, double** par_pos, double** vec_field, double db, int order, int i_min, int i_max){
	double x[3];
	double v;
	for(int i=i_min; i< i_max; i++){
		for (int j=0; j < par_num; j++){
			for(int k=0; k<3;k++) x[k]=par_pos[k][i*(par_num + 2)+j];
			for(int k=0; k<3;k++){
				v = 0;
				assign_fc(vec_field[k], x, v, mesh_num, order, false);
				par_pos[k][i*(par_num + 2)+j] = get_per(par_pos[k][i*(par_num + 2)+j] + v*db, mesh_num);
			}
		}
	}
}

/* void upd_pos_vel0_one(int par_num, int mesh_num, double** par_pos, double** vec_field, double db, int order, int i_min, int i_max){
	double x[3];
	double v;
	for(int i=i_min; i< i_max; i++){
		x[0] = i / par_num;
		x[1] = i % par_num;
		for (int j=0; j < par_num; j++){
			x[2] = j; // initial coordinates for i-th particle
			for(int k=0; k<3;k++){
				v = 0;
				assign_fc(vec_field[k], x, v, mesh_num, order, false);
				par_pos[k][i*(par_num + 2)+j] = get_per(par_pos[k][i*(par_num + 2)+j] + v*db, mesh_num);
			}
		}
	}
}
*/

void upd_pos_mid_one(int par_num, int mesh_num, double*** par_pos, double** vec_field, double db, int order, int i_min, int i_max){
	// modified MIDPOINT method, need to store previous velocitites
	// par_pos[0][j][i] :: j-th position coordinate of i-th particle
	// par_pos[1][j][i] :: j-th velocity coordinate (previous) of i-th particle
	// assuming that vec_field is evaluated at b+db/2
	double x[3];
	double v;
	for(int i=i_min; i< i_max; i++){
		for (int j=0; j < par_num; j++){
			for(int k=0; k<3;k++) x[k]=par_pos[0][k][i*(par_num + 2)+j] + db/2.*par_pos[1][k][i*(par_num + 2)+j];
			for(int k=0; k<3;k++){
				v = 0;
				assign_fc(vec_field[k], x, v, mesh_num, order, false);
				par_pos[0][k][i*(par_num + 2)+j] = get_per(par_pos[0][k][i*(par_num + 2)+j] + v*db, mesh_num);
				par_pos[1][k][i*(par_num + 2)+j] = v;
			}
		}
	}
}

void upd_pos_leapfrog_one(int par_num, int mesh_num, double*** par_pos, double** force_field, double b_half, double db, int order, int i_min, int i_max){
	double x[3]; // position at half step
	double v, f;
	for(int i=i_min; i< i_max; i++){
		for (int j=0; j < par_num; j++){
			for(int k=0; k<3;k++) x[k]=par_pos[0][k][i*(par_num + 2)+j] + db/2.*par_pos[1][k][i*(par_num + 2)+j];
		//	for(int k=0; k<3;k++) x[k]=par_pos[0][k][i*(par_num + 2)+j];
			for(int k=0; k<3;k++){
				f = 0;
				assign_fc(force_field[k], x, f, mesh_num, order, false);				
				v = par_pos[1][k][i*(par_num + 2)+j];
				f = -3/(2.*b_half)*(v-f);
				v = v + f*db;
				par_pos[0][k][i*(par_num + 2)+j] = get_per(x[k] + v*db/2., mesh_num);
		//		par_pos[0][k][i*(par_num + 2)+j] = get_per(x[k] + v*db, mesh_num);
				par_pos[1][k][i*(par_num + 2)+j] = v;
			}
		}
	}
}

void upd_pos_mff_one(int par_num, int mesh_num, double** par_pos, double** force_field, double db, int order, int i_min, int i_max){
	double x[3], v_x[3];
	double v, v_n;
	double check, rel;
	for(int i=i_min; i< i_max; i++){
		for (int j=0; j < par_num; j++){
			for(int k=0; k<3;k++) x[k]=par_pos[k][i*(par_num + 2)+j];
			v = par_pos[3][i*(par_num + 2)+j];
			v_n = 0;
			for(int k=0; k<3;k++){
				v_x[k] = 0;
				assign_fc(force_field[k], x, v_x[k], mesh_num, order, false);
				v_n += pow (v_x[k], 2);
			}
			v_n = sqrt(v_n);
			check = 0;
			for(int k=0; k<3;k++){
				check += pow(v*v_x[k]/v_n, 2.);
				par_pos[k][i*(par_num + 2)+j] = get_per(x[k] + v*v_x[k]/v_n*db, mesh_num);
			}
			check = sqrt(check);
			rel = abs(v-check)/v;
			if (rel > 1E-10){
				printf("Original v = %.12f, but new v = %.12f. Relative error = %.12f\nDirection = (", v, check, abs(v-check)/v);
				for (int p = 0; p<3; p++) printf(" %f ", v_x[p]);
				printf(")\tNormalization = %f\n\n", v_n);
			}
		}
	}
}

void upd_pos_nt(int par_num, int mesh_num, double** par_pos, double** vec_field, double db, int order, int nt){
	int i_min = 0;
	int i_max = pow(par_num, 2.);
	int i_inc = i_max / nt;
	i_max = i_max % nt + i_inc;
	thread* th_dens = new thread[nt];
	for (int i = 0; i < nt; i++){
		th_dens[i] = thread(upd_pos_one, par_num, mesh_num, par_pos, vec_field, db, order, i_min, i_max);
		i_min = i_max;
		i_max += i_inc;
	}
	for (int i = 0; i < nt; i++) th_dens[i].join();
	delete[] th_dens;	
}

/* void upd_pos_vel0_nt(int par_num, int mesh_num, double** par_pos, double** vec_field, double db, int order, int nt){
	// updating positions based on initial position of particle, i.e. ZA
	int i_min = 0;
	int i_max = pow(par_num, 2.);
	int i_inc = i_max / nt;
	i_max = i_max % nt + i_inc;
	thread* th_dens = new thread[nt];
	for (int i = 0; i < nt; i++){
		th_dens[i] = thread(upd_pos_vel0_one, par_num, mesh_num, par_pos, vec_field, db, order, i_min, i_max);
		i_min = i_max;
		i_max += i_inc;
	}
	for (int i = 0; i < nt; i++) th_dens[i].join();
	delete[] th_dens;	
}
*/

void upd_pos_mid_nt(int par_num, int mesh_num, double*** par_pos, double** vec_field, double db, int order, int nt){
	int i_min = 0;
	int i_max = pow(mesh_num, 2.);
	int i_inc = i_max / nt;
	i_max = i_max % nt + i_inc;
	thread* th_dens = new thread[nt];
	for (int i = 0; i < nt; i++){
		th_dens[i] = thread(upd_pos_mid_one, par_num, mesh_num, par_pos, vec_field, db, order, i_min, i_max);
		i_min = i_max;
		i_max += i_inc;
	}
	for (int i = 0; i < nt; i++) th_dens[i].join();
	// printf(">>> All threads in routine <upd_pos_mid_nt> joined.\n");
	delete[] th_dens;	
}

void upd_pos_leapfrog_nt(int par_num, int mesh_num, double*** par_pos, double** vec_field, double b_half, double db, int order, int nt){
	int i_min = 0;
	int i_max = pow(par_num, 2.);
	int i_inc = i_max / nt;
	i_max = i_max % nt + i_inc;
	thread* th_dens = new thread[nt];
	for (int i = 0; i < nt; i++){
		th_dens[i] = thread(upd_pos_leapfrog_one, par_num, mesh_num, par_pos, vec_field, b_half, db, order, i_min, i_max);
		i_min = i_max;
		i_max += i_inc;
	}
	for (int i = 0; i < nt; i++) th_dens[i].join();
	delete[] th_dens;	
}

void upd_pos_mff_nt(int par_num, int mesh_num, double** par_pos, double** vec_field, double db, int order, int nt){
	int i_min = 0;
	int i_max = pow(par_num, 2.);
	int i_inc = i_max / nt;
	i_max = i_max % nt + i_inc;
	thread* th_dens = new thread[nt];
	for (int i = 0; i < nt; i++){
		th_dens[i] = thread(upd_pos_mff_one, par_num, mesh_num, par_pos, vec_field, db, order, i_min, i_max);
		i_min = i_max;
		i_max += i_inc;
	}
	for (int i = 0; i < nt; i++) th_dens[i].join();
	delete[] th_dens;	
}

void gen_pot_k_one(int mesh_num, double* potential, int i_min, int i_max){
	int k_vec[3];
	double k2;	
	for(int i=i_min; i < i_max;i++){				
		get_coord_fftw_complex(mesh_num, i, k_vec);
		k2 = get_k_sq(mesh_num, k_vec); // get k2 ncluding symetry around N/2
		if (k2 == 0){
			potential[2*i] = 0;
			potential[2*i+1] = 0;
		} else{
			potential[2*i] /= -(k2*pow(2.*PI/mesh_num, 2.));
			potential[2*i+1] /= -(k2*pow(2.*PI/mesh_num, 2.));
		}
	}
}

void gen_pot_k_nt(int mesh_num, double* potential, int nt){
	int i_min = 0;
	int i_max = pow(mesh_num, 2)*(mesh_num/2 + 1);	
	int i_inc = i_max / nt;
	i_max = i_max % nt + i_inc;
	thread* th_dens = new thread[nt];
	for (int i = 0; i < nt; i++){
		th_dens[i] = thread(gen_pot_k_one, mesh_num, potential, i_min, i_max);
		i_min = i_max;
		i_max += i_inc;
	}
	for (int i = 0; i < nt; i++) th_dens[i].join();
	delete[] th_dens;
}


void gen_displ_pos_ng_one(int par_num, int Ng, double** vec_field, double** par_pos, double b, int i_min, int i_max){
	int p_vec[3];
	double s_q, p_q, x_q;
	const int mesh_num = Ng*par_num;
	int pos;
	
	for(int i=i_min; i< i_max; i++){
		p_vec[0] = (i / par_num) * Ng;
		p_vec[1] = (i % par_num) * Ng;
		for (int j=0; j < par_num; j++){
			p_vec[2] = j * Ng; // initial coordinates for i-th particle
			pos = get_pos(mesh_num, p_vec); // pos on mesh
			for(int k=0; k<3;k++){
				s_q = vec_field[k][pos]*b;
				p_q = p_vec[k];
				x_q = get_per(p_q + s_q, mesh_num); // final position including periodicity
				par_pos[k][i*(par_num + 2)+j] = x_q;
			}
		}
	}
}

void gen_displ_pos_ng_nt(int par_num, int Ng, double** vec_field, double** par_pos, double b, int nt){
	int i_min = 0;
	int i_max = pow(par_num, 2.);
	int i_inc = i_max / nt;
	i_max = i_max % nt + i_inc;
	thread* th_dens = new thread[nt];
	for (int i = 0; i < nt; i++){
		th_dens[i] = thread(gen_displ_pos_ng_one, par_num, Ng, vec_field, par_pos, b, i_min, i_max);
		i_min = i_max;
		i_max += i_inc;
	}
	for (int i = 0; i < nt; i++) th_dens[i].join();
	delete[] th_dens;	
}

void gen_init_con_one(int par_num, int Ng, double** vec_field, double*** par_pos, int i_min, int i_max){
	int p_vec[3];
	const int mesh_num = Ng*par_num;
	int pos;
	
	for(int i=i_min; i< i_max; i++){
		p_vec[0] = (i / par_num) * Ng;
		p_vec[1] = (i % par_num) * Ng;
		for (int j=0; j < par_num; j++){
			p_vec[2] = j * Ng; // initial coordinates for i-th particle	
			pos = get_pos(mesh_num, p_vec);			
			for(int k=0; k<3;k++){
				par_pos[0][k][i*(par_num + 2)+j] = p_vec[k];
				par_pos[1][k][i*(par_num + 2)+j] = vec_field[k][pos];
			}
		}
	}
}

void gen_init_con_nt(int par_num, int Ng, double** vec_field, double*** par_pos, int nt){
	int i_min = 0;
	int i_max = pow(par_num, 2.);
	int i_inc = i_max / nt;
	i_max = i_max % nt + i_inc;
	thread* th_dens = new thread[nt];
	for (int i = 0; i < nt; i++){
		th_dens[i] = thread(gen_init_con_one, par_num, Ng, vec_field, par_pos, i_min, i_max);
		i_min = i_max;
		i_max += i_inc;
	}
	for (int i = 0; i < nt; i++) th_dens[i].join();
	delete[] th_dens;	
}

void gen_mff_ic_one(int par_num, int Ng, double** vec_field, double** par_pos, int i_min, int i_max){
	int p_vec[3];
	const int mesh_num = Ng*par_num;
	double v;
	int pos;
	for(int i=i_min; i< i_max; i++){
		p_vec[0] = (i / par_num) * Ng;
		p_vec[1] = (i % par_num) * Ng;
		for (int j=0; j < par_num; j++){
			p_vec[2] = j * Ng; // initial coordinates for i-th particle
			v = 0;
			pos = get_pos(mesh_num, p_vec);
			for(int k=0; k<3;k++){
				par_pos[k][i*(par_num + 2)+j] = p_vec[k];
				v += pow(vec_field[k][pos], 2.);
			}
			par_pos[3][i*(par_num + 2)+j] = sqrt(v);
		}
	}
}

void gen_mff_ic_nt(int par_num, int Ng, double** vec_field, double** par_pos, int nt){
	int i_min = 0;
	int i_max = pow(par_num, 2.);
	int i_inc = i_max / nt;
	i_max = i_max % nt + i_inc;
	thread* th_dens = new thread[nt];
	for (int i = 0; i < nt; i++){
		th_dens[i] = thread(gen_mff_ic_one, par_num, Ng, vec_field, par_pos, i_min, i_max);
		i_min = i_max;
		i_max += i_inc;
	}
	for (int i = 0; i < nt; i++) th_dens[i].join();
	delete[] th_dens;	
}

/* ALLOCATION FUNCTIONS */

void** alloc_zeldovich(int mesh_num, int par_num, int bin_num, int nt){
	try{
		void** alloc_zel = new void* [7];
		
		/* Allocate long arrays */
		const long mem_double =
			3*par_num*par_num*(par_num+2) + // particle positions
			3*mesh_num*mesh_num*(mesh_num+2) + // velocity field
			mesh_num*mesh_num*(mesh_num+2); // work space for power spectrum
		
		double* alloc_double = (double*) fftw_malloc (sizeof(double)*mem_double);
		int pos = 0;
		
		double** vel_field = new double* [3];
		for (int i = 0; i < 3; i++){
			vel_field[i] = alloc_double + pos;
			pos += mesh_num*mesh_num*(mesh_num+2);
		}
		alloc_zel[0] = (void*) vel_field;
		
		double** par_pos = new double* [3];
		for (int i = 0; i < 3; i++){
			par_pos[i] = alloc_double + pos;
			pos += par_num*par_num*(par_num+2);
		}
		alloc_zel[1] = (void*) par_pos;
		
		/* Allocate small arrays */
		double* power_aux = alloc_double + pos;
		pos += mesh_num*mesh_num*(mesh_num+2);
		alloc_zel[2] = (void*) power_aux;
		
		fftw_complex* alloc_fftw_c = (fftw_complex*) fftw_malloc (sizeof(fftw_complex)*2*bin_num);
		alloc_zel[3]=(void*)alloc_fftw_c;
		alloc_zel[4]=(void*)(alloc_fftw_c+bin_num);
		
		if (fftw_init_threads() == 0){
			printf("Errors while initializing multi-thread!\n");
			return NULL;
		}
		fftw_plan_with_nthreads(nt); // FFTW_ESTIMATE FFTW_MEASURE
		fftw_plan p_F = fftw_plan_dft_r2c_3d(mesh_num, mesh_num, mesh_num, alloc_double, reinterpret_cast<fftw_complex*>(alloc_double), FFTW_ESTIMATE);
		fftw_plan p_B = fftw_plan_dft_c2r_3d(mesh_num, mesh_num, mesh_num, reinterpret_cast<fftw_complex*>(alloc_double), alloc_double, FFTW_ESTIMATE);
		alloc_zel[5]=(void*)p_F;
		alloc_zel[6]=(void*)p_B;
		
		printf("Allocated %s of memory.\n", humanSize(sizeof(double)*mem_double));
		return alloc_zel;
	}
	catch(bad_alloc&){
		printf("ERROR! Allocation of memory failed!\n");
		return NULL;
	}
	return NULL;
}

void dealloc_zeldovich(void** arrays){
	fftw_destroy_plan((fftw_plan)arrays[5]);
	fftw_destroy_plan((fftw_plan)arrays[6]);
	fftw_cleanup_threads();
	fftw_free((fftw_complex*)arrays[3]);	
	fftw_free(((double**)arrays[0])[0]);	
	delete[] (double**)arrays[0];
	delete[] (double**)arrays[1];
	delete[] arrays;
}

void** alloc_adhesion(int mesh_num, int par_num, int bin_num, int nt){
	try{
		void** alloc_adh = new void* [8];
		
		/* Allocate long arrays */
		const long mem_double =
			3*par_num*par_num*(par_num+2) + // particle positions
			3*mesh_num*mesh_num*(mesh_num+2) + // velocity field
			mesh_num*mesh_num*(mesh_num+2) + // work space for power spectrum
			mesh_num*mesh_num*(mesh_num+2); // initial expotential
		
		double* alloc_double = (double*) fftw_malloc (sizeof(double)*mem_double);
		int pos = 0;
		
		double** vel_field = new double* [3];
		for (int i = 0; i < 3; i++){
			vel_field[i] = alloc_double + pos;
			pos += mesh_num*mesh_num*(mesh_num+2);
		}
		alloc_adh[0] = (void*) vel_field;
		
		double** par_pos = new double* [3];
		for (int i = 0; i < 3; i++){
			par_pos[i] = alloc_double + pos;
			pos += par_num*par_num*(par_num+2);
		}
		alloc_adh[1] = (void*) par_pos;
		
		double* power_aux = alloc_double + pos;
		pos += mesh_num*mesh_num*(mesh_num+2);
		alloc_adh[2] = (void*) power_aux;
		
		double* expotential_0 = alloc_double + pos;
		pos += mesh_num*mesh_num*(mesh_num+2);
		alloc_adh[3] = (void*) expotential_0;
		
		/* Allocate small arrays */
		fftw_complex* alloc_fftw_c = (fftw_complex*) fftw_malloc (sizeof(fftw_complex)*2*bin_num);
		alloc_adh[4]=(void*)alloc_fftw_c;
		alloc_adh[5]=(void*)(alloc_fftw_c+bin_num);
		
		if (fftw_init_threads() == 0){
			printf("Errors while initializing multi-thread!\n");
			return NULL;
		}
		fftw_plan_with_nthreads(nt); // FFTW_ESTIMATE FFTW_MEASURE
		fftw_plan p_F = fftw_plan_dft_r2c_3d(mesh_num, mesh_num, mesh_num, alloc_double, reinterpret_cast<fftw_complex*>(alloc_double), FFTW_ESTIMATE);
		fftw_plan p_B = fftw_plan_dft_c2r_3d(mesh_num, mesh_num, mesh_num, reinterpret_cast<fftw_complex*>(alloc_double), alloc_double, FFTW_ESTIMATE);
		alloc_adh[6]=(void*)p_F;
		alloc_adh[7]=(void*)p_B;
		
		printf("Allocated %s of memory.\n", humanSize(sizeof(double)*mem_double));
		return alloc_adh;
	}
	catch(bad_alloc&){
		printf("ERROR! Allocation of memory failed!\n");
		return NULL;
	}
	return NULL;
}

void dealloc_adhesion(void** arrays){
	fftw_destroy_plan((fftw_plan)arrays[6]);
	fftw_destroy_plan((fftw_plan)arrays[7]);
	fftw_cleanup_threads();
	fftw_free((fftw_complex*)arrays[4]);
	fftw_free(((double**)arrays[0])[0]);
	delete[] (double**)arrays[0];
	delete[] (double**)arrays[1];
	delete[] arrays;
}

void** alloc_frozen_pot(int mesh_num, int par_num, int bin_num, int nt){
	try{
		void** alloc_fp = new void* [7];
		
		/* Allocate long arrays */
		const long mem_double =
			3*par_num*par_num*(par_num+2) + // particle positions
			3*par_num*par_num*(par_num+2) + // particle velocities
			3*mesh_num*mesh_num*(mesh_num+2) + // force field
			mesh_num*mesh_num*(mesh_num+2); // work space for power spectrum
		
		double* alloc_double = (double*) fftw_malloc (sizeof(double)*mem_double);
		int pos = 0;
		
		double** force_field = new double* [3];
		for (int i = 0; i < 3; i++){
			force_field[i] = alloc_double + pos;
			pos += mesh_num*mesh_num*(mesh_num+2);
		}
		alloc_fp[0] = (void*) force_field;
		
		double*** par_pos = new double** [2];
		par_pos[0] = new double* [3];
		par_pos[1] = new double* [3];
		for (int i = 0; i < 3; i++){
			par_pos[0][i] = alloc_double + pos;
			pos += par_num*par_num*(par_num+2);
		}
		for (int i = 0; i < 3; i++){
			par_pos[1][i] = alloc_double + pos;
			pos += par_num*par_num*(par_num+2);
		}
		alloc_fp[1] = (void*) par_pos;
		
		/* Allocate small arrays */
		double* power_aux = alloc_double + pos;
		pos += mesh_num*mesh_num*(mesh_num+2);
		alloc_fp[2] = (void*) power_aux;
		
		fftw_complex* alloc_fftw_c = (fftw_complex*) fftw_malloc (sizeof(fftw_complex)*2*bin_num);
		alloc_fp[3]=(void*)alloc_fftw_c;
		alloc_fp[4]=(void*)(alloc_fftw_c+bin_num);
		
		if (fftw_init_threads() == 0){
			printf("Errors while initializing multi-thread!\n");
			return NULL;
		}
		fftw_plan_with_nthreads(nt); // FFTW_ESTIMATE FFTW_MEASURE
		fftw_plan p_F = fftw_plan_dft_r2c_3d(mesh_num, mesh_num, mesh_num, alloc_double, reinterpret_cast<fftw_complex*>(alloc_double), FFTW_ESTIMATE);
		fftw_plan p_B = fftw_plan_dft_c2r_3d(mesh_num, mesh_num, mesh_num, reinterpret_cast<fftw_complex*>(alloc_double), alloc_double, FFTW_ESTIMATE);
		alloc_fp[5]=(void*)p_F;
		alloc_fp[6]=(void*)p_B;
		
		printf("Allocated %s of memory.\n", humanSize(sizeof(double)*mem_double));
		return alloc_fp;
	}
	catch(bad_alloc&){
		printf("ERROR! Allocation of memory failed!\n");
		return NULL;
	}
	return NULL;
}

void dealloc_frozen_pot(void** arrays){
	fftw_destroy_plan((fftw_plan)arrays[5]);
	fftw_destroy_plan((fftw_plan)arrays[6]);
	fftw_cleanup_threads();
	fftw_free((fftw_complex*)arrays[3]);	
	fftw_free(((double**)arrays[0])[0]);	
	delete[] (double**)arrays[0];
	delete[] ((double***)arrays[1])[0];
	delete[] ((double***)arrays[1])[1];
	delete[] (double***)arrays[1];
	delete[] arrays;
}

void** alloc_mod_frozen_flow(int mesh_num, int par_num, int bin_num, int nt){
	try{
		void** alloc_zel = new void* [7];
		
		/* Allocate long arrays */
		const long mem_double =
			3*par_num*par_num*(par_num+2) + // particle positions
			par_num*par_num*(par_num+2) + // velocity magnitude
			3*mesh_num*mesh_num*(mesh_num+2) + // velocity field
			mesh_num*mesh_num*(mesh_num+2); // work space for power spectrum
		
		double* alloc_double = (double*) fftw_malloc (sizeof(double)*mem_double);
		int pos = 0;
		
		double** vel_field = new double* [3];
		for (int i = 0; i < 3; i++){
			vel_field[i] = alloc_double + pos;
			pos += mesh_num*mesh_num*(mesh_num+2);
		}
		alloc_zel[0] = (void*) vel_field;
		
		double** par_pos = new double* [4];
		for (int i = 0; i < 4; i++){
			par_pos[i] = alloc_double + pos;
			pos += par_num*par_num*(par_num+2);
		}
		alloc_zel[1] = (void*) par_pos;
		
		/* Allocate small arrays */
		double* power_aux = alloc_double + pos;
		pos += mesh_num*mesh_num*(mesh_num+2);
		alloc_zel[2] = (void*) power_aux;
		
		fftw_complex* alloc_fftw_c = (fftw_complex*) fftw_malloc (sizeof(fftw_complex)*2*bin_num);
		alloc_zel[3]=(void*)alloc_fftw_c;
		alloc_zel[4]=(void*)(alloc_fftw_c+bin_num);
		
		if (fftw_init_threads() == 0){
			printf("Errors while initializing multi-thread!\n");
			return NULL;
		}
		fftw_plan_with_nthreads(nt); // FFTW_ESTIMATE FFTW_MEASURE
		fftw_plan p_F = fftw_plan_dft_r2c_3d(mesh_num, mesh_num, mesh_num, alloc_double, reinterpret_cast<fftw_complex*>(alloc_double), FFTW_ESTIMATE);
		fftw_plan p_B = fftw_plan_dft_c2r_3d(mesh_num, mesh_num, mesh_num, reinterpret_cast<fftw_complex*>(alloc_double), alloc_double, FFTW_ESTIMATE);
		alloc_zel[5]=(void*)p_F;
		alloc_zel[6]=(void*)p_B;
		
		printf("Allocated %s of memory.\n", humanSize(sizeof(double)*mem_double));
		return alloc_zel;
	}
	catch(bad_alloc&){
		printf("ERROR! Allocation of memory failed!\n");
		return NULL;
	}
	return NULL;
}





