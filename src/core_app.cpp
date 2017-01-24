
#include "stdafx.h"
#include "core.h"
#include "core_mesh.h"
#include "core_power.h"
#include <fftw3.h>
#include "CBRNG_Random.h"

#define CORR
#define N_MAX 1

using namespace std;
const double PI = acos(-1.);

static void set_unpert_pos_one_par(Vec_3D<int>& unpert_pos, int par_index, int par_per_dim, int Ng)
{
	unpert_pos.x = (par_index / (par_per_dim * par_per_dim)) * Ng;
	unpert_pos.y = ((par_index / par_per_dim) % par_per_dim) * Ng;
	unpert_pos.z = (par_index % par_per_dim) * Ng;
}

static void set_velocity_one_par(const Vec_3D<int>& unpert_pos, Vec_3D<double>& displ_field, const vector<Mesh> &vel_field)
{
	for (int i = 0; i < 3; i++) displ_field[i] = vel_field[i](unpert_pos);
}

void set_unpert_pos(const Sim_Param &sim, Particle_x* particles)
{
	Vec_3D<int> unpert_pos;
	const int par_per_dim = sim.mesh_num / sim.Ng;
	
	#pragma omp parallel for private(unpert_pos)
	for(int i=0; i< sim.par_num; i++)
	{
		set_unpert_pos_one_par(unpert_pos, i, par_per_dim, sim.Ng);		
		particles[i] = Particle_x(Vec_3D<double>(unpert_pos));
	}
}

void set_unpert_pos_w_vel(const Sim_Param &sim, Particle_v* particles, const vector<Mesh> &vel_field)
{
	Vec_3D<int> unpert_pos;
	Vec_3D<double> velocity;
	const int par_per_dim = sim.mesh_num / sim.Ng;
	
	#pragma omp parallel for private(unpert_pos, velocity)
	for(int i=0; i< sim.par_num; i++)
	{
		set_unpert_pos_one_par(unpert_pos, i, par_per_dim, sim.Ng);
		set_velocity_one_par(unpert_pos, velocity, vel_field);
		particles[i] = Particle_v(Vec_3D<double>(unpert_pos), velocity);
	}
}

void set_pert_pos(const Sim_Param &sim, double db, Particle_x* particles, const vector< Mesh> &vel_field)
{
	Vec_3D<int> unpert_pos;
	Vec_3D<double> displ_field;
	Vec_3D<double> pert_pos;
	
	const int par_per_dim = sim.mesh_num / sim.Ng;
	
	#pragma omp parallel for private(unpert_pos, displ_field, pert_pos)
	for(int i=0; i< sim.par_num; i++)
	{
		set_unpert_pos_one_par(unpert_pos, i, par_per_dim, sim.Ng);
		set_velocity_one_par(unpert_pos, displ_field, vel_field);
		pert_pos = Vec_3D<double>(unpert_pos) + displ_field*db;
		get_per(pert_pos, sim.mesh_num);
		particles[i] = Particle_x(pert_pos);		
	}
}

static void gen_gauss_white_noise(const Sim_Param &sim, Mesh* rho)
{
	// Get keys for each slab in the x axis that this rank contains
	vector<unsigned long> slab_keys;
	slab_keys.resize(rho->N);
	GetSlabKeys(slab_keys.data(), 0, rho->N, sim.seed);
	
	unsigned long ikey, index;
	double rn1, rn2, rn;
		
	#pragma omp parallel for private(ikey, index, rn1, rn2, rn)
	for(long i=0; i<rho->N; ++i) 
	{
		ikey = slab_keys[i];
		for(long j=0; j<rho->N; ++j) 
		{
			for(long k=0; k<rho->N+2; ++k) 
			{
				index = j*rho->N + k;
				GetRandomDoublesWhiteNoise(rn1, rn2, ikey, index);

				rn = rn1*rn1 + rn2*rn2;
				(*rho)(i, j, k) = rn2 * sqrt(-2.0*log(rn)/rn);
			}
		}
    }
	 
	#ifdef CORR
	double t_mean = mean(rho->real(), rho->length);
	double t_std_dev = std_dev(rho->real(), rho->length, t_mean);
	printf("\t[mean = %.12f, stdDev = %.12f]\t-->", t_mean, t_std_dev);
	(*rho)-=t_mean;
	(*rho)/=t_std_dev;
	#endif
	
	double tmp = mean(rho->real(), rho->length);
	printf("\t[mean = %.12f, stdDev = %.12f]\n", tmp, std_dev(rho->real(), rho->length, tmp));
}

static void gen_rho_w_pow_k(const Sim_Param &sim, Mesh* rho)
{
	double k;
	#pragma omp parallel for private(k)
	for(int i=0; i < rho->length / 2;i++)
	{
		k = 2.*PI/sim.box_size*sqrt(get_k_sq(rho->N, i));
		(*rho)[2*i] *= sqrt(lin_pow_spec(sim.power, k));
		(*rho)[2*i+1] *= sqrt(lin_pow_spec(sim.power, k));
	}
}

void gen_rho_dist_k(const Sim_Param &sim, Mesh* rho, const fftw_plan &p_F)
	/**
	Generate density distributions \rho(k) in k-space.
	At first, a gaussian white noise (mean = 0, stdDev = 1) is generated,
	then it is convoluted with given power spectrum.
	**/
{
	printf("Generating gaussian white noise...\n");
	gen_gauss_white_noise(sim, rho);
	
	printf("Generating gaussian white noise in k-sapce...\n");
	fftw_execute_dft_r2c(p_F, *rho);
	
	printf("Generating density distributions with given power spectrum...\n");
	gen_rho_w_pow_k(sim, rho);
}

void pwr_spec_k(const Sim_Param &sim, const Mesh &rho_k, Mesh* power_aux)
{
	/* Computing the power spectrum P(k) */	
	/* Preserve values in rho_k */
	
	double w_k;
	Vec_3D<int> k_vec;
	
	#pragma omp parallel for private(w_k, k_vec)
	for(int i=0; i < rho_k.length/2;i++)
	{
		w_k = 1.;
		get_k_vec(rho_k.N, i, k_vec);
		for (int j = 0; j < 3; j++) if (k_vec[j] != 0) w_k *= pow(sin(PI*k_vec[j]/sim.mesh_num)/(PI*k_vec[j]/sim.mesh_num), sim.order + 1);
		(*power_aux)[2*i+1] = (rho_k[2*i]*rho_k[2*i] + rho_k[2*i+1]*rho_k[2*i+1])/(w_k*w_k);
		(*power_aux)[2*i] = 2.*PI/sim.box_size*k_vec.norm(); // physical k
	}
}

void gen_pow_spec_binned(const Sim_Param &sim, const Mesh &power_aux, vector<fftw_complex>* pwr_spec_binned)
{
	double log_bin = pow(sim.k_max / sim.k_min, 1./sim.bin_num);
	double k;
	int bin;
	printf("Computing binned power spectrum P(k)...\n");
		 
	#pragma omp parallel for
	for (int j = 0; j < sim.bin_num; j++){
		(*pwr_spec_binned)[j][0] = 0.;
		(*pwr_spec_binned)[j][1] = 0.;
	}
	
	for (int i = 0; i < power_aux.length / 2; i++){
		k = power_aux[2*i];
		if ((k <=sim.k_max) && (k>=sim.k_min)){
			bin = (int)(log(k/sim.k_min)/log(log_bin));
			(*pwr_spec_binned)[bin][1] += power_aux[2*i+1]; // P(k)
			(*pwr_spec_binned)[bin][0]++;
		}
	}
		
	const double k_min = sim.k_min*sqrt(log_bin);
	#pragma omp parallel for private(k)
	for (int j = 0; j < sim.bin_num; j++){
		if ((*pwr_spec_binned)[j][0]) (*pwr_spec_binned)[j][1] /= (*pwr_spec_binned)[j][0];
		k = k_min*pow(log_bin, j);
		(*pwr_spec_binned)[j][0] = k;
	}
}

void gen_pot_k(Mesh* rho_k)
{
	printf("Computing potential in k-space...\n");
	double k2;
	#pragma omp parallel for private(k2)
	for(int i=0; i < rho_k->length/2;i++){				
		k2 = get_k_sq(rho_k->N, i);
		if (k2 == 0){
			(*rho_k)[2*i] = 0;
			(*rho_k)[2*i+1] = 0;
		} else{
			(*rho_k)[2*i] /= -(k2*pow(2.*PI/rho_k->N, 2.));
			(*rho_k)[2*i+1] /= -(k2*pow(2.*PI/rho_k->N, 2.));
		}
	}
}

static double CIC_opt(int index, int N)
{
	// optimalization for CIC
	int k_vec[3];
	double k_n[3];
	double U2, U_n, G_n, k2n;
	
	get_k_vec(N, index, k_vec);
	G_n = 0;
	U2 = 1;
	for (int j = 0; j < 3; j++) U2 *= 1./3.*(1+2*cos(PI*k_vec[j]/N));
	for (int n1 = -N_MAX; n1 < N_MAX + 1; n1++)
	{
		k_n[0] = 2 * PI / N*k_vec[0] + 2 * PI*n1;
		for (int n2 = -N_MAX; n2 < N_MAX + 1; n2++)
		{
			k_n[1] = 2 * PI / N*k_vec[1] + 2 * PI*n2;
			for (int n3 = -N_MAX; n3 < N_MAX + 1; n3++)
			{
				k_n[2] = 2 * PI / N*k_vec[2] + 2 * PI*n3;
				U_n = 1.;
				k2n = 0;
				for(int j=0; j<3; j++)
				{
					if (k_n[j] != 0) U_n *= sin(k_n[j] / 2.) / (k_n[j] / 2.);
					k2n += pow(k_n[j],2);
				}
				if (k2n != 0)
				{
					for(int j=0; j<3; j++)
					{										
						G_n += 2 * PI / N*k_vec[j]* // D(k)
						k_n[j]/k2n* // R(k_n)
						pow(U_n, 2.); // W(k) for CIC
					}
				}
			}
		}
	}
	if ((G_n != G_n) || (U2 != U2))
	{
		printf("Gn = %f\tU2 = %f, k = (%i, %i, %i) \n", G_n, U2, k_vec[0], k_vec[1], k_vec[2]);
		return 1.;
	}
	return G_n/U2;
}

void gen_displ_k(vector<Mesh>* vel_field, const Mesh& pot_k)
{
	printf("Computing displacement in k-space...\n");
	double opt;
	int k_vec[3];
	double potential_tmp[2];
	
	#pragma omp parallel for private(opt, k_vec, potential_tmp)
	for(int i=0; i < pot_k.length/2;i++)
	{
		potential_tmp[0] = pot_k[2*i]; // prevent overwriting if vel_field[0] == pot_k
		potential_tmp[1] = pot_k[2*i+1]; // prevent overwriting if vel_field[0] == pot_k
		opt = CIC_opt(i, pot_k.N);
		get_k_vec(pot_k.N, i, k_vec);		
		for(int j=0; j<3;j++)
		{
			(*vel_field)[j][2*i] = k_vec[j]*potential_tmp[1]*(2.*PI/pot_k.N)*opt; // 2*PI/N comes from derivative WITH RESPECT to the mesh coordinates
			(*vel_field)[j][2*i+1] = -k_vec[j]*potential_tmp[0]*(2.*PI/pot_k.N)*opt; // 2*PI/N comes from derivative WITH RESPECT to the mesh coordinates
		}
	}
}

void get_rho_from_par(Particle_x* particles, Mesh* rho, const Sim_Param &sim)
{
	printf("Computing the density field from particle positions...\n");
	double m = pow(sim.Ng, 3);
	
	#pragma omp parallel for
	for (int i = 0; i < rho->length; i++)
	{
		(*rho)[i]=-1.;
	}
	
	#pragma omp parallel for
	for (int i = 0; i < sim.par_num; i++)
	{
		assign_to(rho, particles[i].position, m, sim.order);
	}
}

void gen_dens_binned(const Mesh& rho, vector<int> &dens_binned, const Sim_Param &sim)
{
	printf("Computing binned density field...\n");
	unsigned bin;
	double rho_avg;
	dens_binned.assign(dens_binned.size(), 0);
	
	for (int i = 0; i < sim.mesh_num; i+=sim.Ng)
	{
		for (int j = 0; j < sim.mesh_num; j+=sim.Ng)
		{
			for (int k = 0; k < sim.mesh_num; k+=sim.Ng)
			{
				// Need to go through all mesh cells [i, i+Ng-1]*[j, j+Ng-1], [k, k+Ng, -1]
				rho_avg = 0;
				for (int ii = i; ii < i+sim.Ng; ii++)
				{
					for (int jj = j; jj  < j+sim.Ng; jj++)
					{
						for (int kk = k; kk < k+sim.Ng; kk++)
						{
							rho_avg+=rho(ii, jj, kk);
						}
					}
				}
				rho_avg /= pow(sim.Ng, 3);
				bin = (int)((rho_avg+1)/0.2);
				if (bin >= dens_binned.capacity()) dens_binned.resize(bin+1);
				dens_binned[bin]++;
			}
		}
	}
}