
#include "stdafx.h"
#include "core.h"
#include "core_mesh.h"
#include "core_power.h"
#include <fftw3.h>
#include "CBRNG_Random.h"

#define N_MAX 1

using namespace std;
const double PI = acos(-1.);

static void set_unpert_pos_one_par(Vec_3D<int>& unpert_pos, const int par_index, const int par_per_dim, const int Ng)
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
    const int Ng = sim.Ng;
	
	#pragma omp parallel for private(unpert_pos)
	for(int i=0; i< sim.par_num; i++)
	{
		set_unpert_pos_one_par(unpert_pos, i, par_per_dim, Ng);		
		particles[i] = Particle_x(Vec_3D<double>(unpert_pos));
	}
}

void set_unpert_pos_w_vel(const Sim_Param &sim, Particle_v* particles, const vector<Mesh> &vel_field)
{
	Vec_3D<int> unpert_pos;
	Vec_3D<double> velocity;
	const int par_per_dim = sim.mesh_num / sim.Ng;
    const int Ng = sim.Ng;
    
	#pragma omp parallel for private(unpert_pos, velocity)
	for(int i=0; i< sim.par_num; i++)
	{
		set_unpert_pos_one_par(unpert_pos, i, par_per_dim, Ng);
		set_velocity_one_par(unpert_pos, velocity, vel_field);
		particles[i] = Particle_v(Vec_3D<double>(unpert_pos), velocity);
	}
}

void set_pert_pos(const Sim_Param &sim, const double db, Particle_x* particles, const vector< Mesh> &vel_field)
{
	Vec_3D<int> unpert_pos;
	Vec_3D<double> displ_field;
	Vec_3D<double> pert_pos;
	
    const int par_per_dim = sim.mesh_num / sim.Ng;
    const int Ng = sim.Ng;
    const int Nm = sim.mesh_num;
	
	#pragma omp parallel for private(unpert_pos, displ_field, pert_pos)
	for(int i=0; i< sim.par_num; i++)
	{
		set_unpert_pos_one_par(unpert_pos, i, par_per_dim, Ng);
		set_velocity_one_par(unpert_pos, displ_field, vel_field);
		pert_pos = Vec_3D<double>(unpert_pos) + displ_field*db;
		get_per(pert_pos, Nm);
		particles[i] = Particle_x(pert_pos);		
	}
}

void set_pert_pos_w_vel(const Sim_Param &sim, const double db, Particle_v* particles, const vector< Mesh> &vel_field)
{
	Vec_3D<int> unpert_pos;
	Vec_3D<double> velocity;
	Vec_3D<double> pert_pos;
	
	const int par_per_dim = sim.mesh_num / sim.Ng;
	const int Ng = sim.Ng;
    const int Nm = sim.mesh_num;

	#pragma omp parallel for private(unpert_pos, velocity, pert_pos)
	for(int i=0; i< sim.par_num; i++)
	{
		set_unpert_pos_one_par(unpert_pos, i, par_per_dim, Ng);
		set_velocity_one_par(unpert_pos, velocity, vel_field);
		pert_pos = Vec_3D<double>(unpert_pos) + velocity*db;
		get_per(pert_pos, Nm);
		particles[i] = Particle_v(pert_pos, velocity);		
	}
}

void upd_pos_first_order(const Sim_Param &sim, const double db, Particle_x* particles, const vector< Mesh> &vel_field)
{
	// Euler method
	
	Vec_3D<double> v;
	const int order = sim.order;
    const int Nm = sim.mesh_num;

	#pragma omp parallel for private(v)
	for (int i = 0; i < sim.par_num; i++)
	{
		v.assign(0., 0., 0.);
		assign_from(vel_field, particles[i].position, &v, order);
		particles[i].position += v*db;
		get_per(particles[i].position, Nm);
	}
}

void upd_pos_second_order(const Sim_Param &sim, const double db, const double b, Particle_v* particles, const vector< Mesh> &force_field)
{
	// Leapfrog method for frozen-flow
	
    Vec_3D<double> f_half;
    const int order = sim.order;
    const int Nm = sim.mesh_num;

	#pragma omp parallel for private(f_half)
	for (int i = 0; i < sim.par_num; i++)
	{
		particles[i].position += particles[i].velocity*(db/2.);
		f_half.assign(0., 0., 0.);
		assign_from(force_field, particles[i].position, &f_half, order);
	//	if (i % (sim.par_num / 13) == 0) printf("particle num = %i, fl = (%f, %f, %f)\n", i, f_half.x, f_half.y, f_half.z);
		f_half = (particles[i].velocity - f_half)*(-3/(2.*(b-db/2.))); // <- EOM
		
		particles[i].velocity += f_half*db;
		particles[i].position += particles[i].velocity*(db/2.);
		get_per(particles[i].position, Nm);
	}
}

static double force_ref(const double r, const double a){
	// Reference force for an S_2-shaped particle
	double z = 2 * r / a;
	if (z > 2) return 1 / (r*r);
	else if (z > 1) return (12 / (z*z) - 224 + 896 * z - 840 * z*z + 224 * pow(z, 3) +
							70 * pow(z, 4) - 48 * pow(z, 5) + 7 * pow(z, 6)) / (35 * a*a);
	else return (224 * z - 224 * pow(z, 3) + 70 * pow(z, 4) + 48 * pow(z, 5) - 21 * pow(z, 7)) / (35 * a*a);
}

static double force_tot(double r){
	return 1 / (r*r);
}

void force_short(const Sim_Param &sim, const LinkedList& linked_list, Particle_v *particles,
				 const Vec_3D<double> &position, Vec_3D<double>* force)
{	// Calculate short range force in position, force is added
	int p;
    Vec_3D<int> y;
	Vec_3D<double> dr_vec;
	double dr;
    const double m = pow(sim.Ng, 3);
    const int Nm = sim.mesh_num;
    const int a = sim.a;

	const Vec_3D<int> z =(Vec_3D<int>)(position / sim.Hc); // chain position of particle
	for (y[0] = z[0] -1; y[0] < z[0] + 2; y[0]++)
	{
		for (y[1] = z[1] - 1; y[1] < z[1] + 2; y[1]++){
		
			for (y[2] = z[2] - 1; y[2] < z[2] + 2; y[2]++)
			{
				p = linked_list.HOC(y);
				while (p != -1){
					dr_vec = get_sgn_distance(particles[p].position, position, Nm);
					dr = dr_vec.norm();
			//		dr = get_distance(particles[p].position, position, sim.mesh_num);
					if ((dr < sim.rs) && (dr != 0)) // Short range force is set 0 for separation larger than cutoff radius
			//		if (dr != 0) // Short range force is set 0 for separation larger than cutoff radius
					{ 
						(*force) += (m*(force_tot(dr) - force_ref(dr, a))/(dr*4*PI))*dr_vec;
					}
					p = linked_list.LL[p];
				}
			}
		}
	}
}

void upd_pos_second_order_w_short_force(const Sim_Param &sim, LinkedList* linked_list, const double db,
										const double b, Particle_v* particles, const vector< Mesh> &force_field)
{
	// Leapfrog method for modified frozen-flow
	Vec_3D<double> f_half;
    const int order = sim.order;
    const int Nm = sim.mesh_num;

	// three loops to compute distances and velocities at fixed positon
	
	#pragma omp parallel for private(f_half)
	for (int i = 0; i < sim.par_num; i++)
	{
		particles[i].position += particles[i].velocity*(db/2.);
	}
	
	printf("Creating linked list...\n");
	linked_list->get_linked_list(particles);
//	double fs, fl;
	#pragma omp parallel for private(f_half)
	for (int i = 0; i < sim.par_num; i++)
	{
		f_half.assign(0., 0., 0.);
		
		// long-range force
		assign_from(force_field, particles[i].position, &f_half, order);
//		fl = f_half.norm();
//		fs = fl;
		// short range force
		force_short(sim, *linked_list, particles, particles[i].position, &f_half);
//		fs -= f_half.norm();
//		if (i % (sim.par_num / 13) == 0) printf("particle num = %i, fl = %f, fs = %f\n",i, fl, fs);
		f_half = (particles[i].velocity - f_half)*(-3/(2.*(b-db/2.))); // <- EOM
		
		particles[i].velocity += f_half*db;
	}
	
	#pragma omp parallel for private(f_half)
	for (int i = 0; i < sim.par_num; i++)
	{
		particles[i].position += particles[i].velocity*(db/2.);
		get_per(particles[i].position, Nm);
	}
}

static void gen_gauss_white_noise(const Sim_Param &sim, Mesh* rho)
{
	// Get keys for each slab in the x axis that this rank contains
	vector<unsigned long> slab_keys;
	slab_keys.resize(rho->N1);
	GetSlabKeys(slab_keys.data(), 0, rho->N1, sim.seed);
	
	unsigned long ikey, index;
	double rn1, rn2, rn;
		
	#pragma omp parallel for private(ikey, index, rn1, rn2, rn)
	for(long i=0; i<rho->N1; ++i) 
	{
		ikey = slab_keys[i];
		for(long j=0; j<rho->N2; ++j) 
		{
			for(long k=0; k<rho->N3; ++k) 
			{
                #ifdef REAL_NOISE
                index = j*rho->N2 + k; // N2 to have the same code as HACC
                #else
                index = j*rho->N3 + k; // N3 to have each number unique
                #endif
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
    printf("\t[min = %.12f, max = %.12f]\n", min(rho->real(), rho->length), max(rho->real(), rho->length));
}

static void gen_rho_w_pow_k(const Sim_Param &sim, Mesh* rho)
{
    double k;
    const double L = sim.box_size;
    const double k0 = 2.*PI/L;
    const int N = rho->N;
	#pragma omp parallel for private(k)
	for(int i=0; i < rho->length / 2;i++)
	{
        k = k0*sqrt(get_k_sq(N, i));
        (*rho)[2*i] *= sqrt(lin_pow_spec(sim.power, k));
        (*rho)[2*i+1] *= sqrt(lin_pow_spec(sim.power, k));

        #ifndef OLD_NORM
        (*rho)[2*i] /= pow(L, 3/2.);
        (*rho)[2*i+1] /= pow(L, 3/2.);
        #endif

        #ifndef REAL_NOISE
        (*rho)[2*i] /= sqrt(2.);
        (*rho)[2*i+1] /= sqrt(2.);
            #ifdef FFTW_SYM
            (*rho)[2*i] *= pow(N, 3/2.);
            (*rho)[2*i+1] *= pow(N, 3/2.);
            #endif
        #else
            #ifndef FFTW_SYM
            (*rho)[2*i] *= pow(N, 3/2.);
            (*rho)[2*i+1] *= pow(N, 3/2.);
            #endif
        #endif
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
    
    #ifdef REAL_NOISE
	printf("Generating gaussian white noise in k-sapce...\n");
    fftw_execute_dft_r2c(p_F, *rho);

    double t_mean = mean(rho->real(), rho->length);
	double t_std_dev = std_dev(rho->real(), rho->length, t_mean);
    printf("\t[mean = %.12f, stdDev = %.12f]\n", t_mean, t_std_dev);
    printf("\t[min = %.12f, max = %.12f]\n", min(rho->real(), rho->length), max(rho->real(), rho->length));
    #endif

	printf("Generating density distributions with given power spectrum...\n");
	gen_rho_w_pow_k(sim, rho);
}

void pwr_spec_k(const Sim_Param &sim, const Mesh &rho_k, Mesh* power_aux)
{
	/* Computing the power spectrum P(k) */	
	/* Preserve values in rho_k */
	
	double w_k;
    Vec_3D<int> k_vec;
    const int order = sim.order;
    const int Nm = sim.mesh_num;
    const int NM = sim.mesh_num_pwr;
    const int L = sim.box_size;
    const double k0 = 2.*PI/L;

	#pragma omp parallel for private(w_k, k_vec)
	for(int i=0; i < rho_k.length/2;i++)
	{
		w_k = 1.;
		get_k_vec(rho_k.N, i, k_vec);
		for (int j = 0; j < 3; j++) if (k_vec[j] != 0) w_k *= pow(sin(PI*k_vec[j]/NM)/(PI*k_vec[j]/NM), order + 1);
        (*power_aux)[2*i+1] = (rho_k[2*i]*rho_k[2*i] + rho_k[2*i+1]*rho_k[2*i+1])/(w_k*w_k);
        #ifndef OLD_NORM
        (*power_aux)[2*i+1] *= pow(L, 3.);
        #endif

		(*power_aux)[2*i] = k0*k_vec.norm(); // physical k
	}
}

void gen_pow_spec_binned(const Sim_Param &sim, const Mesh &power_aux, vector<double_2>* pwr_spec_binned)
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
    
    #pragma omp parallel for private(k, bin)
	for (int i = 0; i < power_aux.length / 2; i++){
		k = power_aux[2*i];
		if ((k <=sim.k_max) && (k>=sim.k_min)){
            bin = (int)(log(k/sim.k_min)/log(log_bin));
            #pragma omp atomic
            (*pwr_spec_binned)[bin][1] += power_aux[2*i+1]; // P(k)
            #pragma omp atomic
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

void gen_pot_k(const Mesh& rho_k, Mesh* pot_k)
{   /*
    pot_k can be Mesh of differen (bigger) size rho_k,
    !!!> ALL physical FACTORS ARE therefore TAKEN FROM rho_k <!!!
    */
	printf("Computing potential in k-space...\n");
    double k2;
    const int N = rho_k.N; // for case when pot_k is different mesh than vel_field
    const int l_half = rho_k.length/2;

	#pragma omp parallel for private(k2)
	for(int i=0; i < l_half;i++){				
		k2 = get_k_sq(N, i);
		if (k2 == 0){
			(*pot_k)[2*i] = 0;
			(*pot_k)[2*i+1] = 0;
		} else{
			(*pot_k)[2*i] = -rho_k[2*i]/(k2*pow(2.*PI/N, 2.));
			(*pot_k)[2*i+1] = -rho_k[2*i+1]/(k2*pow(2.*PI/N, 2.));
		}
	}
}

void gen_pot_k(Mesh* rho_k){ gen_pot_k(*rho_k, rho_k); }

static double S2_shape(const double k2, const double a)
{
	if (a == 0) return 1.;
	
	double t = sqrt(k2)*a / 2;
	return 12 / pow(t, 4)*(2 - 2 * cos(t) - t*sin(t));
}

static double CIC_opt(const int index, const int N, const double a)
{
	// optimalization for CIC and S2 shaped particle
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
						k_n[j]/k2n*pow(S2_shape(k2n, a), 2)* // R(k_n)
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

void gen_displ_k_S2(vector<Mesh>* vel_field, const Mesh& pot_k, const double a)
{   /*
    pot_k can be Mesh of differen (bigger) size than each vel_field,
    !!!> ALL physical FACTORS ARE therefore TAKEN FROM vel_field[0] <!!!
    */
	if (a == -1) printf("Computing displacement in k-space...\n");
	else if (a == 0) printf("Computing displacement in k-space with CIC opt...\n");
	else printf("Computing force in k-space for S2 shaped particles with CIC opt...\n");

	double opt;
	int k_vec[3];
    double potential_tmp[2];
    
    const int N = (*vel_field)[0].N; // for case when pot_k is different mesh than vel_field
    const int l_half = (*vel_field)[0].length/2;
	
	#pragma omp parallel for private(opt, k_vec, potential_tmp)
	for(int i=0; i < l_half;i++)
	{
		potential_tmp[0] = pot_k[2*i]; // prevent overwriting if vel_field[0] == pot_k
		potential_tmp[1] = pot_k[2*i+1]; // prevent overwriting if vel_field[0] == pot_k
		if (a == -1) opt = 1.;
		else opt = CIC_opt(i, N, a);
		get_k_vec(N, i, k_vec);		
		for(int j=0; j<3;j++)
		{
			// 2*PI/N comes from derivative WITH RESPECT to the mesh coordinates
			(*vel_field)[j][2*i] = k_vec[j]*potential_tmp[1]*(2.*PI/N)*opt;
			(*vel_field)[j][2*i+1] = -k_vec[j]*potential_tmp[0]*(2.*PI/N)*opt;
		}
	}
}

void gen_displ_k(vector<Mesh>* vel_field, const Mesh& pot_k) {gen_displ_k_S2(vel_field, pot_k, -1);}

void gen_displ_k_cic(vector<Mesh>* vel_field, const Mesh& pot_k) {gen_displ_k_S2(vel_field, pot_k, 0.);}

void gen_dens_binned(const Mesh& rho, vector<int> &dens_binned, const Sim_Param &sim)
{
	printf("Computing binned density field...\n");
	unsigned bin;
    double rho_avg;
    const int Ng_pwr = sim.Ng_pwr;

	dens_binned.assign(dens_binned.size(), 0);
    
    #pragma omp parallel for private(bin, rho_avg)
	for (int i = 0; i < rho.N; i+=Ng_pwr)
	{
		for (int j = 0; j < rho.N; j+=Ng_pwr)
		{
			for (int k = 0; k < rho.N; k+=Ng_pwr)
			{
				// Need to go through all mesh cells [i, i+Ng-1]*[j, j+Ng-1], [k, k+Ng, -1]
				rho_avg = 0;
				for (int ii = i; ii < i+Ng_pwr; ii++)
				{
					for (int jj = j; jj  < j+Ng_pwr; jj++)
					{
						for (int kk = k; kk < k+Ng_pwr; kk++)
						{
							rho_avg+=rho(ii, jj, kk);
						}
					}
				}
				rho_avg /= pow(Ng_pwr, 3);
                bin = (int)((rho_avg+1)/0.2);
                if (bin >= dens_binned.size()) bin = dens_binned.size() - 1;
                // if (bin >= dens_binned.capacity()) dens_binned.resize(bin+1);
                #pragma omp atomic
                dens_binned[bin]++;
                //dens_binned[bin] += pow(sim.Ng_pwr, 3);
			}
		}
	}
}