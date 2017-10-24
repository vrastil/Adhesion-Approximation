
#include "stdafx.h"
#include "core.h"
#include "core_mesh.h"
#include "core_power.h"
#include "CBRNG_Random.h"

#define N_MAX 1

using namespace std;

static void set_unpert_pos_one_par(Vec_3D<int>& unpert_pos, const int par_index, const int par_per_dim, const int Ng)
{
	unpert_pos[0] = (par_index / (par_per_dim * par_per_dim)) * Ng;
	unpert_pos[1] = ((par_index / par_per_dim) % par_per_dim) * Ng;
	unpert_pos[2] = (par_index % par_per_dim) * Ng;
}

static void set_velocity_one_par(const Vec_3D<int> unpert_pos, Vec_3D<double>& displ_field, const vector<Mesh> &vel_field)
{
	for (int i = 0; i < 3; i++) displ_field[i] = vel_field[i](unpert_pos);
}

void set_unpert_pos(const Sim_Param &sim, Particle_x* particles)
{
	Vec_3D<int> unpert_pos;
    const int par_per_dim = sim.mesh_num / sim.Ng;
    const int Ng = sim.Ng;
	
	#pragma omp parallel for private(unpert_pos)
	for(unsigned i=0; i< sim.par_num; i++)
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
	for(unsigned i=0; i< sim.par_num; i++)
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
	for(unsigned i=0; i< sim.par_num; i++)
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
	for(unsigned i=0; i< sim.par_num; i++)
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
	for (unsigned i = 0; i < sim.par_num; i++)
	{
		v.fill(0.);
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
	for (unsigned i = 0; i < sim.par_num; i++)
	{
		particles[i].position += particles[i].velocity*(db/2.);
		f_half.fill(0.);
		assign_from(force_field, particles[i].position, &f_half, order);
	//	if (i % (sim.par_num / 13) == 0) printf("particle num = %i, fl = (%f, %f, %f)\n", i, f_half[0], f_half[1], f_half[2]);
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
				 const Vec_3D<double> position, Vec_3D<double>* force)
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
	for (unsigned i = 0; i < sim.par_num; i++)
	{
		particles[i].position += particles[i].velocity*(db/2.);
	}
	
	printf("Creating linked list...\n");
	linked_list->get_linked_list(particles);
//	double fs, fl;
	#pragma omp parallel for private(f_half)
	for (unsigned i = 0; i < sim.par_num; i++)
	{
		f_half.fill(0.);
		
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
	for (unsigned i = 0; i < sim.par_num; i++)
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
	for(unsigned long i=0; i<rho->N1; ++i) 
	{
		ikey = slab_keys[i];
		for(unsigned long j=0; j<rho->N2; ++j) 
		{
			for(unsigned long k=0; k<rho->N3; ++k) 
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
	for(unsigned i=0; i < rho->length / 2;i++)
	{
        k = k0*sqrt(get_k_sq(N, i));
        (*rho)[2*i] *= sqrt(lin_pow_spec(k, sim.power));
        (*rho)[2*i+1] *= sqrt(lin_pow_spec(k, sim.power));

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

template <class T>
void get_rho_from_par(T* particles, Mesh* rho, const Sim_Param &sim)
{
    printf("Computing the density field from particle positions...\n");
    const double m = pow(sim.Ng_pwr, 3.);
    
    const double mesh_mod = (double)sim.mesh_num_pwr/sim.mesh_num;

    rho->assign(-1.);
    
    #pragma omp parallel for
    for (unsigned i = 0; i < sim.par_num; i++)
    {
        assign_to(rho, particles[i].position*mesh_mod, m, sim.order);
    }
}

int get_vel_from_par(Particle_v* particles, vector<Mesh>* vel_field, const Sim_Param &sim)
{
    printf("Computing the velocity field from particle positions...\n");
    const double mesh_mod = (double)sim.mesh_num_pwr/sim.mesh_num;
    const double m = pow(sim.Ng_pwr, 3.);
    for(Mesh& field : *vel_field) field.assign(0.);
    
    #pragma omp parallel for
    for (unsigned i = 0; i < sim.par_num; i++)
    {
        assign_to(vel_field, particles[i].position*mesh_mod, particles[i].velocity*(m*mesh_mod), sim.order);
    }
    return 1;
}

int get_vel_from_par(Particle_x* particles, vector<Mesh>* vel_field, const Sim_Param &sim)
{
    printf("WARNING! Trying to compute velocity divergence with particle positions only! Skipping...\n");
    return 0;
}

void pwr_spec_k(const Sim_Param &sim, const Mesh &rho_k, Mesh* power_aux)
{
    /* Computing the power spectrum P(k)/L^3 -- dimensionLESS!
    in real part of power_aux is stored pk, in imaginary dimensionLESS k
	 Preserve values in rho_k */
	
	double w_k;
    Vec_3D<int> k_vec;
    const int order = sim.order;
    const int NM = sim.mesh_num_pwr;

	#pragma omp parallel for private(w_k, k_vec)
	for(unsigned i=0; i < rho_k.length/2;i++)
	{
		w_k = 1.;
		get_k_vec(NM, i, k_vec);
		for (int j = 0; j < 3; j++) if (k_vec[j] != 0) w_k *= pow(sin(PI*k_vec[j]/NM)/(PI*k_vec[j]/NM), order + 1);
        (*power_aux)[2*i] = (rho_k[2*i]*rho_k[2*i] + rho_k[2*i+1]*rho_k[2*i+1])/(w_k*w_k);
		(*power_aux)[2*i+1] = k_vec.norm();
	}
}

void vel_pwr_spec_k(const Sim_Param &sim, const vector<Mesh> &vel_field, Mesh* power_aux)
{
    /* Computing the velocity power spectrum P(k)/L^3 -- dimensionLESS!
    in real part of power_aux is stored pk, in imaginary dimensionLESS k
	 Preserve values in rho_k */
	
	double w_k;
    Vec_3D<int> k_vec;
    const int order = sim.order;
    const int NM = sim.mesh_num_pwr;

    double vel_div_re, vel_div_im, k;

	#pragma omp parallel for private(w_k, k_vec)
	for(unsigned i=0; i < power_aux->length/2;i++)
	{
        w_k = 1.;
        vel_div_re = vel_div_im = 0;
		get_k_vec(NM, i, k_vec);
        for (int j = 0; j < 3; j++){
            k = 2*PI*k_vec[j] / NM;
            if (k != 0) w_k *= pow(sin(k/2.)/(k/2.), order + 1);
            vel_div_re += vel_field[j][2*i]*k; // do not care about Re <-> Im in 2*PI*i/N, norm only
            vel_div_im += vel_field[j][2*i+1]*k;
        } 
        
        (*power_aux)[2*i] = (vel_div_re*vel_div_re + vel_div_im*vel_div_im)/(w_k*w_k);
		(*power_aux)[2*i+1] = k_vec.norm();
	}
}

void gen_cqty_binned(const double x_min, const double x_max,
                    const Mesh &qty_mesh, Data_x_y<double>& qty_binned, const double mod_q, const double mod_x)
{
    /* bin some complex quantity on mesh in logarithmic bins, assuming:
       Q(x) = mod_q*qty_mesh[2*i]
       x = mod_x*qty_mesh[2*i+1]
       [x] = [x_min] = [x_max]
    */
    const double log_bin = pow(x_max/x_min, 1./qty_binned.size());

    qty_binned.fill(0);

    double x;
    int bin;
    #pragma omp parallel for private(x, bin)
    for (unsigned i = 0; i < qty_mesh.length / 2; i++){
        x = qty_mesh[2*i+1];
        if ((x <x_max) && (x>=x_min)){
            bin = (int)(log(x/x_min)/log(log_bin));
            #pragma omp atomic
            qty_binned.y[bin] += qty_mesh[2*i];
            #pragma omp atomic
            qty_binned.x[bin]++;
        }
    }
    const double x_min_ = mod_x*x_min*sqrt(log_bin);
    unsigned i = 0;
    for (unsigned j = 0; j < qty_binned.size(); ){
        if (qty_binned.x[j]){
            qty_binned.y[j] *= mod_q / qty_binned.x[j];
            x = x_min_*pow(log_bin, i);
            qty_binned.x[j] = x;
            j++;
        }else{
            qty_binned.erase(j);
        }
        i++;
	}
}

void gen_rqty_binned(const double x_min, const double x_max, const double x_0,
    const Mesh &qty_mesh, Data_x_y<double>& qty_binned, const double mod_q)
{
/* bin some real quantity on mesh in linear bins
   x_min and x_max are in [Mpc/h], transformation from mesh distances to Mpc/h assumed to be x_0
   Q(x) = qty_mesh[i]
   x = ([i,j,k].norm())*x_0
   
*/
    const double lin_bin = (x_max - x_min) / qty_binned.size();

    qty_binned.fill(0);

    double x;
    int bin;
    #pragma omp parallel for private(x, bin)
    for (unsigned i = 0; i < qty_mesh.N; i++){
        for (unsigned j = 0; i < qty_mesh.N; i++){
            for (unsigned k = 0; i < qty_mesh.N; i++){

                x = x_0*sqrt(i*i+j*j+k*k);
                if ((x <x_max) && (x>=x_min)){
                    bin = (int)((x-x_min)/lin_bin);
                    #pragma omp atomic
                    qty_binned.y[bin] += qty_mesh(i,j,k);
                    #pragma omp atomic
                    qty_binned.x[bin]++;
                }
            }
        }
    }
    const double x_min_ = x_min+lin_bin/2;
    unsigned i = 0;
    for (unsigned j = 0; j < qty_binned.size(); ){
        if (qty_binned.x[j]){
            qty_binned.y[j] *= mod_q / qty_binned.x[j];
            x = x_min_ + lin_bin*i;
            qty_binned.x[j] = x;
            j++;
        }else{
            qty_binned.erase(j);
        }
        i++;
    }
}

void gen_pow_spec_binned(const Sim_Param &sim, const Mesh &power_aux, Data_x_y<double>* pwr_spec_binned)
{
    const int L = sim.box_size;
    #ifndef OLD_NORM
    const double mod_pk = pow(L, 3.); // P(k) - dimensionFULL!
    #else
    const double mod_pk = 1.;
    #endif
    const double mod_k = 2.*PI/L;
    printf("Computing binned power spectrum...\n");
	gen_cqty_binned(1, sim.mesh_num_pwr, power_aux, *pwr_spec_binned, mod_pk, mod_k);
}

void gen_pow_spec_binned_from_extrap(const Sim_Param &sim, const Extrap_Pk &P_k, Data_x_y<double>* pwr_spec_binned)
{
    const double k_max = sim.k_par.k_print.upper;
    const double k_min = sim.k_par.k_print.lower;
    const double log_bin = pow(k_max/k_min, 1./pwr_spec_binned->size());
    double k;

    #pragma omp parallel for private(k)
	for (unsigned j = 0; j < pwr_spec_binned->size(); j++){
		k = k_min*pow(log_bin, j);
		pwr_spec_binned->x[j] = k;
        pwr_spec_binned->y[j] = P_k.eval(k);
    }
}

void gen_corr_func_binned(const Sim_Param &sim, const Mesh &power_aux, Data_x_y<double>* corr_func_binned)
{
    printf("Computing correlation function via FFT of power spectrum on a mesh...\n");
	gen_rqty_binned(1, 200, sim.x_0_pwr(), power_aux, *corr_func_binned, 1);
}

template<class T>
void gen_corr_func_binned_pp(const Sim_Param &sim, T* particles, Data_x_y<double>* corr_func_binned)
{ /* x_min and x_max are in [Mpc/h], transformation from particles distances to Mpc/h assumed to be sim.x_0() */
    const double x_min = sim.x_corr.lower;
    const double x_max = sim.x_corr.upper;
    const double x_0 = sim.x_0();
    printf("Computing correlation function via direct particle-particle sum...\n");

    const double lin_bin = (x_max - x_min) / corr_func_binned->size(); // [lin_bin] = Mpc/h
    const unsigned Np = sim.par_num;
    const unsigned Np_max = Np;
    const unsigned N = sim.mesh_num;
    double x, dV;
    int bin;
    const double V = pow(sim.box_size/2, 3.);
    const int par_factor = (2*Np - Np_max - 1)*Np_max / 2; // number of particles over which the previus sum was done
    const double x_min_ = x_min+lin_bin/2;

    corr_func_binned->fill(0);

    #pragma omp parallel for private(x)
    for (unsigned i = 0; i < Np_max; i++){
        for (unsigned j = i + 1; j < Np; j++)
        {
            x = x_0*get_distance(particles[i].position, particles[j].position, N);
            if ((x <x_max) && (x>=x_min)){
                bin = (int)((x-x_min)/lin_bin);
                #pragma omp atomic
                corr_func_binned->y[bin]++;
            }
        }
    }

    #pragma omp parallel for private(x, dV)
    for (unsigned j = 0; j < corr_func_binned->size(); j++){
        x = x_min_ + lin_bin*j;
        corr_func_binned->x[j] = x;
        dV = 4*PI*x*x*lin_bin + PI/3*pow(lin_bin, 3);
        corr_func_binned->y[j] /= par_factor*dV / V;
    //    corr_func_binned->y[j]--; // definition factor
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
    const unsigned l_half = rho_k.length/2;

	#pragma omp parallel for private(k2)
	for(unsigned i=0; i < l_half;i++){				
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

static double CIC_opt(Vec_3D<double> k_vec, const double a)
{
	double k_n[3];
	double U2, U_n, G_n, k2n;
	
	G_n = 0;
	U2 = 1;
	for (int j = 0; j < 3; j++) U2 *= 1./3.*(1+2*cos(k_vec[j]/2.));
	for (int n1 = -N_MAX; n1 < N_MAX + 1; n1++)
	{
		k_n[0] = k_vec[0] + 2 * PI*n1;
		for (int n2 = -N_MAX; n2 < N_MAX + 1; n2++)
		{
			k_n[1] = k_vec[1] + 2 * PI*n2;
			for (int n3 = -N_MAX; n3 < N_MAX + 1; n3++)
			{
				k_n[2] = k_vec[2] + 2 * PI*n3;
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
						G_n += k_vec[j]* // D(k)
						k_n[j]/k2n*pow(S2_shape(k2n, a), 2)* // R(k_n)
						pow(U_n, 2.); // W(k) for CIC
					}
				}
			}
		}
	}
	if ((G_n != G_n) || (U2 != U2))
	{
		printf("Gn = %f\tU2 = %f, k = (%f, %f, %f) \n", G_n, U2, k_vec[0], k_vec[1], k_vec[2]);
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
    Vec_3D<int> k_vec;
    Vec_3D<double> k_vec_phys;
    double potential_tmp[2];
    
    const int N = (*vel_field)[0].N; // for case when pot_k is different mesh than vel_field
    const unsigned l_half = (*vel_field)[0].length/2;
	
	#pragma omp parallel for private(opt, k_vec, k_vec_phys, potential_tmp)
	for(unsigned i=0; i < l_half;i++)
	{
		potential_tmp[0] = pot_k[2*i]; // prevent overwriting if vel_field[0] == pot_k
        potential_tmp[1] = pot_k[2*i+1]; // prevent overwriting if vel_field[0] == pot_k
        get_k_vec(N, i, k_vec);
        k_vec_phys = Vec_3D<double>(k_vec)*(2.*PI/N);	
        // no optimalization
        if (a == -1) opt = 1.;
        // optimalization for CIC and S2 shaped particle
		else opt = CIC_opt(k_vec_phys, a);
		for(unsigned j=0; j<3;j++)
		{
			// 2*PI/N comes from derivative WITH RESPECT to the mesh coordinates
			(*vel_field)[j][2*i] = k_vec_phys[j]*potential_tmp[1]*opt;
			(*vel_field)[j][2*i+1] = -k_vec_phys[j]*potential_tmp[0]*opt;
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
    const unsigned N = rho.N;

	dens_binned.assign(dens_binned.size(), 0);
    
    #pragma omp parallel for private(bin, rho_avg)
	for (unsigned i = 0; i < N; i+=Ng_pwr)
	{
		for (unsigned j = 0; j < N; j+=Ng_pwr)
		{
			for (unsigned k = 0; k < N; k+=Ng_pwr)
			{
				// Need to go through all mesh cells [i, i+Ng-1]*[j, j+Ng-1], [k, k+Ng, -1]
				rho_avg = 0;
				for (unsigned ii = i; ii < i+Ng_pwr; ii++)
				{
					for (unsigned jj = j; jj  < j+Ng_pwr; jj++)
					{
						for (unsigned kk = k; kk < k+Ng_pwr; kk++)
						{
							rho_avg+=rho(ii, jj, kk);
						}
					}
				}
				rho_avg /= pow(Ng_pwr, 3);
                bin = (int)((rho_avg+1)/0.1);
                if (bin >= dens_binned.size()) bin = dens_binned.size() - 1;
                // if (bin >= dens_binned.capacity()) dens_binned.resize(bin+1);
                #pragma omp atomic
                dens_binned[bin]++;
                //dens_binned[bin] += pow(sim.Ng_pwr, 3);
			}
		}
	}
}

template void get_rho_from_par(Particle_x*, Mesh*, const Sim_Param&);
template void get_rho_from_par(Particle_v*, Mesh*, const Sim_Param&);
template void gen_corr_func_binned_pp(const Sim_Param&, Particle_x*, Data_x_y<double>*);
template void gen_corr_func_binned_pp(const Sim_Param&, Particle_v*, Data_x_y<double>*);