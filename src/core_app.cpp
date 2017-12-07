
#include "stdafx.h"
#include "core.h"
#include "core_mesh.h"
#include "core_power.h"
#include "CBRNG_Random.h"

#ifndef ORDER
#define ORDER 1
#endif

using namespace std;

static void set_unpert_pos_one_par(Vec_3D<int>& unpert_pos, const unsigned par_index, const unsigned par_per_dim, const unsigned Ng)
{
	unpert_pos[0] = (par_index / (par_per_dim * par_per_dim)) * Ng;
	unpert_pos[1] = ((par_index / par_per_dim) % par_per_dim) * Ng;
	unpert_pos[2] = (par_index % par_per_dim) * Ng;
}

static void set_velocity_one_par(const Vec_3D<int> unpert_pos, Vec_3D<double>& displ_field, const vector<Mesh> &vel_field)
{
	for (unsigned i = 0; i < 3; i++) displ_field[i] = vel_field[i](unpert_pos);
}

void set_unpert_pos(const Sim_Param &sim, Particle_x* particles)
{
	Vec_3D<int> unpert_pos;
    const unsigned par_per_dim = sim.box_opt.par_num_1d;
    const unsigned Ng = sim.box_opt.Ng;
    const unsigned Np = sim.box_opt.par_num;
	
	#pragma omp parallel for private(unpert_pos)
	for(unsigned i=0; i< Np; i++)
	{
		set_unpert_pos_one_par(unpert_pos, i, par_per_dim, Ng);		
		particles[i] = Particle_x(Vec_3D<double>(unpert_pos));
	}
}

void set_unpert_pos_w_vel(const Sim_Param &sim, Particle_v* particles, const vector<Mesh> &vel_field)
{
	Vec_3D<int> unpert_pos;
	Vec_3D<double> velocity;
	const unsigned par_per_dim = sim.box_opt.par_num_1d;
    const unsigned Ng = sim.box_opt.Ng;
    
    const unsigned Np = sim.box_opt.par_num;
	#pragma omp parallel for private(unpert_pos, velocity)
	for(unsigned i=0; i< Np; i++)
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
	
    const unsigned par_per_dim = sim.box_opt.par_num_1d;
    const unsigned Ng = sim.box_opt.Ng;
    const unsigned Nm = sim.box_opt.mesh_num;
    const unsigned Np = sim.box_opt.par_num;
	
	#pragma omp parallel for private(unpert_pos, displ_field, pert_pos)
	for(unsigned i=0; i< Np; i++)
	{
		set_unpert_pos_one_par(unpert_pos, i, par_per_dim, Ng);
		set_velocity_one_par(unpert_pos, displ_field, vel_field);
		pert_pos = Vec_3D<double>(unpert_pos) + displ_field*db;
		get_per(pert_pos, Nm);
		particles[i] = Particle_x(pert_pos);		
	}
}

void set_pert_pos_w_vel(const Sim_Param &sim, const double a, Particle_v* particles, const vector< Mesh> &vel_field)
{
	Vec_3D<int> unpert_pos;
	Vec_3D<double> velocity;
	Vec_3D<double> pert_pos;
	
	const unsigned par_per_dim = sim.box_opt.par_num_1d;
	const unsigned Ng = sim.box_opt.Ng;
    const unsigned Nm = sim.box_opt.mesh_num;
    const unsigned Np = sim.box_opt.par_num;

    const double D = growth_factor(a, sim.cosmo); // growth factor
    const double dDda = growth_change(a, sim.cosmo); // dD / da

	#pragma omp parallel for private(unpert_pos, velocity, pert_pos)
	for(unsigned i=0; i< Np; i++)
	{
		set_unpert_pos_one_par(unpert_pos, i, par_per_dim, Ng);
		set_velocity_one_par(unpert_pos, velocity, vel_field);
		pert_pos = Vec_3D<double>(unpert_pos) + velocity*D;
		get_per(pert_pos, Nm);
		particles[i] = Particle_v(pert_pos, velocity*dDda);		
	}
}

void stream_step(const Sim_Param &sim, const double da, Particle_v* particles)
{
    const unsigned Np = sim.box_opt.par_num;
    #pragma omp parallel for
	for (unsigned i = 0; i < Np; i++)
	{
        particles[i].position += particles[i].velocity*da;
    }
}

void kick_step_no_momentum(const Sim_Param &sim, const double a, Particle_v* particles, const vector< Mesh> &vel_field)
{
    // no memory of previus velocity, 1st order ODE
    const unsigned Np = sim.box_opt.par_num;
    Vec_3D<double> vel;
    const double dDda = growth_change(a, sim.cosmo); // dD / da
    
    #pragma omp parallel for private(vel)
    for (unsigned i = 0; i < Np; i++)
	{
        vel.fill(0.);
        assign_from(vel_field, particles[i].position, &vel);
        particles[i].velocity = vel*dDda;
    }
}

void kick_step_w_momentum(const Sim_Param &sim, const double a, const double da, Particle_v* particles, const vector< Mesh> &force_field)
{
    // classical 2nd order ODE
    const unsigned Np = sim.box_opt.par_num;
    Vec_3D<double> force;
    const double D = growth_factor(a, sim.cosmo);
    const double OL = sim.cosmo.Omega_L()*pow(a,3);
    const double Om = sim.cosmo.Omega_m;
    // -3/2a represents usual EOM, the rest are LCDM corrections
    const double f1 = 3/(2.*a)*(Om+2*OL)/(Om+OL);
    const double f2 = 3/(2.*a)*Om/(Om+OL)*D/a;
    
    #pragma omp parallel for private(force)
    for (unsigned i = 0; i < Np; i++)
	{
        force.fill(0.);
        assign_from(force_field, particles[i].position, &force);
        force = force*f2 - particles[i].velocity*f1;		
        particles[i].velocity += force*da;
    }
}

double force_ref(const double r, const double a){
	// Reference force for an S_2-shaped particle
	double z = 2 * r / a;
	if (z > 2) return 1 / (r*r);
	else if (z > 1) return (12 / (z*z) - 224 + 896 * z - 840 * z*z + 224 * pow(z, 3) +
							70 * pow(z, 4) - 48 * pow(z, 5) + 7 * pow(z, 6)) / (35 * a*a);
	else return (224 * z - 224 * pow(z, 3) + 70 * pow(z, 4) + 48 * pow(z, 5) - 21 * pow(z, 7)) / (35 * a*a);
}

double force_tot(const double r, const double e2){
	return 1 / (r*r+e2);
}

void force_short(const Sim_Param &sim, const double D, const LinkedList& linked_list, Particle_v *particles,
				 const Vec_3D<double> position, Vec_3D<double>* force, Interp_obj* fs_interp)
{	// Calculate short range force in position, force is added
    #define FORCE_SHORT_NO_INTER
	int p;
	Vec_3D<double> dr_vec;
    double dr2;
    double dr; // <-- #ifdef FORCE_SHORT_NO_INTER
    const double m = pow(sim.box_opt.Ng, 3) / D;
    const unsigned Nm = sim.box_opt.mesh_num;
    const double rs2 = pow(sim.app_opt.rs, 2);
    const double e2 = pow(sim.box_opt.Ng*0.1, 2); // <-- #ifdef FORCE_SHORT_NO_INTER

    IT<3> it(position, sim.app_opt.Hc);
    do{
        p = linked_list.HOC(it.vec);
        while (p != -1){
            dr_vec = get_sgn_distance(particles[p].position, position, Nm);
            dr2 = dr_vec.norm2();
            if ((dr2 < rs2) && (dr2 != 0)) // Short range force is set 0 for separation larger than cutoff radius
            {
                #ifdef FORCE_SHORT_NO_INTER
                dr = sqrt(dr2);
                (*force) += dr_vec*(force_tot(dr, e2) - force_ref(dr, sim.app_opt.a))*m/(dr*4*PI);
                #else
                (*force) += dr_vec*(m/sqrt(dr2)*fs_interp->eval(dr2));
                #endif
            }
            p = linked_list.LL[p];
        }
    } while( it.iter() );
}

void kick_step_w_pp(const Sim_Param &sim, const double a, const double da, Particle_v* particles, const vector< Mesh> &force_field,
                    LinkedList* linked_list, Interp_obj* fs_interp)
{    // 2nd order ODE with long & short range potential
    const unsigned Np = sim.box_opt.par_num;
    Vec_3D<double> force;
    const double D = growth_factor(a, sim.cosmo);
    const double OL = sim.cosmo.Omega_L()*pow(a,3);
    const double Om = sim.cosmo.Omega_m;
    // -3/2a represents usual EOM, the rest are LCDM corrections
    const double f1 = 3/(2.*a)*(Om+2*OL)/(Om+OL);
    const double f2 = 3/(2.*a)*Om/(Om+OL)*D/a;
    
    printf("Creating linked list...\n");
	linked_list->get_linked_list(particles);

    cout << "Computing short and long range parts of the potential...\n";
    #pragma omp parallel for private(force)
    for (unsigned i = 0; i < Np; i++)
	{
        force.fill(0.);
        assign_from(force_field, particles[i].position, &force); // long-range force
        force_short(sim, D, *linked_list, particles, particles[i].position, &force, fs_interp); // short range force

        force = force*f2 - particles[i].velocity*f1;		
        particles[i].velocity += force*da;
    }
}

void upd_pos_first_order(const Sim_Param &sim, const double da, const double a, Particle_v* particles, const vector< Mesh> &vel_field)
{
    /// Leapfrog method for frozen-flow / adhesion
    stream_step(sim, da/2., particles);
    kick_step_no_momentum(sim, a-da/2., particles, vel_field);
    stream_step(sim, da/2., particles);
    get_per(particles, sim.box_opt.par_num, sim.box_opt.mesh_num);
}

void upd_pos_second_order(const Sim_Param &sim, const double da, const double a, Particle_v* particles, const vector< Mesh> &force_field)
{
    // Leapfrog method for frozen-potential
    stream_step(sim, da/2., particles);
    kick_step_w_momentum(sim, a-da/2., da, particles, force_field);
    stream_step(sim, da/2., particles);
    get_per(particles, sim.box_opt.par_num, sim.box_opt.mesh_num);
}

void upd_pos_second_order_w_pp(const Sim_Param &sim, const double da, const double a, Particle_v* particles, const vector< Mesh> &force_field,
                               LinkedList* linked_list, Interp_obj* fs_interp)
{
    // Leapfrog method for modified frozen-potential
    stream_step(sim, da/2., particles);
    kick_step_w_pp(sim, a-da/2., da, particles, force_field, linked_list, fs_interp);
    stream_step(sim, da/2., particles);
    get_per(particles, sim.box_opt.par_num, sim.box_opt.mesh_num);
}

static void gen_gauss_white_noise(const Sim_Param &sim, Mesh* rho)
{
	// Get keys for each slab in the x axis that this rank contains
	vector<unsigned long> slab_keys;
	slab_keys.resize(rho->N1);
	GetSlabKeys(slab_keys.data(), 0, rho->N1, sim.run_opt.seed);
	
	unsigned long ikey, index;
    double rn1, rn2, rn, tmp;
    const unsigned N = rho->N;
		
	#pragma omp parallel for private(ikey, index, rn1, rn2, rn, tmp)
	for(unsigned long i=0; i < N; ++i)
	{
		ikey = slab_keys[i];
		for(unsigned long j=0; j < N; ++j) 
		{
            #ifndef NOISE_HALF
			for(unsigned long k=0; k< N; ++k)
			{
                index = j*N + k; 
				GetRandomDoublesWhiteNoise(rn1, rn2, rn, ikey, index);
                tmp = sqrt(-2.0*log(rn)/rn);
                (*rho)(i, j, k) = rn2 * tmp;
            }
            #else
            for(unsigned long k=0; k < N/2; ++k) // go over half, use both random numbers
			{
                index = j*N + k;
				GetRandomDoublesWhiteNoise(rn1, rn2, rn, ikey, index);
                tmp = sqrt(-2.0*log(rn)/rn);
                (*rho)(i, j, 2*k) = rn2 * tmp;
                (*rho)(i, j, 2*k+1) = rn1 * tmp;
            }
            #endif
		}
    }     
    double t_mean;
	#ifdef CORR
	t_mean = mean(rho->real(), rho->length);
	double t_std_dev = std_dev(rho->real(), rho->length, t_mean);
	printf("\t[mean = %.12f, stdDev = %.12f]\t-->", t_mean, t_std_dev);
	(*rho)-=t_mean;
	(*rho)/=t_std_dev;
	#endif
	
	t_mean = mean(rho->real(), rho->length);
    printf("\t[mean = %.12f, stdDev = %.12f]\n", t_mean, std_dev(rho->real(), rho->length, t_mean));
    printf("\t[min = %.12f, max = %.12f]\n", min(rho->real(), rho->length), max(rho->real(), rho->length));
}

static void gen_rho_w_pow_k(const Sim_Param &sim, Mesh* rho)
{
    double k;
    const double L = sim.box_opt.box_size;
    const double k0 = 2.*PI/L;
    const int phase = sim.run_opt.phase ? 1 : -1;
    const unsigned N = rho->N;
    const unsigned len = rho->length / 2;
    const double mod = phase * pow(N / L, 3/2.); // pair sim, gaussian real -> fourier factor, dimension trans. Pk -> Pk*
    
	#pragma omp parallel for private(k)
	for(unsigned i=0; i < len; i++)
	{
        k = k0*sqrt(get_k_sq(N, i));
        (*rho)[2*i] *= mod*sqrt(lin_pow_spec(1, k, sim.cosmo));
        (*rho)[2*i+1] *= mod*sqrt(lin_pow_spec(1, k, sim.cosmo));
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
    fftw_execute_dft_r2c(p_F, *rho);

	printf("Generating density distributions with given power spectrum...\n");
	gen_rho_w_pow_k(sim, rho);
}

template <class T>
void get_rho_from_par(T* particles, Mesh* rho, const Sim_Param &sim)
{
    printf("Computing the density field from particle positions...\n");
    const double m = pow(sim.box_opt.Ng_pwr, 3.);
    const double mesh_mod = (double)sim.box_opt.mesh_num_pwr/sim.box_opt.mesh_num;
    const unsigned Np = sim.box_opt.par_num;

    rho->assign(-1.);
    
    #pragma omp parallel for
    for (unsigned i = 0; i < Np; i++)
    {
        assign_to(rho, particles[i].position*mesh_mod, m);
    }
}

int get_vel_from_par(Particle_v* particles, vector<Mesh>* vel_field, const Sim_Param &sim)
{
    printf("Computing the velocity field from particle positions...\n");
    const double mesh_mod = (double)sim.box_opt.mesh_num_pwr/sim.box_opt.mesh_num;
    const double m = pow(sim.box_opt.Ng_pwr, 3.);
    const unsigned Np = sim.box_opt.par_num;

    for(Mesh& field : *vel_field){
        field.assign(0.);
    }
    #pragma omp parallel for
    for (unsigned i = 0; i < Np; i++)
    {
        assign_to(vel_field, particles[i].position*mesh_mod, particles[i].velocity*(m*mesh_mod));
    }
    return 1;
}

int get_vel_from_par(Particle_x* particles, vector<Mesh>* vel_field, const Sim_Param &sim)
{
    printf("WARNING! Trying to compute velocity divergence with particle positions only! Skipping...\n");
    return 0;
}

void pwr_spec_k(const Mesh &rho_k, Mesh* power_aux)
{
    /* Computing the power spectrum P(k)/L^3 -- dimensionLESS!

    > in real part [even] of power_aux is stored pk, in imaginary [odd] dimensionLESS k
	> preserve values in rho_k
    > as power_aux can be Mesh of different (bigger) size than rho_k, all sizes / lengths are taken from rho_k
    */
	
	double w_k;
    Vec_3D<int> k_vec;
    const unsigned NM = rho_k.N;
    const unsigned half_length = rho_k.length / 2;

	#pragma omp parallel for private(w_k, k_vec)
	for(unsigned i=0; i < half_length;i++)
	{
		w_k = 1.;
		get_k_vec(NM, i, k_vec);
		for (int j = 0; j < 3; j++) if (k_vec[j] != 0) w_k *= pow(sin(k_vec[j]*PI/NM)/(k_vec[j]*PI/NM), ORDER + 1);
        (*power_aux)[2*i] = (rho_k[2*i]*rho_k[2*i] + rho_k[2*i+1]*rho_k[2*i+1])/(w_k*w_k);
		(*power_aux)[2*i+1] = k_vec.norm();
	}
}

void pwr_spec_k_init(const Mesh &rho_k, Mesh* power_aux)
{
    /* same as above but now there is NO w_k correction */

    Vec_3D<int> k_vec;
    const unsigned NM = rho_k.N;
    const unsigned half_length = rho_k.length / 2;

	#pragma omp parallel for private(k_vec)
	for(unsigned i=0; i < half_length;i++)
	{
		get_k_vec(NM, i, k_vec);
        (*power_aux)[2*i] = rho_k[2*i]*rho_k[2*i] + rho_k[2*i+1]*rho_k[2*i+1];
		(*power_aux)[2*i+1] = k_vec.norm();
	}
}

void vel_pwr_spec_k(const vector<Mesh> &vel_field, Mesh* power_aux)
{
    /* Computing the velocity power spectrum divergence P(k)/L^3 -- dimensionLESS!

    > in real part [even] of power_aux is stored pk, in imaginary [odd] dimensionLESS k
	> preserve values in rho_k
    > as power_aux can be Mesh of different (bigger) size than rho_k, all sizes / lengths are taken from rho_k
    */
	
	double w_k;
    Vec_3D<int> k_vec;

    const unsigned NM = vel_field[0].N;
    const unsigned half_length = vel_field[0].length / 2;

    double vel_div_re, vel_div_im, k; // temporary store of Pk in case vel_field[0] = power_aux

	#pragma omp parallel for private(w_k, k_vec, k, vel_div_re, vel_div_im)
	for(unsigned i=0; i < half_length; i++)
	{
        w_k = 1.;
        vel_div_re = vel_div_im = 0;
		get_k_vec(NM, i, k_vec);
        for (int j = 0; j < 3; j++){
            k = k_vec[j]*2*PI / NM;
            if (k != 0) w_k *= pow(sin(k/2.)/(k/2.), ORDER + 1);
            vel_div_re += vel_field[j][2*i]*k; // do not care about Re <-> Im in 2*PI*i/N, norm only
            vel_div_im += vel_field[j][2*i+1]*k;
        } 
        
        (*power_aux)[2*i] = (vel_div_re*vel_div_re + vel_div_im*vel_div_im)/(w_k*w_k);
		(*power_aux)[2*i+1] = k_vec.norm();
	}
}

void gen_cqty_binned(const double x_min, const double x_max, const unsigned bins_per_decade,
                    const Mesh &qty_mesh, const unsigned half_length, Data_Vec<double,2>& qty_binned, const double mod_q, const double mod_x)
{
    /* bin some complex quantity on mesh in logarithmic bins, assuming:
       Q(x) = mod_q*qty_mesh[2*i]
       x = mod_x*qty_mesh[2*i+1]
       [mesh[2*i+1]] = [x_min] = [x_max]

       return binned data in qty_binned {x, <Q(x)>}

       Note: passing length of the array for case when my mesh is bigger than data stored in there
             overloaded function exists when this is not the case
    */

    unsigned req_size = (unsigned)ceil(bins_per_decade*log10(x_max/x_min));
    qty_binned.resize(req_size);
    qty_binned.fill(0);
    vector<unsigned> tmp(req_size, 0); // for counts in bins

    double x;
    int bin;
    
    /* compute sum x, Q(x), Q^2(x) in bins */
    #pragma omp parallel for private(x, bin)
    for (unsigned i = 0; i < half_length; i++){
        x = qty_mesh[2*i+1];
        if ((x <x_max) && (x>=x_min)){
            bin = (int)((log10(x) - log10(x_min)) * bins_per_decade);
            #pragma omp atomic
            qty_binned[0][bin] += x;
            #pragma omp atomic
            qty_binned[1][bin] += qty_mesh[2*i];
            #pragma omp atomic
            tmp[bin]++;
        }
    }

    /* compute average x, Q(x) in bins */
    unsigned count;
    for (unsigned j = 0; j < qty_binned.size(); ){
        count = tmp[j];
        if (count){
            qty_binned[0][j] *= mod_x / count;
            qty_binned[1][j] *= mod_q / count;
            j++;
        }else{
            qty_binned.erase(j);
            tmp.erase(tmp.begin() + j);
        }
	}
}


void gen_cqty_binned(const double x_min, const double x_max, const unsigned bins_per_decade,
                    const Mesh &qty_mesh, Data_Vec<double, 2>& qty_binned, const double mod_q, const double mod_x)
{
    gen_cqty_binned(x_min, x_max, bins_per_decade, qty_mesh, qty_mesh.length / 2, qty_binned, mod_q, mod_x);
}
/*
void gen_rqty_binned(const double x_min, const double x_max, const double x_0,
    const Mesh &qty_mesh, Data_Vec<double, 2>& qty_binned, const double mod_q)
{
/* bin some real quantity on mesh in linear bins
   x_min and x_max are in [Mpc/h], transformation from mesh distances to Mpc/h assumed to be x_0
   Q(x) = qty_mesh[i]
   x = ([i,j,k].norm())*x_0
   
*//*
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
                    qty_binned[1][bin] += qty_mesh(i,j,k);
                    #pragma omp atomic
                    qty_binned[0][bin]++;
                }
            }
        }
    }
    const double x_min_ = x_min+lin_bin/2;
    unsigned i = 0;
    for (unsigned j = 0; j < qty_binned.size(); ){
        if (qty_binned[0][j]){
            qty_binned[1][j] *= mod_q / qty_binned[0][j];
            x = x_min_ + lin_bin*i;
            qty_binned[0][j] = x;
            j++;
        }else{
            qty_binned.erase(j);
        }
        i++;
    }
}
*/

void gen_pow_spec_binned(const Sim_Param &sim, const Mesh &power_aux, Data_Vec<double, 2>* pwr_spec_binned)
{
    const double mod_pk = pow(sim.box_opt.box_size, 3.); // P(k) -> dimensionFULL!
    const double mod_k = 2.*PI/sim.box_opt.box_size;
    printf("Computing binned power spectrum...\n");
	gen_cqty_binned(1, sim.box_opt.mesh_num_pwr,  sim.out_opt.bins_per_decade, power_aux, *pwr_spec_binned, mod_pk, mod_k);
}

void gen_pow_spec_binned_init(const Sim_Param &sim, const Mesh &power_aux, const unsigned half_length, Data_Vec<double, 2>* pwr_spec_binned)
{
    /* same as above but now  power_aux is storing only data [0...mesh_num], NOT mesh_num_pwr */
    const double mod_pk = pow(sim.box_opt.box_size, 3.); // P(k) -> dimensionFULL!
    const double mod_k = 2.*PI/sim.box_opt.box_size;
    printf("Computing binned initial power spectrum...\n");
	gen_cqty_binned(1, sim.box_opt.mesh_num,  sim.out_opt.bins_per_decade, power_aux, half_length, *pwr_spec_binned, mod_pk, mod_k);
}

void gen_pow_spec_binned_from_extrap(const Sim_Param &sim, const Extrap_Pk &P_k, Data_Vec<double, 2>* pwr_spec_binned)
{
    const double k_max = sim.other_par.k_print.upper;
    const double k_min = sim.other_par.k_print.lower;
    const double log_bin = 1./ sim.out_opt.bins_per_decade;
    double k;
    unsigned req_size = (unsigned)ceil( sim.out_opt.bins_per_decade*log10(k_max/k_min));
    pwr_spec_binned->resize(req_size);

    #pragma omp parallel for private(k)
	for (unsigned j = 0; j < pwr_spec_binned->size(); j++){
        k = k_min*exp10(j*log_bin);
		(*pwr_spec_binned)[0][j] = k;
        (*pwr_spec_binned)[1][j] = P_k(k);
    }
}


void gen_pot_k(const Mesh& rho_k, Mesh* pot_k)
{   /*
    pot_k can be Mesh of differen (bigger) size rho_k,
    !!!> ALL physical FACTORS ARE therefore TAKEN FROM rho_k <!!!
    */
	printf("Computing potential in k-space...\n");
    double k2;
    const unsigned N = rho_k.N; // for case when pot_k is different mesh than vel_field
    const double d2_k = pow(2.*PI/N, 2.); // factor from second derivative with respect to the mesh coordinates
    const unsigned l_half = rho_k.length/2;

	#pragma omp parallel for private(k2)
	for(unsigned i=0; i < l_half;i++){				
		k2 = get_k_sq(N, i);
		if (k2 == 0){
			(*pot_k)[2*i] = 0;
			(*pot_k)[2*i+1] = 0;
		} else{
			(*pot_k)[2*i] = -rho_k[2*i]/(k2*d2_k);
			(*pot_k)[2*i+1] = -rho_k[2*i+1]/(k2*d2_k);
		}
	}
}

void gen_pot_k(Mesh* rho_k){ gen_pot_k(*rho_k, rho_k); }

static double S2_shape(const double k2, const double a)
{
	if (a == 0) return 1.;
	
    double t = sqrt(k2)*a / 2.;
    if (t == 0) return 1.;
	return 12 / pow(t, 4)*(2 - 2 * cos(t) - t*sin(t));
}

static double CIC_opt(Vec_3D<double> k_vec, const double a)
{
#define N_MAX 1
#ifndef N_MAX
    double s2 = pow(S2_shape(k_vec.norm2(), a), 2);
    for(int j=0; j<3; j++)
    {
        if (k_vec[j] != 0) s2 /= pow(sin(k_vec[j] / 2.) / (k_vec[j] / 2.), 2); //W (k) for CIC (order 1)
    }
    return s2;
#else
	double k_n[3];
	double U2, U_n, G_n, k2n;
	
	G_n = 0;
	U2 = 1;
	for (int j = 0; j < 3; j++) U2 *= 1./3.*(1+2*pow(cos(k_vec[j]/2.), 2)); // inf sum of U_n^2 for CIC
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
                U_n = pow(U_n, 2.); // W(k) for CIC (order 1)
				if (k2n != 0)
				{
					for(int j=0; j<3; j++)
					{										
						G_n += k_vec[j]* // i.D(k)
						k_n[j]/k2n*pow(S2_shape(k2n, a), 2)* // R*(k_n) /i
						U_n; // U_n
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
#endif
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
    
    const unsigned N = (*vel_field)[0].N; // for case when pot_k is different mesh than vel_field
    const double d_k = 2.*PI/N;  // 2*PI/N comes from derivative WITH RESPECT to the mesh coordinates
    const unsigned l_half = (*vel_field)[0].length/2;
	
	#pragma omp parallel for private(opt, k_vec, k_vec_phys, potential_tmp)
	for(unsigned i=0; i < l_half;i++)
	{
		potential_tmp[0] = pot_k[2*i]; // prevent overwriting if vel_field[0] == pot_k
        potential_tmp[1] = pot_k[2*i+1]; // prevent overwriting if vel_field[0] == pot_k
        get_k_vec(N, i, k_vec);
        k_vec_phys = Vec_3D<double>(k_vec)*d_k;	
        // no optimalization
        if (a == -1) opt = 1.;
        // optimalization for CIC and S2 shaped particle
        else opt = CIC_opt(k_vec_phys, a);
		for(unsigned j=0; j<3;j++)
		{
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
    const unsigned Ng_pwr = sim.box_opt.Ng_pwr;
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
                //dens_binned[bin] += pow(sim.box_opt.Ng_pwr, 3);
			}
		}
	}
}

template void get_rho_from_par(Particle_x*, Mesh*, const Sim_Param&);
template void get_rho_from_par(Particle_v*, Mesh*, const Sim_Param&);