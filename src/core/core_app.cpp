/**
 * @brief implementation for common functions for all types of approximations
 * 
 * @file core_app.cpp
 * @author Michal Vrastil
 * @date 2018-07-11
 */

#include <algorithm>
#include "core_app.h"
#include "core_mesh.h"
#include "CBRNG_Random.h"

#ifndef ORDER
#define ORDER 1
#endif

template<typename T>
static T mean(const std::vector<T>& data)
{
    T tmp(0);
	
	#pragma omp parallel for reduction(+:tmp)
	for (auto it = data.begin(); it < data.end(); ++it) tmp += *it;
	
	return tmp / data.size();
}

static FTYPE_t mean(const Mesh& data)
{
    return mean(data.data);
}

template<typename T>
static T std_dev(const std::vector<T>& data, T mean)
{
    T tmp(0);
	
	#pragma omp parallel for reduction(+:tmp)
	for (auto it = data.begin(); it < data.end(); ++it) tmp += pow2(*it-mean);
	
	return sqrt(tmp / data.size());
}

static FTYPE_t std_dev(const Mesh& data, FTYPE_t mean)
{
    return std_dev(data.data, mean);
}

template<typename T>
static T min(const std::vector<T>& data)
{
    return *std::min_element(data.begin(), data.end());
}

static FTYPE_t min(const Mesh& data){
    return min(data.data);
}

template<typename T>
static T max(const std::vector<T>& data)
{
    return *std::max_element(data.begin(), data.end());
}

static FTYPE_t max(const Mesh& data){
    return max(data.data);
}

static void set_unpert_pos_one_par(Vec_3D<size_t>& unpert_pos, const size_t par_index, const size_t par_per_dim, const size_t Ng)
{
	unpert_pos[0] = (par_index / (par_per_dim * par_per_dim)) * Ng;
	unpert_pos[1] = ((par_index / par_per_dim) % par_per_dim) * Ng;
	unpert_pos[2] = (par_index % par_per_dim) * Ng;
}

static void set_velocity_one_par(const Vec_3D<size_t>& unpert_pos, Vec_3D<FTYPE_t>& displ_field, const std::vector<Mesh> &vel_field)
{
	for (size_t i = 0; i < 3; i++) displ_field[i] = vel_field[i](unpert_pos);
}

void set_unpert_pos(const Sim_Param &sim, std::vector<Particle_x<FTYPE_t>>& particles)
{
	Vec_3D<size_t> unpert_pos;
    const size_t par_per_dim = sim.box_opt.par_num_1d;
    const size_t Ng = sim.box_opt.Ng;
    const size_t Np = sim.box_opt.par_num;
	
	#pragma omp parallel for private(unpert_pos)
	for(size_t i=0; i< Np; i++)
	{
		set_unpert_pos_one_par(unpert_pos, i, par_per_dim, Ng);		
		particles[i] = Particle_x<FTYPE_t>(unpert_pos);
	}
}

void set_unpert_pos_w_vel(const Sim_Param &sim, std::vector<Particle_v<FTYPE_t>>& particles, const std::vector<Mesh> &vel_field)
{
	Vec_3D<size_t> unpert_pos;
	Vec_3D<FTYPE_t> velocity;
	const size_t par_per_dim = sim.box_opt.par_num_1d;
    const size_t Ng = sim.box_opt.Ng;
    
    const size_t Np = sim.box_opt.par_num;
	#pragma omp parallel for private(unpert_pos, velocity)
	for(size_t i=0; i< Np; i++)
	{
		set_unpert_pos_one_par(unpert_pos, i, par_per_dim, Ng);
		set_velocity_one_par(unpert_pos, velocity, vel_field);
		particles[i] = Particle_v<FTYPE_t>(unpert_pos, velocity);
	}
}

void set_pert_pos(const Sim_Param &sim, const FTYPE_t db, std::vector<Particle_x<FTYPE_t>>& particles, const std::vector< Mesh> &vel_field)
{
    BOOST_LOG_TRIVIAL(debug) << "Setting initial positions of particles...";
	Vec_3D<size_t> unpert_pos;
	Vec_3D<FTYPE_t> displ_field;
	Vec_3D<FTYPE_t> pert_pos;
	
    const size_t par_per_dim = sim.box_opt.par_num_1d;
    const size_t Ng = sim.box_opt.Ng;
    const size_t Nm = sim.box_opt.mesh_num;
    const size_t Np = sim.box_opt.par_num;
	
	#pragma omp parallel for private(unpert_pos, displ_field, pert_pos)
	for(size_t i=0; i< Np; i++)
	{
		set_unpert_pos_one_par(unpert_pos, i, par_per_dim, Ng);
		set_velocity_one_par(unpert_pos, displ_field, vel_field);
		pert_pos = displ_field*db + unpert_pos;
		get_per(pert_pos, Nm);
		particles[i] = Particle_x<FTYPE_t>(pert_pos);		
	}
}

void set_pert_pos(const Sim_Param &sim, const FTYPE_t a, std::vector<Particle_v<FTYPE_t>>& particles, const std::vector< Mesh> &vel_field)
{
    BOOST_LOG_TRIVIAL(debug) << "Setting initial positions and velocitis of particles...";
	Vec_3D<size_t> unpert_pos;
	Vec_3D<FTYPE_t> velocity;
	Vec_3D<FTYPE_t> pert_pos;
	
	const size_t par_per_dim = sim.box_opt.par_num_1d;
	const size_t Ng = sim.box_opt.Ng;
    const size_t Nm = sim.box_opt.mesh_num;
    const size_t Np = sim.box_opt.par_num;

    const FTYPE_t D = growth_factor(a, sim.cosmo); // growth factor
    const FTYPE_t dDda = growth_change(a, sim.cosmo); // dD / da

	#pragma omp parallel for private(unpert_pos, velocity, pert_pos)
	for(size_t i=0; i< Np; i++)
	{
		set_unpert_pos_one_par(unpert_pos, i, par_per_dim, Ng);
		set_velocity_one_par(unpert_pos, velocity, vel_field);
		pert_pos = velocity*D + unpert_pos;
		get_per(pert_pos, Nm);
		particles[i] = Particle_v<FTYPE_t>(pert_pos, velocity*dDda);		
	}
}

static void gen_gauss_white_noise(const Sim_Param &sim, Mesh& rho)
{
	// Get keys for each slab in the x axis that this rank contains
	std::vector<size_t> slab_keys;
	slab_keys.resize(rho.N1);
	GetSlabKeys(slab_keys.data(), 0, rho.N1, sim.run_opt.seed);
	
	size_t ikey, index;
    FTYPE_t rn1, rn2, rn, tmp;
    const size_t N = rho.N;
		
	#pragma omp parallel for private(ikey, index, rn1, rn2, rn, tmp)
	for(size_t i=0; i < N; ++i)
	{
		ikey = slab_keys[i];
		for(size_t j=0; j < N; ++j) 
		{
            #ifndef NOISE_HALF
			for(size_t k=0; k< N; ++k)
			{
                index = j*N + k; 
				GetRandomDoublesWhiteNoise(rn1, rn2, rn, ikey, index);
                tmp = sqrt(-2*log(rn)/rn);
                rho(i, j, k) = rn2 * tmp;
            }
            #else
            for(size_t k=0; k < N/2; ++k) // go over half, use both random numbers
			{
                index = j*N + k;
				GetRandomDoublesWhiteNoise(rn1, rn2, rn, ikey, index);
                tmp = sqrt(-2*log(rn)/rn);
                rho(i, j, 2*k) = rn2 * tmp;
                rho(i, j, 2*k+1) = rn1 * tmp;
            }
            #endif
		}
    }     
    FTYPE_t t_mean;
	#ifdef CORR
	t_mean = mean(rho);
	FTYPE_t t_std_dev = std_dev(rho, t_mean);
	BOOST_LOG_TRIVIAL(debug) << "\t[mean = " << t_mean << ", stdDev = " << t_std_dev << "]\t-->";
	rho-=t_mean;
	rho/=t_std_dev;
	#endif
	
	t_mean = mean(rho);
    BOOST_LOG_TRIVIAL(debug) << "\t[mean = " << t_mean << ", stdDev = " << std_dev(rho, t_mean) << "]\t<--";
    BOOST_LOG_TRIVIAL(debug) << "\t[min = " << min(rho) << ", max = " << max(rho) << "]";
}

static FTYPE_t truncation_fce(FTYPE_t k, FTYPE_t k2_G)
{
    return exp(-k*k/k2_G);
}

static void gen_rho_w_pow_k(const Sim_Param &sim, Mesh& rho)
{
    // extract const variables
    const FTYPE_t L = sim.box_opt.box_size;
    const FTYPE_t k0 = 2*PI/L;
    const int phase = sim.run_opt.phase ? 1 : -1;
    const size_t N = rho.N;
    const size_t len = rho.length / 2;
    const FTYPE_t mod = phase * pow(N / L, 3/2.); // pair sim, gaussian real -> fourier factor, dimension trans. Pk -> Pk*
    const bool truncation = sim.cosmo.truncated_pk;

    if(truncation)
    {
        const FTYPE_t k2_G = sim.cosmo.k2_G;
        #pragma omp parallel for
        for(size_t i=0; i < len; i++)
        {
            FTYPE_t k = k0*sqrt(get_k_sq(N, i));
            FTYPE_t tr = truncation_fce(k, k2_G);
            rho[2*i] *= mod*sqrt(lin_pow_spec(1, k, sim.cosmo))*tr;
            rho[2*i+1] *= mod*sqrt(lin_pow_spec(1, k, sim.cosmo))*tr;
        }
    }
    else
    {
        #pragma omp parallel for
        for(size_t i=0; i < len; i++)
        {
            FTYPE_t k = k0*sqrt(get_k_sq(N, i));
            rho[2*i] *= mod*sqrt(lin_pow_spec(1, k, sim.cosmo));
            rho[2*i+1] *= mod*sqrt(lin_pow_spec(1, k, sim.cosmo));
        }
    }
    
}

/**
 * @brief Generate density distributions \f$\rho(k)\f$ in k-space.
 * 
 * @param sim 
 * @param rho 
 * @param p_F 
 * 
 * At first, a gaussian white noise (mean = 0, stdDev = 1) is generated,
 * then it is convoluted with given power spectrum.
 */
void gen_rho_dist_k(const Sim_Param &sim, Mesh& rho, const FFTW_PLAN_TYPE &p_F)
{
	BOOST_LOG_TRIVIAL(debug) << "Generating gaussian white noise...";
	gen_gauss_white_noise(sim, rho);
    fftw_execute_dft_r2c(p_F, rho);

	BOOST_LOG_TRIVIAL(debug) << "Generating density distributions with given power spectrum...";
	gen_rho_w_pow_k(sim, rho);
}

template <class T>
void get_rho_from_par(const std::vector<T>& particles, Mesh& rho, const Sim_Param &sim)
{
    BOOST_LOG_TRIVIAL(debug) << "Computing the density field from particle positions...";

    const size_t Np = sim.box_opt.par_num;
    if (particles.size() != Np){
        std::string msg = "Number of particles (" + std::to_string(particles.size()) + ") does not correspond with simulation parameters (" + std::to_string(Np) + ")!";
        throw std::range_error(msg);
    }
    const FTYPE_t m = pow((FTYPE_t)rho.N, 3) / Np;
    const FTYPE_t mesh_mod = (FTYPE_t)rho.N/sim.box_opt.mesh_num;
    

    rho.assign(-1.);
    
    #pragma omp parallel for
    for (size_t i = 0; i < Np; i++)
    {
        assign_to(rho, particles[i].position*mesh_mod, m);
    }
}

bool get_vel_from_par(const std::vector<Particle_v<FTYPE_t>>& particles, std::vector<Mesh>& vel_field, const Sim_Param &sim)
{
    BOOST_LOG_TRIVIAL(debug) << "Computing the velocity field from particle positions...";
    const FTYPE_t mesh_mod = (FTYPE_t)sim.box_opt.mesh_num_pwr/sim.box_opt.mesh_num;
    const FTYPE_t m = pow((FTYPE_t)sim.box_opt.Ng_pwr, 3);
    const size_t Np = sim.box_opt.par_num;

    for(Mesh& field : vel_field){
        field.assign(0.);
    }
    #pragma omp parallel for
    for (size_t i = 0; i < Np; i++)
    {
        assign_to(vel_field, particles[i].position*mesh_mod, particles[i].velocity*(m*mesh_mod));
    }
    return true;
}

bool get_vel_from_par(const std::vector<Particle_x<FTYPE_t>>& particles, std::vector<Mesh>& vel_field, const Sim_Param &sim)
{
    BOOST_LOG_TRIVIAL(debug) << "WARNING! Trying to compute velocity divergence with particle positions only! Skipping...";
    return false;
}

void pwr_spec_k(const Mesh &rho_k, Mesh& power_aux)
{
    /* Computing the power spectrum P(k)/L^3 -- dimensionLESS!

    > in real part [even] of power_aux is stored pk, in imaginary [odd] dimensionLESS k
	> preserve values in rho_k
    > as power_aux can be Mesh of different (bigger) size than rho_k, all sizes / lengths are taken from rho_k
    */
	
	FTYPE_t w_k;
    Vec_3D<int> k_vec;
    const size_t NM = rho_k.N;
    const size_t half_length = rho_k.length / 2;

	#pragma omp parallel for private(w_k, k_vec)
	for(size_t i=0; i < half_length;i++)
	{
		w_k = 1.;
		get_k_vec(NM, i, k_vec);
		for (unsigned int j = 0; j < 3; j++) if (k_vec[j] != 0) w_k *= pow(sin(k_vec[j]*PI/NM)/(k_vec[j]*PI/NM), ORDER + 1);
        power_aux[2*i] = (rho_k[2*i]*rho_k[2*i] + rho_k[2*i+1]*rho_k[2*i+1])/(w_k*w_k);
		power_aux[2*i+1] = k_vec.norm();
	}
}

void pwr_spec_k_init(const Mesh &rho_k, Mesh& power_aux)
{
    /* same as above but now there is NO w_k correction */

    Vec_3D<int> k_vec;
    const size_t NM = rho_k.N;
    const size_t half_length = rho_k.length / 2;

	#pragma omp parallel for private(k_vec)
	for(size_t i=0; i < half_length;i++)
	{
		get_k_vec(NM, i, k_vec);
        power_aux[2*i] = pow2(rho_k[2*i]) + pow2(rho_k[2*i+1]);
		power_aux[2*i+1] = k_vec.norm();
	}
}

void vel_pwr_spec_k(const std::vector<Mesh> &vel_field, Mesh& power_aux)
{
    /* Computing the velocity power spectrum divergence P(k)/L^3 -- dimensionLESS!

    > in real part [even] of power_aux is stored pk, in imaginary [odd] dimensionLESS k
	> preserve values in rho_k
    > as power_aux can be Mesh of different (bigger) size than rho_k, all sizes / lengths are taken from rho_k
    */
	
	FTYPE_t w_k;
    Vec_3D<int> k_vec;

    const size_t NM = vel_field[0].N;
    const size_t half_length = vel_field[0].length / 2;

    FTYPE_t vel_div_re, vel_div_im, k; // temporary store of Pk in case vel_field[0] = power_aux

	#pragma omp parallel for private(w_k, k_vec, k, vel_div_re, vel_div_im)
	for(size_t i=0; i < half_length; i++)
	{
        w_k = 1.;
        vel_div_re = vel_div_im = 0;
		get_k_vec(NM, i, k_vec);
        for (unsigned int j = 0; j < 3; j++){
            k = k_vec[j]*2*PI / NM;
            if (k != 0) w_k *= pow(sin(k/2)/(k/2), ORDER + 1);
            vel_div_re += vel_field[j][2*i]*k; // do not care about Re <-> Im in 2*PI*i/N, norm only
            vel_div_im += vel_field[j][2*i+1]*k;
        } 
        
        power_aux[2*i] = (vel_div_re*vel_div_re + vel_div_im*vel_div_im)/(w_k*w_k);
		power_aux[2*i+1] = k_vec.norm();
	}
}

void gen_cqty_binned(const FTYPE_t x_min, const FTYPE_t x_max, const size_t bins_per_decade,
                    const Mesh &qty_mesh, const size_t half_length, Data_Vec<FTYPE_t,2>& qty_binned, const FTYPE_t mod_q, const FTYPE_t mod_x)
{
    /* bin some complex quantity on mesh in logarithmic bins, assuming:
       Q(x) = mod_q*qty_mesh[2*i]
       x = mod_x*qty_mesh[2*i+1]
       [mesh[2*i+1]] = [x_min] = [x_max]

       return binned data in qty_binned {x, <Q(x)>}

       Note: passing length of the array for case when my mesh is bigger than data stored in there
             overloaded function exists when this is not the case
    */

    size_t req_size = (size_t)ceil(bins_per_decade*log10(x_max/x_min));
    qty_binned.resize(req_size);
    qty_binned.fill(0);
    std::vector<size_t> tmp(req_size, 0); // for counts in bins

    FTYPE_t x;
    size_t bin;
    
    /* compute sum x, Q(x), Q^2(x) in bins */
    #pragma omp parallel for private(x, bin)
    for (size_t i = 0; i < half_length; i++){
        x = qty_mesh[2*i+1];
        if ((x <x_max) && (x>=x_min)){
            bin = (size_t)((log10(x) - log10(x_min)) * bins_per_decade);
            #pragma omp atomic
            qty_binned[0][bin] += x;
            #pragma omp atomic
            qty_binned[1][bin] += qty_mesh[2*i];
            #pragma omp atomic
            tmp[bin]++;
        }
    }

    /* compute average x, Q(x) in bins */
    size_t count;
    for (size_t j = 0; j < qty_binned.size(); ){
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


void gen_cqty_binned(const FTYPE_t x_min, const FTYPE_t x_max, const size_t bins_per_decade,
                    const Mesh &qty_mesh, Data_Vec<FTYPE_t, 2>& qty_binned, const FTYPE_t mod_q, const FTYPE_t mod_x)
{
    gen_cqty_binned(x_min, x_max, bins_per_decade, qty_mesh, qty_mesh.length / 2, qty_binned, mod_q, mod_x);
}

void gen_pow_spec_binned(const Sim_Param &sim, const Mesh &power_aux, Data_Vec<FTYPE_t, 2>& pwr_spec_binned)
{
    const FTYPE_t mod_pk = pow(sim.box_opt.box_size, 3); // P(k) -> dimensionFULL!
    const FTYPE_t mod_k = 2*PI/sim.box_opt.box_size;
    BOOST_LOG_TRIVIAL(debug) << "Computing binned power spectrum...";
	gen_cqty_binned(1, sim.box_opt.mesh_num_pwr,  sim.out_opt.bins_per_decade, power_aux, pwr_spec_binned, mod_pk, mod_k);
}

void gen_pow_spec_binned_init(const Sim_Param &sim, const Mesh &power_aux, const size_t half_length, Data_Vec<FTYPE_t, 2>& pwr_spec_binned)
{
    /* same as above but now  power_aux is storing only data [0...mesh_num], NOT mesh_num_pwr */
    const FTYPE_t mod_pk = pow(sim.box_opt.box_size, 3); // P(k) -> dimensionFULL!
    const FTYPE_t mod_k = 2*PI/sim.box_opt.box_size;
    BOOST_LOG_TRIVIAL(debug) << "Computing binned initial power spectrum...";
	gen_cqty_binned(1, sim.box_opt.mesh_num,  sim.out_opt.bins_per_decade, power_aux, half_length, pwr_spec_binned, mod_pk, mod_k);
}

template<class P, typename T, size_t N> // P = everything callable P_k(k), T = float-type, N = number
void gen_pow_spec_binned_from_extrap(const Sim_Param &sim, const P &P_k, Data_Vec<T, N>& pwr_spec_binned)
{
    const T k_max = sim.other_par.k_print.upper;
    const T k_min = sim.other_par.k_print.lower;
    const T log_bin = T(1) / sim.out_opt.bins_per_decade;
    T k;
    size_t req_size = (size_t)ceil( sim.out_opt.bins_per_decade*log10(k_max/k_min));
    pwr_spec_binned.resize(req_size);

    #pragma omp parallel for private(k)
	for (size_t j = 0; j < pwr_spec_binned.size(); j++){
        k = k_min*pow(T(10), j*log_bin);
		pwr_spec_binned[0][j] = k;
        pwr_spec_binned[1][j] = P_k(k);
    }
}


void gen_pot_k(const Mesh& rho_k, Mesh& pot_k)
{   /*
    pot_k can be Mesh of differen (bigger) size rho_k,
    !!!> ALL physical FACTORS ARE therefore TAKEN FROM rho_k <!!!
    */
	BOOST_LOG_TRIVIAL(debug) << "Computing potential in k-space...";
    FTYPE_t k2;
    const size_t N = rho_k.N; // for case when pot_k is different mesh than vel_field
    const FTYPE_t d2_k = pow2(2*PI/N); // factor from second derivative with respect to the mesh coordinates
    const size_t l_half = rho_k.length/2;

	#pragma omp parallel for private(k2)
	for(size_t i=0; i < l_half;i++){				
		k2 = get_k_sq(N, i);
		if (k2 == 0){
			pot_k[2*i] = 0;
			pot_k[2*i+1] = 0;
		} else{
			pot_k[2*i] = -rho_k[2*i]/(k2*d2_k);
			pot_k[2*i+1] = -rho_k[2*i+1]/(k2*d2_k);
		}
	}
}

void gen_pot_k(Mesh& rho_k){ gen_pot_k(rho_k, rho_k); }

static FTYPE_t S2_shape(const FTYPE_t k2, const FTYPE_t a)
{
	if (a == 0) return 1.;
	
    FTYPE_t t = sqrt(k2)*a / 2;
    if (t == 0) return 1.;
	return 12 / pow(t, 4)*(2 - 2 * cos(t) - t*sin(t));
}

static FTYPE_t CIC_opt(Vec_3D<FTYPE_t> k_vec, const FTYPE_t a)
{
#define N_MAX 1
#ifndef N_MAX
    FTYPE_t s2 = pow2(S2_shape(k_vec.norm2(), a));
    for(unsigned int j=0; j<3; j++)
    {
        if (k_vec[j] != 0) s2 /= pow2(sin(k_vec[j] / 2) / (k_vec[j] / 2)); //W (k) for CIC (order 1)
    }
    return s2;
#else
	FTYPE_t k_n[3];
	FTYPE_t U2, U_n, G_n, k2n;
	
	G_n = 0;
	U2 = 1;
	for (unsigned int j = 0; j < 3; j++) U2 *= (1+2*pow2(cos(k_vec[j]/2)))/3; // inf sum of U_n^2 for CIC
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
				for(unsigned int j=0; j<3; j++)
				{
					if (k_n[j] != 0) U_n *= sin(k_n[j] / 2) / (k_n[j] / 2);
					k2n += pow2(k_n[j]);
                }
                U_n = pow2(U_n); // W(k) for CIC (order 1)
				if (k2n != 0)
				{
					for(unsigned int j=0; j<3; j++)
					{										
						G_n += k_vec[j]* // i.D(k)
						k_n[j]/k2n*pow2(S2_shape(k2n, a))* // R*(k_n) /i
						U_n; // U_n
					}
				}
			}
		}
	}
	if ((G_n != G_n) || (U2 != U2))
	{
		BOOST_LOG_TRIVIAL(warning) << "Gn = " << G_n << "\tU2 = " << U2 << ", k = (" << k_vec[0] << ", " << k_vec[1] << ", " << k_vec[2] << ")";
		return 1.;
	}
    return G_n/U2;
#endif
}

void gen_displ_k_S2(std::vector<Mesh>& vel_field, const Mesh& pot_k, const FTYPE_t a)
{   /*
    pot_k can be Mesh of differen (bigger) size than each vel_field,
    !!!> ALL physical FACTORS ARE therefore TAKEN FROM vel_field[0] <!!!
    */
	if (a == -1) BOOST_LOG_TRIVIAL(debug) << "Computing displacement in k-space...";
	else if (a == 0) BOOST_LOG_TRIVIAL(debug) << "Computing displacement in k-space with CIC opt...";
	else BOOST_LOG_TRIVIAL(debug) << "Computing force in k-space for S2 shaped particles with CIC opt...";

	FTYPE_t opt;
    Vec_3D<int> k_vec;
    Vec_3D<FTYPE_t> k_vec_phys;
    FTYPE_t potential_tmp[2];
    
    const size_t N = vel_field[0].N; // for case when pot_k is different mesh than vel_field
    const FTYPE_t d_k = 2*PI/N;  // 2*PI/N comes from derivative WITH RESPECT to the mesh coordinates
    const size_t l_half = vel_field[0].length/2;
	
	#pragma omp parallel for private(opt, k_vec, k_vec_phys, potential_tmp)
	for(size_t i=0; i < l_half;i++)
	{
		potential_tmp[0] = pot_k[2*i]; // prevent overwriting if vel_field[0] == pot_k
        potential_tmp[1] = pot_k[2*i+1]; // prevent overwriting if vel_field[0] == pot_k
        get_k_vec(N, i, k_vec);
        k_vec_phys = d_k*k_vec;	
        // no optimalization
        if (a == -1) opt = 1.;
        // optimalization for CIC and S2 shaped particle
        else opt = CIC_opt(k_vec_phys, a);
		for(size_t j=0; j<3;j++)
		{
			vel_field[j][2*i] = k_vec_phys[j]*potential_tmp[1]*opt;
			vel_field[j][2*i+1] = -k_vec_phys[j]*potential_tmp[0]*opt;
		}
	}
}

void gen_displ_k(std::vector<Mesh>& vel_field, const Mesh& pot_k) {gen_displ_k_S2(vel_field, pot_k, -1);}

void gen_displ_k_cic(std::vector<Mesh>& vel_field, const Mesh& pot_k) {gen_displ_k_S2(vel_field, pot_k, 0.);}

void gen_dens_binned(const Mesh& rho, std::vector<size_t> &dens_binned, const Sim_Param &sim)
{
	BOOST_LOG_TRIVIAL(debug) << "Computing binned density field...";
	size_t bin;
    FTYPE_t rho_avg;
    const size_t Ng_pwr = sim.box_opt.Ng_pwr;
    const size_t N = rho.N;

	dens_binned.assign(dens_binned.size(), 0);
    
    #pragma omp parallel for private(bin, rho_avg)
	for (size_t i = 0; i < N; i+=Ng_pwr)
	{
		for (size_t j = 0; j < N; j+=Ng_pwr)
		{
			for (size_t k = 0; k < N; k+=Ng_pwr)
			{
				// Need to go through all mesh cells [i, i+Ng-1]*[j, j+Ng-1], [k, k+Ng, -1]
				rho_avg = 0;
				for (size_t ii = i; ii < i+Ng_pwr; ii++)
				{
					for (size_t jj = j; jj  < j+Ng_pwr; jj++)
					{
						for (size_t kk = k; kk < k+Ng_pwr; kk++)
						{
							rho_avg+=rho(ii, jj, kk);
						}
					}
				}
				rho_avg /= pow((FTYPE_t)Ng_pwr, 3);
                bin = (size_t)((rho_avg+1)/0.1);
                if (bin >= dens_binned.size()) bin = dens_binned.size() - 1;
                // if (bin >= dens_binned.capacity()) dens_binned.resize(bin+1);
                #pragma omp atomic
                dens_binned[bin]++;
                //dens_binned[bin] += pow(sim.box_opt.Ng_pwr, 3);
			}
		}
	}
}

template void get_rho_from_par(const std::vector<Particle_x<FTYPE_t>>&, Mesh&, const Sim_Param&);
template void get_rho_from_par(const std::vector<Particle_v<FTYPE_t>>&, Mesh&, const Sim_Param&);
template void gen_pow_spec_binned_from_extrap(const Sim_Param&, const Extrap_Pk<FTYPE_t, 2>&, Data_Vec<FTYPE_t, 2>&);