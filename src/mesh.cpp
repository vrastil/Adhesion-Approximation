#include "mesh.h"
#include "stdafx.h"
#include "cmd_line.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>
#include "CBRNG_Random.h"
#include "threadpool.hpp"
#include <boost/filesystem.hpp>

// #define CORR
#define N_MAX 1

using namespace std;
using namespace boost::threadpool;
namespace fs = boost::filesystem;

const double PI = acos(-1.);

static const char *humanSize(uint64_t bytes){
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

// *******************
// * CLASS :: c_Mesh *
// *******************

template <class c_Type>
c_Mesh<c_Type>::c_Mesh(int n1, int n2, int n3):
	N1(n1), N2(n2), N3(n3), len(n1*n2*n3), p_data(len) { printf("Creating mesh of length %i...\n", len); }

template <class c_Type>
c_Mesh<c_Type>::c_Mesh(int n):
	N1(n), N2(n), N3(n),len(n*n*n), p_data(len) {}

template <class c_Type>
c_Mesh<c_Type>::c_Mesh(int n1, int n2, int n3, c_Type init_val):
	N1(n1), N2(n2), N3(n3), len(n1*n2*n3), p_data(len, init_val) {}

template <class c_Type>
c_Mesh<c_Type>::c_Mesh(int n, c_Type init_val):
	N1(n), N2(n), N3(n), len(n*n*n), p_data(len, init_val) {}

template <class c_Type>
c_Mesh<c_Type>& c_Mesh<c_Type>::operator+=(const c_Type& rhs)
{
	for (int i = 0; i < len; i++) this->p_data[i]+=rhs;
	return *this;
}

template <class c_Type>
c_Mesh<c_Type>& c_Mesh<c_Type>::operator*=(const c_Type& rhs)
{
	for (int i = 0; i < len; i++) this->p_data[i]*=rhs;
	return *this;
}

template <class c_Type>
c_Mesh<c_Type>& c_Mesh<c_Type>::operator/=(const c_Type& rhs)
{
	for (int i = 0; i < len; i++) this->p_data[i]/=rhs;
	return *this;
}

template <class c_Type>
int c_Mesh<c_Type>::get_kx(int index) const
{
	int tmp = get_x_c(index);
	if (tmp < N1/2.) return tmp;
	else if (tmp > N1/2.) return tmp - N1;
	else return 0;
}

template <class c_Type>
int c_Mesh<c_Type>::get_ky(int index) const
{
	int tmp = get_y_c(index);
	if (tmp < N2/2.) return tmp;
	else if (tmp > N2/2.) return tmp - N2;
	else return 0;
}

template <class c_Type>
int c_Mesh<c_Type>::get_kz(int index) const
{
	int tmp = get_z_c(index);
	if (tmp < N3/2. - 1) return tmp;
	else if (tmp > N3/2. -1) return tmp - N3;
	else return 0;
}

template <class c_Type>
void c_Mesh<c_Type>::get_k_vec(int index, int* k_vec) const
{
	k_vec[0] = get_kx(index);
	k_vec[1] = get_ky(index);
	k_vec[2] = get_kz(index);
}

template <class c_Type>
int c_Mesh<c_Type>::get_k_sq(int index) const
{
	int tmp, k_2 = 0;
	
	tmp = get_x_c(index);
	k_2+=((tmp <= N1/2.) ? pow(tmp,2) : pow(tmp - N1, 2.));
	tmp = get_y_c(index);
	k_2+=((tmp <= N2/2.) ? pow(tmp,2) : pow(tmp - N2, 2.));
	tmp = get_z_c(index);
	k_2+=((tmp <= N3/2.) ? pow(tmp,2) : pow(tmp - N3, 2.));
	
	return k_2;
}

// ************************
// * CLASS :: c_Sim_Param *
// ************************

int c_Sim_Param::init(int ac, char* av[])
{
	int err = handle_cmd_line(ac, av, this);
	if (err) {is_init = 0; return err;}
	else {
		is_init = 1;
		par_num = pow(mesh_num / Ng, 3);
		power.k2_G *= power.k2_G;
		power.eval_pwr();
		b_in = 1./(z_in + 1);
		b_out = 1./(z_out + 1);
		k_min = 2.*PI/box_size;
		k_max = 2.*PI*mesh_num/box_size;
		return err;
	}
}

void c_Sim_Param::print_info()
{
	if (is_init) 
	{
		printf("\n");
		printf("Num_par:\t%i\n", Ng);
		printf("Num_mesh:\t%i^3\n", mesh_num);
		printf("Box size:\t%i Mpc/h\n", box_size);
		printf("Starting redshift:\t%G\n", z_in);
		printf("The primordial power spectrum 'P(k)=A*k^ns' has amplitude A = %G and spectral index ns = %G.\n", power.A, power.ns);
		if (power.k2_G == 0) printf("Smoothing length was not set.\n");
		else printf("Smoothing wavenumber is %G h/Mpc.\n", sqrt(power.k2_G));
		printf("'viscozity' for adhesion approximation is %G px^2.\n", nu);
		printf("The program will try to use %i threads.\n", nt);
		cout << "Output will be written to folder '"<< out_dir << "'\n";
		printf("\n");
	}
	else printf("WARNING! Simulation parameters are not initialized!\n");
}

// *********************
// * CLASS :: c_Part_x *
// *********************

c_Part_x::c_Part_x (double x, double y, double z)
{
	position = { x, y, z };
}

void c_Part_x::get_per(int per)
{
	for (unsigned i = 0; i < 3; i++)
	{
		if (position[i] > per) position[i] = fmod(position[i], per);
		else if (position[i] < 0) position[i] = fmod(position[i], per) + per;
	}
}

void c_Part_x::drift(vector<double> velocity, double db)
{
	for(int i = 0; i < 3; i++) position[i]+=velocity[i]*db;
}

// *********************
// * CLASS :: c_Part_v *
// *********************

c_Part_v::c_Part_v(double x, double y, double z): c_Part_x(x, y, z), velocity(3) {}


c_Part_v::c_Part_v (double x, double y, double z, double vx, double vy, double vz):
	c_Part_x(x, y, z)
{
	velocity = { vx, vy, vz };
}

void c_Part_v::kick(vector<double> force, double db)
{
	for(int i = 0; i < 3; i++) velocity[i]+=force[i]*db;
}

// *****************************
// * CLASS :: c_Pow_Spec_Param *
// *****************************

static const double Omega_0 = 1;
static const double h = 0.67;

static double transfer_function_2(double k, void* parameters)
{
	if (k == 0) return 1.;
	double q = k / (Omega_0*h);
	double T_k =	log(1+2.34*q)/(2.34*q)*
					pow(1 + 3.89*q + pow(16.1*q, 2.) + pow(5.4*q, 3.) + pow(6.71*q, 4.)
					, -1./4.);
	return pow(T_k, 2.);
}

static double power_spectrum_T(double k, double* parameters)
{
	if (k == 0) return 0;
	double A = parameters[0];
	double ns = parameters[1];
	return A*pow(k, ns)*transfer_function_2(k, parameters);
}

static double power_spectrum_scale_free(double k, double* parameters)
{
	if (k == 0) return 0;
	double A = parameters[0];
	double ns = parameters[1];
	return A*pow(k, ns);
}

static double flat_power_spectrum(double k, double* parameters)
{
	double A = parameters[0];
	return A;
}

static double single_power_spectrum_T(double k, double* parameters)
{
	if ((k > 0.01) and (k < 0.04)) return power_spectrum_T(k, parameters);
	else return 0.;
}

static double power_spectrum(double k, double* parameters)
{
	double i_pwr_type = parameters[2];
	switch ((int)i_pwr_type)
	{
		case 0: return power_spectrum_T(k, parameters);
		case 1: return power_spectrum_scale_free(k, parameters);
		case 2: return flat_power_spectrum(k, parameters);
		case 3: return single_power_spectrum_T(k, parameters);
		default: return power_spectrum_T(k, parameters);
	}
}


static double power_spectrum_s8(double k, void* parameters)
{
	return	k*k/(2.*PI*PI)* // spherical factor
			power_spectrum(k, (double*)parameters)* // P(k)
			pow(3.*gsl_sf_bessel_j1(k*8)/(k*8), 2.); // window function R = 8 Mpc/h
}
	
void c_Pow_Spec_Param::eval_pwr()
{
	switch ((int)i_pwr_type){
		case 0: pwr_type = power_law_T; break;
		case 1: pwr_type = power_law; break;
		case 2: pwr_type = flat; break;
		case 3: pwr_type = single; break;
		default: pwr_type = power_law_T;
	}
}

void c_Pow_Spec_Param::norm_pwr()
{
	/* Normalize the power spectrum */
	printf("Computing normalization of the given power spectrum...\n");
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (1000);
	double result, error;
	
	double parameters[3] = {A, ns, i_pwr_type};
	gsl_function F;
	F.function = &power_spectrum_s8;
	F.params = parameters;
	printf("Before <gsl_integration_qagiu>\n");
	gsl_integration_qagiu (&F, 0, 0, 1e-7, 1000, w, &result, &error); 
	
	gsl_integration_workspace_free (w);
	A = s8*s8/result;
}

double c_Pow_Spec_Param::lin_pow_spec(double k) const
{
	double parameters[3] = {A, ns, i_pwr_type};
	return power_spectrum(k, parameters);
}

// ***************************
// * CLASS :: c_App_Var_base *
// ***************************

c_App_Var_base::c_App_Var_base(const c_Sim_Param &sim, string app_str):
	b(sim.b_in), b_out(sim.b_out), db(sim.b_in), z_suffix_const(app_str),
	app_field(3, c_Mesh<double>(sim.mesh_num, sim.mesh_num, sim.mesh_num + 2)),
	power_aux (sim.mesh_num, sim.mesh_num, sim.mesh_num + 2),
	pwr_spec_binned(sim.bin_num), pwr_spec_binned_0(sim.bin_num),
	pool(sim.nt), track(4, sim.par_num)
{
	// FFTW PREPARATION
	err = !fftw_init_threads();
	fftw_plan_with_nthreads(sim.nt);
	p_F = fftw_plan_dft_r2c_3d(sim.mesh_num, sim.mesh_num, sim.mesh_num, power_aux.data(),
		reinterpret_cast<fftw_complex*>(power_aux.data()), FFTW_ESTIMATE);
	p_B = fftw_plan_dft_c2r_3d(sim.mesh_num, sim.mesh_num, sim.mesh_num, reinterpret_cast<fftw_complex*>(power_aux.data()),
		power_aux.data(), FFTW_ESTIMATE);
}

c_App_Var_base::~c_App_Var_base()
{
	// FFTW CLEANUP
	fftw_destroy_plan(p_F);
	fftw_destroy_plan(p_B);
	fftw_cleanup_threads();
}

string c_App_Var_base::z_suffix()
{
	z_suffix_num.str("");
	z_suffix_num << fixed << setprecision(2) << z();
	return z_suffix_const + "z" + z_suffix_num.str();
}


void c_App_Var_base::upd_time()
{
	step++;
	if ((b_out - b) < db) db = b_out - b;
	else db = 0.01;
	b += db;
}

// **********************
// * CLASS :: c_App_Var *
// **********************

c_App_Var::c_App_Var(const c_Sim_Param &sim, string app_str):
	c_App_Var_base(sim, app_str)
{
	particles.reserve(sim.par_num);
	printf("Allocated %s of memory.\n", humanSize
	(
		sizeof(c_Part_v)*particles.capacity()+
		sizeof(double)*(app_field[0].p_data.capacity()*3+power_aux.p_data.capacity())
	));
}

static void set_unpert_pos_th(int i_min, int i_max, const c_Sim_Param &sim, vector<c_Part_x>* particles)
{
	double x, y, z;
	int par_per_dim = sim.mesh_num / sim.Ng;
	
	for(int i=i_min; i< i_max; i++)
	{
		x = i / (par_per_dim * par_per_dim);
		y = (i / par_per_dim) % par_per_dim;
		z = i % par_per_dim;
		
		particles->at(i) = c_Part_x(x, y, z);
	}
}

 void c_App_Var::set_unpert_pos(const c_Sim_Param &sim)
 {
	 if (!particles.empty())
	{
		printf("Reseting positions...\n");
		particles.clear();
		particles.reserve(sim.par_num);
	}
	printf("Setting initial positions of particles...\n");
	auto tmp_func = bind(set_unpert_pos_th, placeholders::_1, placeholders::_2, sim, &particles);
	pool.add_task(0, sim.par_num, tmp_func);
 }
 
// ************************
// * CLASS :: c_App_Var_v *
// ************************

c_App_Var_v::c_App_Var_v(const c_Sim_Param &sim, string app_str):
	c_App_Var_base(sim, app_str), particles(sim.par_num, c_Part_v(0,0,0,0,0,0))
{
	particles.reserve(sim.par_num);
	printf("Allocated %s of memory.\n", humanSize
	(
		sizeof(c_Part_v)*particles.capacity()+
		sizeof(double)*(app_field[0].p_data.capacity()*3+power_aux.p_data.capacity())
	));
}

static void set_unpert_pos_w_vel_th(int i_min, int i_max, const c_Sim_Param &sim, vector<c_Part_v>* particles, const std::vector< c_Mesh<double> > &vel_field)
{
	int x, y, z;
	double vx, vy, vz;
	int par_per_dim = sim.mesh_num / sim.Ng;
	
	for(int i=i_min; i< i_max; i++)
	{
		x = i / (par_per_dim * par_per_dim) * sim.Ng;
		y = (i / par_per_dim) % par_per_dim * sim.Ng;
		z = i % par_per_dim * sim.Ng;
		vx = vel_field[0].get_val(x, y, z);
		vy = vel_field[1].get_val(x, y, z);
		vz = vel_field[2].get_val(x, y, z);
		particles->at(i) = c_Part_v(x, y, z, vx, vy, vz);
	}
}

 void c_App_Var_v::set_unpert_pos_w_vel(const c_Sim_Param &sim, const std::vector< c_Mesh<double> > &vel_field)
 {
/*	 if (!particles.empty())
	{
		printf("Reseting positions and velocities...\n");
		particles.clear();
		particles.reserve(sim.par_num);
	} */
	printf("Setting initial positions and velocities of particles...\n");
	auto tmp_func = bind(set_unpert_pos_w_vel_th, placeholders::_1, placeholders::_2, sim, &particles, vel_field);
	pool.add_task(0, sim.par_num, tmp_func);
 }
 
// *******************
// * CLASS :: c_Pool *
// *******************

c_Pool::c_Pool(int num_thread):
	num_thread(num_thread), th_pool(num_thread)
{}

void c_Pool::add_task(int i_from, int i_to, function<void(int, int)> func)
{
	int i_min = i_from;
    int i_max = i_to; 
    int i_inc = (i_max - i_min) / num_thread;
    i_max = (i_max - i_min) % num_thread + i_inc;
	
    for (int i = 0; i < num_thread; i++)
    {
		auto func_one_thread = bind(func, i_min, i_max);		
        th_pool.schedule(func_one_thread);
        i_min = i_max;
        i_max += i_inc;
    }
    th_pool.wait();
}

// ***********************
// * CLASS :: c_Tracking *
// ***********************

c_Tracking::c_Tracking(int num_track_par, int par_num):
	num_track_par(num_track_par)
{
	printf("Initializing IDs of tracked particles...\n");
	par_ids.reserve(num_track_par*num_track_par);
	int par_dim = pow(par_num, 1/3.);
	int x, y, z;
	double s;
	y = par_dim / 2; // middle of the cube
	s = par_dim / (4.*(num_track_par+1.)); // quarter of the cube
	for (int i=1; i<=num_track_par;i++)
	{
		z = (int)(s*i);
		for (int j=1; j<=num_track_par;j++)
		{
			x = (int)(s*j);
			par_ids.push_back(x*par_dim*par_dim+y*par_dim+z);
		}
	}
}

// ************
// * ROUTINES *
// ************

template <typename t_Type>
t_Type mean(t_Type* p_data, int len)
{
	t_Type tmp = 0;
	for (int i = 0; i < len; i++) tmp += p_data[i];
	return tmp / len;
}

template <typename t_Type>
t_Type std_dev(t_Type* p_data, int len, t_Type t_mean)
{
	t_Type tmp = 0;
	for (int i = 0; i < len; i++) tmp += pow(p_data[i]-t_mean, 2.);
	return sqrt(tmp / len);
}

template <typename t_Type>
double norm(t_Type* p_data, int len)
{
	double tmp = 0;
	for (int i = 0; i < len; i++) tmp += pow(p_data[i], 2.);
	return sqrt(tmp);
}

template <typename t_Type>
double norm(t_Type* p_data)
{
	return norm<t_Type>(p_data, 3);
}

static void normalize_FFT_FORWARD_th(int i_min, int i_max, c_Mesh<double>* rho)
{
	for (int i = i_min; i < i_max; i++) (*rho)[i] /= sqrt(pow(rho->N1, 3));
//	for (int i = i_min; i < i_max; i++) (*rho)[i] /= pow(rho->N1, 3);
}

static void normalize_FFT_BACKWARD_th(int i_min, int i_max, c_Mesh<double>* rho)
{
	for (int i = i_min; i < i_max; i++) (*rho)[i] /= sqrt(pow(rho->N1, 3));
}

static void normalize_FFT_FORWARD(c_Mesh<double>* rho, c_Pool* pool)
{
	auto tmp_func = bind(normalize_FFT_FORWARD_th, placeholders::_1, placeholders::_2, rho);
	pool->add_task(0, rho->len, tmp_func);	
}

static void normalize_FFT_BACKWARD(c_Mesh<double>* rho, c_Pool* pool)
{
	auto tmp_func = bind(normalize_FFT_BACKWARD_th, placeholders::_1, placeholders::_2, rho);
	pool->add_task(0, rho->len, tmp_func);	
}

static double CIC_opt(int index, const c_Mesh<double> &pot)
{
	// optimalization for CIC
	int k_vec[3];
	double k_n[3];
	double U2, U_n, G_n, k2n;
	
	pot.get_k_vec(index, k_vec);
	G_n = 0;
	U2 = 1;
	for (int j = 0; j < 3; j++) U2 *= 1./3.*(1+2*cos(PI*k_vec[j]/pot.N1));
	for (int n1 = -N_MAX; n1 < N_MAX + 1; n1++)
	{
		k_n[0] = 2 * PI / pot.N1*k_vec[0] + 2 * PI*n1;
		for (int n2 = -N_MAX; n2 < N_MAX + 1; n2++)
		{
			k_n[1] = 2 * PI / pot.N1*k_vec[1] + 2 * PI*n2;
			for (int n3 = -N_MAX; n3 < N_MAX + 1; n3++)
			{
				k_n[2] = 2 * PI / pot.N1*k_vec[2] + 2 * PI*n3;
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
						G_n += 2 * PI / pot.N1*k_vec[j]* // D(k)
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

static void gen_gauss_white_noise_th(int i_min, int i_max, c_Mesh<double>* rho, unsigned long *slab_keys)
{
	unsigned long ikey, index;
	double rn1, rn2, rn;
	
	for(int i=i_min; i< i_max; i++)
	{
		ikey = slab_keys[i / rho->N2];
		for (int j=0; j < rho->N3-2; j++)
		{
			index = i % rho->N2 + j; // 2D index in slab
			GetRandomDoublesWhiteNoise(rn1, rn2, ikey, index);
			rn = rn1*rn1 + rn2*rn2;
			(*rho)[rho->index(i, j)] = rn2 * sqrt(-2.0*log(rn)/rn);
		}
	}
}

static void gen_gauss_white_noise(const c_Sim_Param &sim, c_Mesh<double>* rho, c_Pool* pool)
{	
	// Get keys for each slab in the x axis that this rank contains
	
	vector<unsigned long> slab_keys;
	slab_keys.resize(rho->N1);
	GetSlabKeys(slab_keys.data(), 0, rho->N1-2, sim.seed);
	
	// Manage multi-thread
	auto tmp_func = bind(gen_gauss_white_noise_th, placeholders::_1, placeholders::_2, rho, slab_keys.data());
	pool->add_task(0, rho->N1*rho->N2, tmp_func);
	
	#ifdef CORR
	double t_mean = mean(rho->data(), rho->len);
	double t_std_dev = std_dev(rho->data(), rho->len, t_mean);
	printf("\t[mean = %.12f, stdDev = %.12f]\t-->", t_mean, t_std_dev);
	(*rho)-=t_mean;
	(*rho)/=t_std_dev;
	#endif
	
	double tmp = mean(rho->data(), rho->len);
	printf("\t[mean = %.12f, stdDev = %.12f]\n", tmp, std_dev(rho->data(), rho->len, tmp));
}

void fftw_execute_dft_r2c(const fftw_plan &p_F, c_Mesh<double>* rho, c_Pool* pool)
{
	double *rho_x = rho->data();
	fftw_complex* rho_k = reinterpret_cast<fftw_complex*>(rho_x);	
	fftw_execute_dft_r2c(p_F, rho_x, rho_k);
	normalize_FFT_FORWARD(rho, pool);
}

void fftw_execute_dft_c2r(const fftw_plan &p_B, c_Mesh<double>* rho, c_Pool* pool)
{
	double *rho_x = rho->data();
	fftw_complex* rho_k = reinterpret_cast<fftw_complex*>(rho_x);	
	fftw_execute_dft_c2r(p_B, rho_k, rho_x);
	normalize_FFT_BACKWARD(rho, pool);
}

void fftw_execute_dft_r2c_triple(const fftw_plan &p_F, std::vector<c_Mesh<double>>* rho, c_Pool* pool)
{
	for (int i = 0; i < 3; i++) fftw_execute_dft_r2c(p_F, &(*rho)[i], pool);
}

void fftw_execute_dft_c2r_triple(const fftw_plan &p_B, std::vector<c_Mesh<double>>* rho, c_Pool* pool)
{
	for (int i = 0; i < 3; i++) fftw_execute_dft_c2r(p_B, rho->data() + i, pool);
}

static void gen_rho_w_pow_k_th(int i_min, int i_max, const c_Sim_Param &sim, c_Mesh<double>* rho_k)
{
	double k;
	for(int i=i_min; i < i_max;i++)
	{
		k = 2.*PI/sim.box_size*sqrt(rho_k->get_k_sq(i));
		(*rho_k)[2*i] *= sqrt(sim.power.lin_pow_spec(k));
		(*rho_k)[2*i+1] *= sqrt(sim.power.lin_pow_spec(k));
	}
}

static void gen_rho_w_pow_k(const c_Sim_Param &sim, c_Mesh<double>* rho_k, c_Pool* pool)
{
	auto tmp_func = bind(gen_rho_w_pow_k_th, placeholders::_1, placeholders::_2, sim, rho_k);
	pool->add_task(0, rho_k->len / 2, tmp_func);
	(*rho_k)[0]=(*rho_k)[1]=0.; // setting the zero mode to 0
}

static void pwr_spec_k_th(int i_min, int i_max, const c_Sim_Param &sim, const c_Mesh<double> &rho_k, c_Mesh<double>* power_aux)
{
	/* Computing the power spectrum P(k) */	
	/* Preserve values in rho_k */
	
	double w_k;
	int k_vec[3];
	for(int i=i_min; i < i_max;i++)
	{
		w_k = 1.;
		rho_k.get_k_vec(i, k_vec);
		for (int j = 0; j < 3; j++) if (k_vec[j] != 0) w_k *= pow(sin(PI*k_vec[j]/sim.mesh_num)/(PI*k_vec[j]/sim.mesh_num), sim.order + 1);
		(*power_aux)[2*i+1] = (rho_k[2*i]*rho_k[2*i] + rho_k[2*i+1]*rho_k[2*i+1])/(w_k*w_k);
		(*power_aux)[2*i] = 2.*PI/sim.box_size*norm(k_vec); // physical k
	}
}

static void gen_pot_k_th(int i_min, int i_max, c_Mesh<double>* rho_k)
{
	double k2;
	for(int i=i_min; i < i_max;i++){				
		k2 = rho_k->get_k_sq(i);
		if (k2 == 0){
			(*rho_k)[2*i] = 0;
			(*rho_k)[2*i+1] = 0;
		} else{
			(*rho_k)[2*i] /= -(k2*pow(2.*PI/rho_k->N1, 2.));
			(*rho_k)[2*i+1] /= -(k2*pow(2.*PI/rho_k->N1, 2.));
		}
	}
}

static void gen_displ_k_th(int i_min, int i_max, vector<c_Mesh<double>>* vel_field, c_Mesh<double>* pot_k)
{
	double opt;
	int k_vec[3];
	double potential_tmp[2];
	for(int i=i_min; i <i_max;i++)
	{
		potential_tmp[0] = (*pot_k)[2*i]; // prevent overwriting if vel_field[0] == pot_k
		potential_tmp[1] = (*pot_k)[2*i+1]; // prevent overwriting if vel_field[0] == pot_k
		opt = CIC_opt(i, (*pot_k));
		(*pot_k).get_k_vec(i, k_vec);		
		for(int j=0; j<3;j++)
		{
			(*vel_field)[j][2*i] = k_vec[j]*potential_tmp[1]*(2.*PI/(*pot_k).N1)*opt; // 2*PI/N comes from derivative WITH RESPECT to the mesh coordinates
			(*vel_field)[j][2*i+1] = -k_vec[j]*potential_tmp[0]*(2.*PI/(*pot_k).N1)*opt; // 2*PI/N comes from derivative WITH RESPECT to the mesh coordinates
		}

	}
}

static int get_per(int vec, int per)
{
	if (vec > per) vec = vec % per;
	else if (vec < 0) vec = vec % per + per;
	return vec;
}

double wgh_sch(const vector<double> &x, const vector<int> &y, int mesh_num, const int order)
{
	// The weighting scheme used to assign values to the mesh points or vice versa
	// Return value of assigment function on mesh point y from particle in x
	double dx, w = 1;
	int y_per;

	switch (order){
	case 0: {	// NGP: Nearest grid point
				for (int i = 0; i < 3; i++)
				{
					y_per = get_per(y[i], mesh_num);
					if ((int)x[i] != y_per) w *= 0;
				}
				return w;
	}
	case 1: {	// CIC: Cloud in cells
				for (int i = 0; i < 3; i++)
				{
					y_per = get_per(y[i], mesh_num);					
					dx = fmin(fmin(abs(x[i] - y_per), x[i] + mesh_num - y_per), y_per + mesh_num - x[i]);
					if (dx > 1) w *= 0;
					else w *= 1 - dx;
				}
				return w;
	}
	case 2: {	// TSC: Triangular shaped clouds
				for (int i = 0; i < 3; i++)
				{
					y_per = get_per(y[i], mesh_num);
					dx = fmin(fmin(abs(x[i] - y_per), x[i] + mesh_num - y_per), y_per + mesh_num - x[i]);
					if (dx > 1.5) w *= 0;
					else if (dx > 0.5) w *= (1.5 - dx)*(1.5 - dx) / 2.0;
					else w *= 3 / 4.0 - dx*dx;
				}
				return w;
	}
	}
	return 0;
}

static void reset_field_th(int i_min, int i_max, c_Mesh<double>* rho, const double value)
{
	fill(rho->p_data.begin()+i_min, rho->p_data.begin()+i_max, value);
}

static void reset_field(c_Mesh<double>* rho, const double value, c_Pool* pool)
{
	auto tmp_func = bind(reset_field_th, placeholders::_1, placeholders::_2, rho, value);
	pool->add_task(0, rho->len, tmp_func);
}

namespace c11 {
void work_dir_over(string out_dir){
	string wrk_dir;
	fs::path dir(out_dir.c_str());
	if(fs::create_directory(dir)){
        cout << "Directory Created: "<< out_dir << endl;
    }
	
	wrk_dir = out_dir + "par_cut/";
	dir = wrk_dir.c_str();
	if(fs::create_directory(dir)){
        cout << "Directory Created: "<< wrk_dir << endl;
    }
	
	wrk_dir = out_dir + "pwr_diff/";
	dir = wrk_dir.c_str();
	if(fs::create_directory(dir)){
        cout << "Directory Created: "<< wrk_dir << endl;
    }
	
	wrk_dir = out_dir + "pwr_spec/";
	dir = wrk_dir.c_str();
	if(fs::create_directory(dir)){
        cout << "Directory Created: "<< wrk_dir << endl;
    }
	
	wrk_dir = out_dir + "rho_map/";
	dir = wrk_dir.c_str();
	if(fs::create_directory(dir)){
        cout << "Directory Created: "<< wrk_dir << endl;
    }
}
}
// *****************
// * MAIN ROUTINES *
// *****************

void gen_rho_dist_k(const c_Sim_Param &sim, c_Mesh<double>* rho, const fftw_plan &p_F, c_Pool* pool)
	/**
	Generate density distributions \rho(k) in k-space.
	At first, a gaussian white noise (mean = 0, stdDev = 1) is generated,
	then it is convoluted with given power spectrum.
	**/
{
	// Gaussian white noies
	printf("Generating gaussian white noise...");
	gen_gauss_white_noise(sim, rho, pool);
	
	printf("Generating gaussian white noise in k-sapce...\n");
	fftw_execute_dft_r2c(p_F, rho, pool);
	
	/* Generate density distributions with given power spectrum */
	printf("Generating density distributions with given power spectrum...\n");
	gen_rho_w_pow_k(sim, rho, pool);
}

void pwr_spec_k(const c_Sim_Param &sim, const c_Mesh<double> &rho_k, c_Mesh<double>* power_aux, c_Pool* pool)
{
	printf("Computing the power spectrum P(k)...\n");
	auto tmp_func = bind(pwr_spec_k_th, placeholders::_1, placeholders::_2, sim, rho_k, power_aux);
	pool->add_task(0, rho_k.len / 2, tmp_func);
}

void pwr_spec(const c_Sim_Param &sim, c_Mesh<double>* rho, c_Mesh<double>* power_aux, const fftw_plan &p_F, c_Pool* pool)
{
	fftw_execute_dft_r2c(p_F, rho, pool);
	pwr_spec_k(sim, *rho, power_aux, pool);
}

void gen_pow_spec_binned(const c_Sim_Param &sim, const c_Mesh<double> &power_aux, vector<fftw_complex>* pwr_spec_binned)
{
	double log_bin = pow(sim.k_max / sim.k_min, 1./sim.bin_num);
	double k;
	int bin;
	printf("Computing binned power spectrum P(k)...\n");
	for (int j = 0; j < sim.bin_num; j++){
		(*pwr_spec_binned)[j][0] = 0.;
		(*pwr_spec_binned)[j][1] = 0.;
	}	
	for (int i = 0; i < power_aux.len / 2; i++){
		k = power_aux[2*i];
		if ((k <=sim.k_max) && (k>=sim.k_min)){
			bin = (int)(log(k/sim.k_min)/log(log_bin));
			(*pwr_spec_binned)[bin][1] += power_aux[2*i+1]; // P(k)
			(*pwr_spec_binned)[bin][0]++;
		}
	}
		
	k = sim.k_min*sqrt(log_bin);
	for (int j = 0; j < sim.bin_num; j++){
		if ((*pwr_spec_binned)[j][0]) (*pwr_spec_binned)[j][1] /= (*pwr_spec_binned)[j][0];
		(*pwr_spec_binned)[j][0] = k;
		k *= log_bin;
	}
}

void gen_pot_k(c_Mesh<double>* rho_k, c_Pool* pool)
{
	printf("Computing potential in k-space...\n");
	auto tmp_func = bind(gen_pot_k_th, placeholders::_1, placeholders::_2, rho_k);
	pool->add_task(0, rho_k->len / 2, tmp_func);
}

void gen_displ_k(vector<c_Mesh<double>>* vel_field, c_Mesh<double>* pot_k, c_Pool* pool)
{
	printf("Computing displacement in k-space...\n");
	auto tmp_func = bind(gen_displ_k_th, placeholders::_1, placeholders::_2, vel_field, pot_k);
	pool->add_task(0, pot_k->len / 2, tmp_func);
}

void get_rho_from_par(const vector<c_Part_v> &particles, c_Mesh<double>* rho, const int order, c_Pool* pool)
{
	printf("Computing the density field from particle positions...\n");
	double m = pow((double)rho->N1, 3.)/particles.size();
	reset_field(rho, -1., pool);
	for (const auto& one_par : particles) rho->assign_to(one_par.position, m, order);
}

void print_pow_spec(const vector<fftw_complex> &pwr_spec_binned, string out_dir, string suffix)
{
	out_dir += "pwr_spec/";
	string file_name = out_dir + "pwr_spec" + suffix + ".dat";
	FILE* pFile;
	pFile = fopen(file_name.c_str(), "w");
	if (pFile == NULL)
	{
		printf("Error while opening %s\n", file_name.c_str());
		perror("Error");
		return;
	}
	
	cout << "Writing power spectrum into file " << file_name << endl;
	fprintf (pFile, "# This file contains power spectrum P(k) in units [(Mpc/h)^3] depending on wavenumber k in units [h/Mpc].\n");
	fprintf (pFile, "# k [h/Mpc]\tP(k) [(Mpc/h)^3]\n");
	
	for (unsigned j = 0; j < pwr_spec_binned.size(); j++){
		if (pwr_spec_binned[j][1]) fprintf (pFile, "%f\t%f\n",  pwr_spec_binned[j][0], pwr_spec_binned[j][1]);
	}

	fclose (pFile);
}

void print_par_pos_cut_small(const vector<c_Part_v> &particles, int mesh_num, int L, string out_dir, string suffix)
{
	out_dir += "par_cut/";
	string file_name = out_dir + "par_cut" + suffix + ".dat";
	FILE* pFile;
	pFile = fopen(file_name.c_str(), "w");
	if (pFile == NULL)
	{
		printf("Error while opening %s\n", file_name.c_str());
		perror("Error");
		return;
	}
	
	cout << "Writing small cut through the box of particles into file " << out_dir + "par_cut" + suffix + ".dat\n";
	fprintf (pFile, "# This file contains positions of particles in units [Mpc/h].\n");
	double x, y, z, dx;
	for(unsigned i=0; i < particles.size(); i++)
	{
		x = particles[i].position[0];
		y = particles[i].position[1];
		z = particles[i].position[2];			
		dx = abs(y - mesh_num/2.);
		if ((dx < 0.5) && (x < mesh_num/4.) && (z < mesh_num/4.))
		{
			// cut (L/4 x L/4 x 0.5)
			fprintf (pFile, "%f\t%f\n", x/mesh_num*L , z/mesh_num*L);
		}
	}
	fclose (pFile);
}