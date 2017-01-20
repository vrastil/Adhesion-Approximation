
#include "stdafx.h"
#include <fftw3.h>
#include "core.h"
#include "core_cmd.h"
#include "core_mesh.h"

using namespace std;
const double PI = acos(-1.);

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

/**
 * @class:	Mesh
 * @brief:	class handling basic mesh functions, the most important are creating and destroing the underlying data structure
 *			creates a mesh of N*N*(N+2) cells
 */

Mesh::Mesh(int n):N(n), length(n*n*(n+2))
{
	data = new double[length];
//	printf("Normal ctor %p\n", this); 
}

Mesh::Mesh(const Mesh& that): N(that.N), length(that.length)
{
	data = new double[length];
	
	#pragma omp parallel for
	for (int i = 0; i < length; i++) data[i] = that.data[i];
//	printf("Copy ctor %p\n", this);
}

void swap(Mesh& first, Mesh& second)
{
	std::swap(first.length, second.length);
	std::swap(first.N, second.N);
	std::swap(first.data, second.data);
}

Mesh& Mesh::operator=(Mesh& other)
{
	swap(*this, other);
    return *this;
//	printf("Copy assignemnt %p\n", this);
}

Mesh::~Mesh()
{
	delete[] data;
//	printf("dtor %p\n", this);
}

double& Mesh::operator()(Vec_3D<int> pos)
{
	get_per(pos, N);
	return data[pos.x*N*(N+2)+pos.y*(N+2)+pos.z]; 
}

const double& Mesh::operator()(Vec_3D<int> pos) const
{
	get_per(pos, N);
	return data[pos.x*N*(N+2)+pos.y*(N+2)+pos.z];
}

Mesh& Mesh::operator+=(const double& rhs)
{
	#pragma omp parallel for
		for (int i = 0; i < length; i++) this->data[i]+=rhs;
		
	return *this;
}

Mesh& Mesh::operator*=(const double& rhs)
{
	#pragma omp parallel for
		for (int i = 0; i < length; i++) this->data[i]*=rhs;
		
	return *this;
}

Mesh& Mesh::operator/=(const double& rhs)
{
	#pragma omp parallel for
		for (int i = 0; i < length; i++) this->data[i]/=rhs;
		
	return *this;
}

/**
 * @class:	Tracking
 * @brief:	class storing info about tracked particles
 */

Tracking::Tracking(int num_track_par, int par_num):
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

void Tracking::update_track_par(Particle_x* particles)
{
	vector<Particle_x> par_pos_step;
	par_pos_step.reserve(num_track_par*num_track_par);
	for (int i=0; i<num_track_par*num_track_par; i++){
		par_pos_step.push_back(particles[par_ids[i]]);
	}
	par_pos.push_back(par_pos_step);
}

/**
 * @class:	Sim_Param
 * @brief:	class storing simulation parameters
 */
 
int Sim_Param::init(int ac, char* av[])
{
	int err = handle_cmd_line(ac, av, this);
	if (err) {is_init = 0; return err;}
	else {
		is_init = 1;
		if(nt == 0) nt = omp_get_num_threads();
		else omp_set_num_threads(nt);
		par_num = pow(mesh_num / Ng, 3);
		power.k2_G *= power.k2_G;
		b_in = 1./(z_in + 1);
		b_out = 1./(z_out + 1);
		k_min = 2.*PI/box_size;
		k_max = 2.*PI*mesh_num/box_size;
		return err;
	}
}

void Sim_Param::print_info()
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

/**
 * @class:	App_Var_base
 * @brief:	class containing variables for approximations
 */
 
App_Var_base::App_Var_base(const Sim_Param &sim, string app_str):
	b(sim.b_in), b_out(sim.b_out), db(sim.b_in), z_suffix_const(app_str),
	app_field(3, Mesh(sim.mesh_num)),
	power_aux (sim.mesh_num),
	pwr_spec_binned(sim.bin_num), pwr_spec_binned_0(sim.bin_num),
	track(4, sim.par_num)
{
	// FFTW PREPARATION
	err = !fftw_init_threads();
	if (err){
		printf("Errors during multi-thread initialization!\n");
		throw err;
	}
	fftw_plan_with_nthreads(sim.nt);
	p_F = fftw_plan_dft_r2c_3d(sim.mesh_num, sim.mesh_num, sim.mesh_num, power_aux.real(),
		power_aux.complex(), FFTW_ESTIMATE);
	p_B = fftw_plan_dft_c2r_3d(sim.mesh_num, sim.mesh_num, sim.mesh_num, power_aux.complex(),
		power_aux.real(), FFTW_ESTIMATE);
}

App_Var_base::~App_Var_base()
{
	// FFTW CLEANUP
	fftw_destroy_plan(p_F);
	fftw_destroy_plan(p_B);
	fftw_cleanup_threads();
}

string App_Var_base::z_suffix()
{
	z_suffix_num.str("");
	z_suffix_num << fixed << setprecision(2) << z();
	return z_suffix_const + "z" + z_suffix_num.str();
}

void App_Var_base::upd_time()
{
	step++;
	if ((b_out - b) < db) db = b_out - b;
	else db = 0.01;
	b += db;
}

/**
 * @class:	App_Var
 * @brief:	class containing variables for approximations with particle positions only
 */
 
App_Var::App_Var(const Sim_Param &sim, string app_str):
	App_Var_base(sim, app_str)
{
	particles = new Particle_x[sim.par_num];
	printf("Allocated %s of memory.\n", humanSize
	(
		sizeof(Particle_x)*sim.par_num
		+sizeof(double)*(app_field[0].length*3+power_aux.length)
	));
}

App_Var::~App_Var()
{
	delete[] particles;
}

 /**
 * @class:	App_Var_v
 * @brief:	class containing variables for approximations with particle velocities
 */
 
App_Var_v::App_Var_v(const Sim_Param &sim, string app_str):
	App_Var_base(sim, app_str)
{
	particles = new Particle_v[sim.par_num];
	printf("Allocated %s of memory.\n", humanSize
	(
		sizeof(Particle_v)*sim.par_num
		+sizeof(double)*(app_field[0].length*3+power_aux.length)
	));
}

App_Var_v::~App_Var_v()
{
	delete[] particles;
}