/**
 * @file:	core.cpp
 * @brief:	class definitions
 */
 
#include "stdafx.h"
#include <fftw3.h>
#include "core.h"
#include "core_cmd.h"
#include "core_out.h"
#include "core_app.h"
// #include "core_mesh.h"

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
 * @brief:	creates a mesh of N*N*(N+2) cells
 */

Mesh::Mesh(int n): Mesh_base(n, n, n+2), N(n) {}

Mesh::Mesh(const Mesh& that): Mesh_base(that), N(that.N) {}

double& Mesh::operator()(Vec_3D<int> pos)
{
	get_per(pos, N);
	return data[pos.x*N2*N3+pos.y*N3+pos.z]; 
}

const double & Mesh::operator()(Vec_3D<int> pos) const
{
	get_per(pos, N);
	return data[pos.x*N2*N3+pos.y*N3+pos.z];
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

void swap(Mesh& first, Mesh& second)
{
	std::swap(first.length, second.length);
	std::swap(first.N, second.N);
	std::swap(first.data, second.data);
}

Mesh& Mesh::operator=(const Mesh& other)
{
//	printf("Copy assignemnt %p\n", this);
	Mesh temp(other);
	swap(*this, temp);
    return *this;
}

/**
 * @class:	Tracking
 * @brief:	class storing info about tracked particles
 */

Tracking::Tracking(int sqr_num_track_par, int par_num_per_dim):
	sqr_num_track_par(sqr_num_track_par), num_track_par(sqr_num_track_par*sqr_num_track_par)
{
	printf("Initializing IDs of tracked particles...\n");
	par_ids.reserve(num_track_par);
	int x, y, z;
	double s;
	y = par_num_per_dim / 2; // middle of the cube
	s = par_num_per_dim / (4.*(sqr_num_track_par+1.)); // quarter of the cube
	for (int i=1; i<=sqr_num_track_par;i++)
	{
		z = (int)(s*i);
		for (int j=1; j<=sqr_num_track_par;j++)
		{
			x = (int)(s*j);
			par_ids.push_back(x*par_num_per_dim*par_num_per_dim+y*par_num_per_dim+z);
		}
	}
}

void Tracking::update_track_par(Particle_x* particles)
{
	vector<Particle_x> par_pos_step;
	par_pos_step.reserve(num_track_par);
	for (int i=0; i<num_track_par; i++){
		par_pos_step.push_back(particles[par_ids[i]]);
	}
	par_pos.push_back(par_pos_step);
}

void Tracking::update_track_par(Particle_v* particles)
{
	vector<Particle_x> par_pos_step;
	par_pos_step.reserve(num_track_par);
	for (int i=0; i<num_track_par; i++){
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
		// if(nt == 0) nt = omp_get_num_procs();
        if(nt == 0) nt = omp_get_max_threads();
		else omp_set_num_threads(nt);
		par_num = pow(mesh_num / Ng, 3);
		power.k2_G *= power.k2_G;
		b_in = 1./(z_in + 1);
		b_out = 1./(z_out + 1);
		k_min = 2.*PI/box_size;
		k_max = 2.*PI*mesh_num/box_size;
		
//		rs = 1.0;
		a = rs / 0.735;
		M = (int)(mesh_num / rs);
		Hc = double(mesh_num) / M;
		
		return err;
	}
}

void Sim_Param::print_info(string out) const
{
    string file_name;
    FILE* pFile;

    if (out == "") pFile = stdout;
    else
    {
        file_name = out + "sim_param.log";
	    pFile = fopen(file_name.c_str(), "w");
    }
    
	if (is_init) 
	{
		fprintf(pFile, "\n*********************\n");
		fprintf(pFile, "SIMULATION PARAMETERS\n");
		fprintf(pFile, "*********************\n");
		fprintf(pFile, "Ng:\t\t%i\n", Ng);
		fprintf(pFile, "Num_par:\t%G^3\n", pow(par_num, 1/3.));
		fprintf(pFile, "Num_mesh:\t%i^3\n", mesh_num);
		fprintf(pFile, "Box size:\t%i Mpc/h\n", box_size);
		fprintf(pFile, "Redshift:\t%G--->%G\n", z_in, z_out);
		fprintf(pFile, "Pk:\t\t[sigma_8 = %G, As = %G, ns = %G, k_smooth = %G, pwr_type = %i]\n", 
			power.s8, power.A, power.ns, sqrt(power.k2_G), power.pwr_type);
		fprintf(pFile, "AA:\t\t[nu = %G px^2]\n", nu);
		fprintf(pFile, "LL:\t\t[rs = %G, a = %G, M = %i, Hc = %G]\n", rs, a, M, Hc);
		fprintf(pFile, "num_thread:\t%i\n", nt);
        fprintf(pFile,  "Output:\t\t'%s'\n", out_dir.c_str());
	}
	else fprintf(pFile, "WARNING! Simulation parameters are not initialized!\n");

    if (out != "") fclose (pFile);
}


void Sim_Param::print_info() const
{
	Sim_Param::print_info("");
}

/**
 * @class:	App_Var_base
 * @brief:	class containing variables for approximations
 */
 
App_Var_base::App_Var_base(const Sim_Param &sim, string app_str):
	err(0), step(0), print_every(1),
	b(sim.b_in), b_out(sim.b_out), db(sim.db), z_suffix_const(app_str),
	app_field(3, Mesh(sim.mesh_num)),
	power_aux (sim.mesh_num),
	pwr_spec_binned(sim.bin_num), pwr_spec_binned_0(sim.bin_num),
	track(4, sim.mesh_num/sim.Ng),
	dens_binned(20)
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

void App_Var_base::print_x(const Sim_Param &sim, string out_dir_app, Particle_x* particles)
{
	/* Printing positions */
	print_par_pos_cut_small(particles, sim, out_dir_app, z_suffix());
	track.update_track_par(particles);
	print_track_par(track, sim, out_dir_app, z_suffix());

	/* Printing density */
	get_rho_from_par(particles, &power_aux, sim);
	gen_dens_binned(power_aux, dens_binned, sim);
	print_rho_map(power_aux, sim, out_dir_app, z_suffix());
	print_dens_bin(dens_binned, sim.mesh_num, out_dir_app, z_suffix());

	/* Printing power spectrum */
	fftw_execute_dft_r2c(p_F, power_aux);
	pwr_spec_k(sim, power_aux, &power_aux);
	gen_pow_spec_binned(sim, power_aux, &pwr_spec_binned);
	print_pow_spec(pwr_spec_binned, out_dir_app, z_suffix());
	print_pow_spec_diff(pwr_spec_binned, pwr_spec_binned_0, b, out_dir_app, z_suffix());

	upd_supp();
}

void App_Var_base::print_v(const Sim_Param &sim, string out_dir_app, Particle_v* particles)
{
	/* Printing positions */
	print_par_pos_cut_small(particles, sim, out_dir_app, z_suffix());
	track.update_track_par(particles);
	print_track_par(track, sim, out_dir_app, z_suffix());

	/* Printing density */
	get_rho_from_par(particles, &power_aux, sim);
	gen_dens_binned(power_aux, dens_binned, sim);
	print_rho_map(power_aux, sim, out_dir_app, z_suffix());
	print_dens_bin(dens_binned, sim.mesh_num, out_dir_app, z_suffix());

	/* Printing power spectrum */
	fftw_execute_dft_r2c(p_F, power_aux);
	pwr_spec_k(sim, power_aux, &power_aux);
	gen_pow_spec_binned(sim, power_aux, &pwr_spec_binned);
	print_pow_spec(pwr_spec_binned, out_dir_app, z_suffix());
	print_pow_spec_diff(pwr_spec_binned, pwr_spec_binned_0, b, out_dir_app, z_suffix());

	upd_supp();
}

void App_Var_base::upd_time()
{
	step++;
	if ((b_out - b) < db) db = b_out - b;
	b += db;
}

void App_Var_base::upd_supp()
{
	double P_k, P_ZA, supp_tmp = 0;
	int i = 0, j = 0;
	while (i < 10){
		P_k = pwr_spec_binned[j][1];
		P_ZA = pwr_spec_binned_0[j][1] * pow(b, 2.);
		if((P_ZA) && (P_k))
		{
			supp_tmp += (P_k-P_ZA)/P_ZA;
			i++;
		}
		j++;
	}
	supp.push_back(double_2{ {b, supp_tmp / i} });
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

/**
 * @class:	App_Var_AA
 * @brief:	class containing variables for adhesion approximation
 */
 
 App_Var_AA::App_Var_AA(const Sim_Param &sim, string app_str):
	App_Var_base(sim, app_str), expotential (sim.mesh_num)
{
	particles = new Particle_x[sim.par_num];
	printf("Allocated %s of memory.\n", humanSize
	(
		sizeof(Particle_x)*sim.par_num
		+sizeof(double)*(app_field[0].length*3+power_aux.length +expotential.length)
	));
}

App_Var_AA::~App_Var_AA()
{
	delete[] particles;
}

/**
 * @class:	App_Var_FP_mod
 * @brief:	class containing variables for modified Frozen-potential approximation
 */
 
 App_Var_FP_mod::App_Var_FP_mod(const Sim_Param &sim, string app_str):
	App_Var_base(sim, app_str), linked_list (sim.par_num, sim.M, sim.Hc)
{
	particles = new Particle_v[sim.par_num];
	printf("Allocated %s of memory.\n", humanSize
	(
		sizeof(Particle_v)*sim.par_num
		+sizeof(double)*(app_field[0].length*3+power_aux.length)
		+sizeof(int)*(linked_list.HOC.length+linked_list.par_num)
	));
}

App_Var_FP_mod::~App_Var_FP_mod()
{
	delete[] particles;
}

/**
 * @class LinkedList
 * @brief class handling linked lists
 */


LinkedList::LinkedList(int par_num, int m, double hc):
	par_num(par_num), Hc(hc), LL(par_num), HOC(m, m, m) {}
	
void LinkedList::get_linked_list(Particle_v* particles)
{
	HOC.assign(-1);
	for (int i = 0; i < par_num; i++)
	{
		LL[i] = HOC(Vec_3D<int>(particles[i].position/Hc));
		HOC(Vec_3D<int>(particles[i].position/Hc)) = i;
	}
}