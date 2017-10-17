/**
 * @file:	core.cpp
 * @brief:	class definitions
 */
 
#include "stdafx.h"
#include "core.h"
#include "core_cmd.h"
#include "core_out.h"
#include "core_app.h"
#include "core_mesh.h"
#include "core_power.h"

using namespace std;
using json = nlohmann::json;

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
 * @class:	Vec_3D<T>
 * @brief:	class handling basic 3D-vector functions
 */

 template <typename T>
 double Vec_3D<T>::norm() const
 {
     T tmp(0);
     for (const T val : vec)
     {
         tmp += val*val;
     }
     return sqrt(tmp);
 }
 
 template <typename T>
 Vec_3D<T>& Vec_3D<T>::operator+=(const Vec_3D<T>& rhs)
 {
     for(unsigned i = 0; i < 3; ++i)
     {
         vec[i] += rhs[i];
     }
     return *this;
 }
 
 template <typename T>
 Vec_3D<T> operator+(Vec_3D<T> lhs, const Vec_3D<T>& rhs)
 {
     lhs += rhs;
     return lhs;
 }
 
 template <typename T>
 Vec_3D<T>& Vec_3D<T>::operator-=(const Vec_3D<T>& rhs)
 {
     for(unsigned i = 0; i < 3; ++i)
     {
         vec[i] -= rhs[i];
     }
     return *this;
 }
 
 template <typename T>
 Vec_3D<T> operator-(Vec_3D<T> lhs, const Vec_3D<T>& rhs)
 {
     lhs -= rhs;
     return lhs;
 }
 
 template <typename T>
 Vec_3D<T>& Vec_3D<T>::operator+=(T rhs)
 {
     for(T& val : vec)
     {
         val += rhs;
     }
     return *this;
 }
 
 template <typename T>
 Vec_3D<T>& Vec_3D<T>::operator-=(T rhs)
 {
     for(T& val : vec)
     {
         val -= rhs;
     }
     return *this;
 }
 
 template <typename T>
 Vec_3D<T>& Vec_3D<T>::operator*=(T rhs)
 {
     for(T& val : vec)
     {
         val *= rhs;
     }
     return *this;
 }
 
 template <typename T>
 Vec_3D<T> operator*(Vec_3D<T> lhs, T rhs)
 {
     lhs *= rhs;
     return lhs;
 }
 
 template <typename T>
 Vec_3D<T> operator*(T lhs, Vec_3D<T> rhs)
 {
     rhs *= lhs;
     return rhs;
 }
 
 template <typename T>
 Vec_3D<T> operator+(Vec_3D<T> lhs, T rhs)
 {
     lhs += rhs;
     return lhs;
 }
 
 template <typename T>
 Vec_3D<T> operator+(T lhs, Vec_3D<T> rhs)
 {
     rhs += lhs;
     return rhs;
 }
 
 template <typename T>
 Vec_3D<T> operator-(Vec_3D<T> lhs, T rhs)
 {
     lhs -= rhs;
     return lhs;
 }
 
 template <typename T>
 Vec_3D<T> operator-(T lhs, Vec_3D<T> rhs)
 {
     rhs -= lhs;
     return rhs;
 }
 
 template <typename T>
 Vec_3D<T>& Vec_3D<T>::operator/=(T rhs)
 {
     for(T& val : vec)
     {
         val /= rhs;
     }
     return *this;
 }
 
 template <typename T>
 Vec_3D<T> operator/(Vec_3D<T> lhs, T rhs)
 {
     lhs /= rhs;
     return lhs;
 }
 
 template <typename T>
 template<class U>
 Vec_3D<T>::operator Vec_3D<U>() const
 {
     Vec_3D<U> lhs;
     for(unsigned i = 0; i < 3; ++i)
     {
         lhs[i] = static_cast<U>((*this)[i]);
     }
     return lhs;
 }
 
 /**
  * @class:	Mesh_base<T>
  * @brief:	class handling basic mesh functions, the most important are creating and destroing the underlying data structure
  *			creates a mesh of N1*N2*N3 cells
  */
 
 template <typename T>
 Mesh_base<T>::Mesh_base(unsigned n1, unsigned n2, unsigned n3):N1(n1), N2(n2), N3(n3), length(n1*n2*n3)
 {
     data = new T[length];
 //	printf("Normal base ctor %p, N1 = %i, N2 = %i, N3 = %i\n", this, N1, N2, N3); 
 }
 
 template <typename T>
 Mesh_base<T>::Mesh_base(const Mesh_base<T>& that): N1(that.N1), N2(that.N2), N3(that.N3), length(that.length)
 {
     data = new T[length];
     
     #pragma omp parallel for
     for (unsigned i = 0; i < length; i++) data[i] = that.data[i];
 //	printf("Copy base ctor %p\n", this);
 }
 
 template <typename T>
 void swap(Mesh_base<T>& first, Mesh_base<T>& second)
 {
     std::swap(first.length, second.length);
     std::swap(first.N1, second.N1);
     std::swap(first.N2, second.N2);
     std::swap(first.N3, second.N3);
     std::swap(first.data, second.data);
 }
 
 template <typename T>
 Mesh_base<T>& Mesh_base<T>::operator=(const Mesh_base<T>& other)
 {
 //	printf("Copy base assignemnt %p\n", this);
     Mesh_base<T> temp(other);
     swap(*this, temp);
     return *this;
 }
 
 template <typename T>
 Mesh_base<T>::~Mesh_base<T>()
 {
     delete[] data;
 //	printf("dtor base %p\n", this);
 }
 
 template <typename T>
 T& Mesh_base<T>::operator()(Vec_3D<int> pos)
 {
     get_per(pos, N1, N2, N3);
     return data[pos[0]*N2*N3+pos[1]*N3+pos[2]]; 
 }
 
 template <typename T>
 const T& Mesh_base<T>::operator()(Vec_3D<int> pos) const
 {
     get_per(pos, N1, N2, N3);
     return data[pos[0]*N2*N3+pos[1]*N3+pos[2]];
 }
 
 template <typename T>
 Mesh_base<T>& Mesh_base<T>::operator+=(const T& rhs)
 {
     #pragma omp parallel for
         for (unsigned i = 0; i < length; i++) this->data[i]+=rhs;
         
     return *this;
 }
 
 template <typename T>
 Mesh_base<T>& Mesh_base<T>::operator*=(const T& rhs)
 {
     #pragma omp parallel for
         for (unsigned i = 0; i < length; i++) this->data[i]*=rhs;
         
     return *this;
 }
 
 template <typename T>
 Mesh_base<T>& Mesh_base<T>::operator/=(const T& rhs)
 {
     #pragma omp parallel for
         for (unsigned i = 0; i < length; i++) this->data[i]/=rhs;
         
     return *this;
 }
 
 template <typename T>
 void Mesh_base<T>::assign(T val)
 {
     #pragma omp parallel for
     for (unsigned i = 0; i < length; i++) this->data[i]=val;
 }

/**
 * @class:	Mesh
 * @brief:	creates a mesh of N*N*(N+2) cells
 */

Mesh::Mesh(unsigned n): Mesh_base(n, n, n+2), N(n) {}

Mesh::Mesh(const Mesh& that): Mesh_base(that), N(that.N) {}

void Mesh::reset_part(bool part)
{
    /* nullify real (part = 0) or complex (part = 1) part of a field */
    #pragma omp parallel for
    for (unsigned i = part; i < this->length; i+=2){
        data[i] = 0;
    }
}

double& Mesh::operator()(Vec_3D<int> pos)
{
	get_per(pos, N);
	return data[pos[0]*N2*N3+pos[1]*N3+pos[2]]; 
}

const double & Mesh::operator()(Vec_3D<int> pos) const
{
	get_per(pos, N);
	return data[pos[0]*N2*N3+pos[1]*N3+pos[2]];
}

Mesh& Mesh::operator+=(const double& rhs)
{
	#pragma omp parallel for
		for (unsigned i = 0; i < length; i++) this->data[i]+=rhs;
		
	return *this;
}

Mesh& Mesh::operator*=(const double& rhs)
{
	#pragma omp parallel for
		for (unsigned i = 0; i < length; i++) this->data[i]*=rhs;
		
	return *this;
}

Mesh& Mesh::operator/=(const double& rhs)
{
	#pragma omp parallel for
		for (unsigned i = 0; i < length; i++) this->data[i]/=rhs;
		
	return *this;
}

void swap(Mesh& first, Mesh& second)
{
    std::swap(first.N, second.N);
    swap<Mesh_base<double>>(first, second);
}

Mesh& Mesh::operator=(const Mesh& other)
{
//	printf("Copy assignemnt %p\n", this);
	Mesh temp(other);
	swap(*this, temp);
    return *this;
}

/**
 * @class:	Pow_Spec_Param
 * @brief:	class storing parameters for power spectrum
 */

void Pow_Spec_Param::init()
{
    is_init = true;
    // CCL VARIABLES
    switch (pwr_type)
    {
        case ccl_EH: config.transfer_function_method = ccl_eisenstein_hu; break;
        case ccl_BBKS: config.transfer_function_method = ccl_bbks; break;
        default: config.transfer_function_method = ccl_bbks; break;
    }
    config.matter_power_spectrum_method = ccl_linear;
    config.mass_function_method = ccl_tinker;
    int status = 0;
    params = ccl_parameters_create_flat_lcdm(Omega_c, Omega_b, h, sigma8, ns, &status);
    cosmo = ccl_cosmology_create(params, config);
}

Pow_Spec_Param::~Pow_Spec_Param()
{
    if(is_init){
        ccl_cosmology_free(cosmo);
    }
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

 template <class T>
 void Tracking::update_track_par(T* particles)
 {
     std::vector<Particle_x> par_pos_step;
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
 
Sim_Param::Sim_Param(int ac, char* av[])
{
	if (handle_cmd_line(ac, av, this)) is_init = 0;
	else {
		is_init = 1;
		/* MULTITHREADING */
        if(nt == 0) nt = omp_get_max_threads();
        else omp_set_num_threads(nt);

        /* SIMULATION BOX*/
        Ng_pwr = Ng*mesh_num_pwr/mesh_num;
        par_num = pow(mesh_num / Ng, 3);
        a = rs / 0.735;
        M = (int)(mesh_num / rs);
        Hc = double(mesh_num) / M;

        /* TIME */
        b_in = 1./(z_in + 1);
		b_out = 1./(z_out + 1);
        
        /* POWER SPECTRUM */
        power.pwr_type = static_cast<e_power_spec>(power.pwr_type_i);
        power.k2_G *= power.k2_G;
        power.h = power.H0/100;
        power.init();
        norm_pwr(&power);

        /* RANGE : k */
        k_par.k_print.lower = 2*PI/box_size;
        k_par.k_print.upper = 2*PI/box_size*mesh_num_pwr;
        k_par.k_interp.lower = get_max_Pk(this);
        k_par.k_interp.upper = PI*pow(par_num, 1/3.) / (2*box_size); // k_Nyq / 2 (for safety)

        k_par.nyquist["analysis"] = PI*mesh_num_pwr / box_size;
        k_par.nyquist["potential"] = PI*mesh_num / box_size;
        k_par.nyquist["particle"] = PI*pow(par_num, 1/3.) / box_size;

        /* CORRELATION FUNCTION */
        x_corr.lower = 1;
        x_corr.upper = 200;
        corr_int = corr_int_type::QAWF;
    }
}

void Sim_Param::print_info(string out, string app) const
{
    if (is_init) 
	{
        if (out == "")
        {
            printf("\n*********************\n");
            printf("SIMULATION PARAMETERS\n");
            printf("*********************\n");
            printf("Ng:\t\t%i\n", Ng);
            printf("Num_par:\t%G^3\n", pow(par_num, 1/3.));
            printf("Num_mesh:\t%i^3\n", mesh_num);
            printf("Num_mesh_pwr:\t%i^3\n", mesh_num_pwr);
            printf("Box size:\t%.0f Mpc/h\n", box_size);
            printf("Redshift:\t%G--->%G\n", z_in, z_out);
            printf("Pk:\t\t[sigma_8 = %G, As = %G, ns = %G, k_smooth = %G, pwr_type = %i]\n", 
                power.sigma8, power.A, power.ns, sqrt(power.k2_G), power.pwr_type);
            printf("AA:\t\t[nu = %G px^2]\n", nu);
            printf("LL:\t\t[rs = %G, a = %G, M = %i, Hc = %G]\n", rs, a, M, Hc);
            printf("num_thread:\t%i\n", nt);
            printf( "Output:\t\t'%s'\n", out_dir.c_str());
        }
        else
        {
            string file_name = out + "sim_param.json";
            ofstream o(file_name);

            json j = {
                {"mesh_num", mesh_num},
                {"mesh_num_pwr", mesh_num_pwr},
                {"Ng", Ng},
                {"par_num", par_num},
                {"box_size", box_size},
                {"redshift", z_in},
                {"redshift_0", z_out},
                {"time_step", db},
                {"app", app},
                {"pwr_type", power.pwr_type},
                {"A", power.A},
                {"index", power.ns},
                {"sigma8", power.sigma8},
                {"smoothing_k", power.k2_G},
                {"Omega_c", power.Omega_c},
                {"Omega_b", power.Omega_b},
                {"h", power.h},
                {"k_pade", k_par.k_pade},
                {"k_nyquist" ,k_par.nyquist},
                {"viscosity", nu},
                {"cut_radius", rs},
                {"num_thread", nt},
                {"out_dir", out},
                {"results", {}}
            };

            switch (power.config.transfer_function_method)
            { // convert to pyccl transfer_function_types keys
                case ccl_emulator: j["transfer_function_method"] = "emulator"; break;
                case ccl_eisenstein_hu: j["transfer_function_method"] = "eisenstein_hu"; break;
                case ccl_bbks: j["transfer_function_method"] = "bbks"; break;
                case ccl_boltzmann_class: j["transfer_function_method"] = "boltzmann_class"; break;
                case ccl_boltzmann_camb: j["transfer_function_method"] = "boltzmann_camb"; break;
            }
            switch (power.config.matter_power_spectrum_method)
            { // convert to pyccl matter_power_spectrum_types keys
                case ccl_linear: j["matter_power_spectrum_method"] = "linear"; break;
                case ccl_halofit: j["matter_power_spectrum_method"] = "halofit"; break;
                case ccl_halo_model: j["matter_power_spectrum_method"] = "halo_model"; break;
            }
            switch (power.config.mass_function_method)
            { // convert to pyccl mass_function_types keys
                case ccl_tinker: j["mass_function_method"] = "tinker"; break;
                case ccl_tinker10: j["mass_function_method"] = "tinker10"; break;
                case ccl_watson: j["mass_function_method"] = "watson"; break;
                case ccl_angulo: j["mass_function_method"] = "angulo"; break;
            }

            o << setw(2) << j << endl;
        }
    } else printf("WARNING! Simulation parameters are not initialized!\n");
}


void Sim_Param::print_info() const
{
	Sim_Param::print_info("", "");
}

/**
 * @class:	App_Var_base
 * @brief:	class containing variables for approximations
 */
 
App_Var_base::App_Var_base(const Sim_Param &sim, string app_str):
	err(0), step(0), print_every(sim.print_every),
	b(sim.b_in), b_init(1.), b_out(sim.b_out), db(sim.db), z_suffix_const(app_str),
	app_field(3, Mesh(sim.mesh_num)),
	power_aux (sim.mesh_num_pwr),
	pwr_spec_binned(sim.bin_num), pwr_spec_binned_0(sim.bin_num), corr_func_binned(sim.bin_num / 2),
	track(4, sim.mesh_num/sim.Ng),
	dens_binned(500), is_init_pwr_spec_0(false)
{
	// FFTW PREPARATION
	err = !fftw_init_threads();
	if (err){
		printf("Errors during multi-thread initialization!\n");
		throw err;
	}
	fftw_plan_with_nthreads(sim.nt);
	p_F = fftw_plan_dft_r2c_3d(sim.mesh_num, sim.mesh_num, sim.mesh_num, app_field[0].real(),
        app_field[0].complex(), FFTW_ESTIMATE);
	p_B = fftw_plan_dft_c2r_3d(sim.mesh_num, sim.mesh_num, sim.mesh_num, app_field[0].complex(),
        app_field[0].real(), FFTW_ESTIMATE);
    p_F_pwr = fftw_plan_dft_r2c_3d(sim.mesh_num_pwr, sim.mesh_num_pwr, sim.mesh_num_pwr, power_aux.real(),
		power_aux.complex(), FFTW_ESTIMATE);
	p_B_pwr = fftw_plan_dft_c2r_3d(sim.mesh_num_pwr, sim.mesh_num_pwr, sim.mesh_num_pwr, power_aux.complex(),
		power_aux.real(), FFTW_ESTIMATE);
}

App_Var_base::~App_Var_base()
{
	// FFTW CLEANUP
	fftw_destroy_plan(p_F);
    fftw_destroy_plan(p_B);
    fftw_destroy_plan(p_F_pwr);
	fftw_destroy_plan(p_B_pwr);
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
	b += db;
}
 
template <class T> 
void App_Var_base::print(const Sim_Param &sim, std::string out_dir_app, T* particles)
{
    /* Printing positions */
    // print_par_pos_cut_small(particles, sim, out_dir_app, z_suffix());
    // print_track_par(track, sim, out_dir_app, z_suffix());

    /* Printing density */
    get_rho_from_par(particles, &power_aux, sim);
    gen_dens_binned(power_aux, dens_binned, sim);    
    // print_rho_map(power_aux, sim, out_dir_app, z_suffix());
    // print_dens_bin(dens_binned, sim.mesh_num, out_dir_app, z_suffix());

    /* Printing power spectrum */
    fftw_execute_dft_r2c(p_F_pwr, power_aux);
    pwr_spec_k(sim, power_aux, &power_aux);
    pwr_spec_binned.resize(sim.bin_num);
    gen_pow_spec_binned(sim, power_aux, &pwr_spec_binned);
    print_pow_spec(pwr_spec_binned, out_dir_app, "_par" + z_suffix());
    if (!is_init_pwr_spec_0){
        pwr_spec_binned_0 = pwr_spec_binned;
        b_init = b;
        is_init_pwr_spec_0 = true;
    }
    print_pow_spec_diff(pwr_spec_binned, pwr_spec_binned_0, b / b_init, out_dir_app, z_suffix());

    /* Print extrapolated power spectrum */
    Extrap_Pk P_k(pwr_spec_binned, sim);
    gen_pow_spec_binned_from_extrap(sim, P_k, &pwr_spec_binned);
    print_pow_spec(pwr_spec_binned, out_dir_app, "_extrap" + z_suffix());

    /* Printing correlation function */
    corr_func_binned.resize(sim.bin_num);
    switch(sim.corr_int)
    {
        case corr_int_type::FFT:
            power_aux.reset_im(); // P(k) is a real function
            fftw_execute_dft_c2r(p_B_pwr, power_aux);
            gen_corr_func_binned(sim, power_aux, &corr_func_binned);
            print_corr_func(corr_func_binned, out_dir_app, "_fft" + z_suffix());
            break;
        case corr_int_type::PP:
            gen_corr_func_binned_pp(sim, particles, &corr_func_binned);
            print_corr_func(corr_func_binned, out_dir_app, "_pp" + z_suffix());
            break;
        case corr_int_type::FFTLOG:
            break;
        case corr_int_type::QAGI:
            gen_corr_func_binned_gsl_qagi(sim, P_k, &corr_func_binned);
            print_corr_func(corr_func_binned, out_dir_app, "_gsl_qagi" + z_suffix());
            break;
        case corr_int_type::QAWO:
            gen_corr_func_binned_gsl_qawo(sim, P_k, &corr_func_binned);
            print_corr_func(corr_func_binned, out_dir_app, "_gsl_qawo" + z_suffix());
            break;
        case corr_int_type::QAWF:
            gen_corr_func_binned_gsl_qawf(sim, P_k, &corr_func_binned);
            print_corr_func(corr_func_binned, out_dir_app, "_gsl_qawf" + z_suffix());
            break;
    }
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

void App_Var::print(const Sim_Param &sim, std::string out_dir_app) {App_Var_base::print(sim, out_dir_app, particles);}

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

void App_Var_v::print(const Sim_Param &sim, std::string out_dir_app) {App_Var_base::print(sim, out_dir_app, particles);}

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

void App_Var_AA::print(const Sim_Param &sim, std::string out_dir_app) {App_Var_base::print(sim, out_dir_app, particles);}

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

void App_Var_FP_mod::print(const Sim_Param &sim, std::string out_dir_app) {App_Var_base::print(sim, out_dir_app, particles);}

/**
 * @class LinkedList
 * @brief class handling linked lists
 */


LinkedList::LinkedList(unsigned par_num, int m, double hc):
	par_num(par_num), Hc(hc), LL(par_num), HOC(m, m, m) {}
	
void LinkedList::get_linked_list(Particle_v* particles)
{
	HOC.assign(-1);
	for (unsigned i = 0; i < par_num; i++)
	{
		LL[i] = HOC(Vec_3D<int>(particles[i].position/Hc));
		HOC(Vec_3D<int>(particles[i].position/Hc)) = i;
	}
}

template class Vec_3D<int>;
template class Vec_3D<double>;
template Vec_3D<double> operator+ (Vec_3D<double>, const Vec_3D<double>&);
template Vec_3D<double> operator- (Vec_3D<double>, const Vec_3D<double>&);
template Vec_3D<double> operator* (Vec_3D<double>, double);
template Vec_3D<double> operator* (double, Vec_3D<double>);
template Vec_3D<double> operator/ (Vec_3D<double>, double);
template Vec_3D<int>::operator Vec_3D<double>() const;
template Vec_3D<double>::operator Vec_3D<int>() const;
template class Mesh_base<int>;
template class Mesh_base<double>;
template void Tracking::update_track_par(Particle_x* particles);
template void Tracking::update_track_par(Particle_v* particles);
template void App_Var_base::print(const Sim_Param&, std::string, Particle_x*);
template void App_Var_base::print(const Sim_Param&, std::string, Particle_v*);

#ifdef TEST
#include "test_core.cpp"
#endif