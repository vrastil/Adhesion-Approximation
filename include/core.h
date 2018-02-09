/**
 * @file:	core.h
 * @brief:	class definitions
 */

#pragma once
#include "stdafx.h"
#include "vec_3d.hpp"

/*********************************************
 * single / double / long double definitions *
 ********************************************/

#ifndef PRECISION
#define PRECISION 2 // default double precision
#endif

#if PRECISION == 1
typedef float FTYPE;
#define MAKE_FFTW_NAME(FUNC_NAME) fftwf_ ## FUNC_NAME
#elif PRECISION == 2
typedef double FTYPE;
#define MAKE_FFTW_NAME(FUNC_NAME) fftw_ ## FUNC_NAME
#elif PRECISION == 3
typedef long double FTYPE;
#define MAKE_FFTW_NAME(FUNC_NAME) fftwl_ ## FUNC_NAME
#endif

#define FFTW_PLAN_TYPE MAKE_FFTW_NAME(plan)
#define FFTW_DEST_PLAN MAKE_FFTW_NAME(destroy_plan)
#define FFTW_COMPLEX_TYPE MAKE_FFTW_NAME(complex)
#define FFTW_PLAN_R2C MAKE_FFTW_NAME(plan_dft_r2c_3d)
#define FFTW_PLAN_C2R MAKE_FFTW_NAME(plan_dft_c2r_3d)
#define FFTW_PLAN_OMP MAKE_FFTW_NAME(plan_with_nthreads)
#define FFTW_PLAN_OMP_INIT MAKE_FFTW_NAME(init_threads)
#define FFTW_PLAN_OMP_CLEAN MAKE_FFTW_NAME(cleanup_threads)
#define FFTW_EXEC_R2C MAKE_FFTW_NAME(execute_dft_r2c)
#define FFTW_EXEC_C2R MAKE_FFTW_NAME(execute_dft_c2r)


constexpr FTYPE PI = (FTYPE)M_PI;

inline float pow(float base, unsigned exp)
{
    float result = 1.f;
    while (exp)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }
    return result;
}

template <typename T> inline T pow2(T base){ return base*base; } //< most often used

template <typename T> inline int sgn(T val)
{
	return (T(0) < val) - (val < T(0));
}

// ignore following in SWIG wrapper
#ifndef SWIG
#include "mesh.hpp"

/**
 * @class:	Particle_x
 * @brief:	class handling particles (position only)
 * @acces:	operator [] to get position coordinates
 */

class Particle_x
{
public:
	// CONSTRUCTORS
	Particle_x(){};
    // template<typename T>
	// Particle_x(T x, T y, T z):
	// position(x,y,z) {};
    template<typename T>
	Particle_x(Vec_3D<T> position):
	position(position) {};
	
	// VARIABLES
	Vec_3D<FTYPE> position;
	
	// OPERATORS
	FTYPE &operator[](int i){ return position[i]; }
	const FTYPE& operator[](int i) const{ return position[i]; }
};

/**
 * @class:	Particle_v
 * @brief:	class handling particles (with velocitites)
 * @acces:	operator [] to get position coordinates
 * 			operator () to get velocity coordinates
 */

class Particle_v : public Particle_x
{
public:
	// CONSTRUCTORS
	Particle_v(){};
	// Particle_v(FTYPE x, FTYPE y, FTYPE z, FTYPE vx, FTYPE vy, FTYPE vz):
	// 	Particle_x(x,y,z), velocity(vx,vy,vz) {};
    template<typename T, typename U>
	Particle_v(Vec_3D<T> position, Vec_3D<U> velocity):
		Particle_x(position), velocity(velocity) {};
	
	// VARIABLES
	Vec_3D<FTYPE> velocity;

	// OPERATORS
	FTYPE &operator()(int i){ return velocity[i]; }
	const FTYPE& operator()(int i) const{ return velocity[i]; }
};
#endif

class Cosmo_Param
{
public:
    // CONSTRUCTORS, DESTRUCTOR
    void init(); //< lazy constructor
    Cosmo_Param();
    ~Cosmo_Param();

    // CCL VARIABLES
    ccl_configuration config;
    ccl_cosmology* cosmo;

    // COSMOLOGY (flat LCDM)
    FTYPE A = 1, ns, k2_G, sigma8;
    FTYPE Omega_m, Omega_b, H0, h;
    FTYPE Omega_c() const { return Omega_m - Omega_b; }
    FTYPE Omega_L() const { return 1 - Omega_m; }

    // PRECOMPUTED VALUES
    FTYPE D_norm;

    // DEALING WITH GSL 'void* param'
    explicit operator void*() const;
};
void to_json(nlohmann::json&, const Cosmo_Param&);
void from_json(const nlohmann::json&, Cosmo_Param&);

/**
 * @class:	Tracking
 * @brief:	class storing info about tracked particles
 */

class Tracking
{
public:
	// CONSTRUCTOR
	Tracking(int sqr_num_track_par, int par_num_per_dim);
	
	// VARIABLES
	int sqr_num_track_par, num_track_par; // square root of number of tracking particles
	std::vector<int> par_ids;
	std::vector<std::vector<Particle_x>> par_pos;
	
	// METHODS
	const unsigned num_step() const{return par_pos.size();};
	template <class T>  void update_track_par(const std::vector<T>& particles);
};

/**
 * @brief An enum class that defines which integration scheme use in computation of correlation function.

Q - quadrature routine

N - non-adaptive integrator
A - adaptive integrator

G - general integrand (user-defined)
W - weight function with integrand

S - singularities can be more readily integrated
P - points of special difficulty can be supplied
I - infinite range of integration
O - oscillatory weight function, cos or sin
F - Fourier integral
C - Cauchy principal value
*/
enum class corr_int_type { QAGI, QAWO, QAWF, FFT, FFTLOG, PP };

/* SIMULATION BOX*/
struct Box_Opt {
    void init();
    /* cmd args */
    unsigned par_num_1d, mesh_num, mesh_num_pwr;
    FTYPE box_size;
    /* derived param*/
    unsigned par_num, Ng, Ng_pwr;
};


/* INTEGRATION */
struct Integ_Opt {
    void init();
    /* cmd args */
    FTYPE z_in, z_out, db;
    /* derived param*/
    FTYPE b_in, b_out;
};


/* OUTPUT */
struct Out_Opt {
    void init();
    /* cmd args */
    unsigned print_every, bins_per_decade, points_per_10_Mpc;
    std::vector<FTYPE> print_z; //< for which redshifts print output on top of print_every (optional)
    std::string out_dir; //< where to save output of the simulation
    bool print_par_pos, print_dens, print_pwr, print_extrap_pwr, print_corr, print_vel_pwr;
    /* derived param*/
    bool get_rho, get_pwr, get_pk_extrap;
};


/* APPROXIMATIONS */
struct Comp_App {
    /* cmd args */
    bool ZA, FF, FP, AA, FP_pp; //< approximations
    bool chi; //< modified gravities
};


/* APPROXIMATIONS */
struct App_Opt {
    void init(const Box_Opt&);
    /* cmd args */
    FTYPE nu, rs;
    /* derived param*/
    FTYPE Hc, a, nu_dim;
    unsigned M;
};


/* RUN */
struct Run_Opt {
    void init();
    bool simulate();
    /* cmd args */
    unsigned nt, mlt_runs;
    unsigned long seed;
    bool pair;        
    /* other*/
    bool phase;
};

// define Range outside because of SWIG
struct Range { FTYPE lower, upper; };

/* OTHER PARAMETERS */
struct Other_par {
    void init(const Box_Opt&);
    // k-range where to use (linear) interpolation and k-range in which print 'pwr_spec_extrap_*'
    ///range in which compute the correlation function 
    Range k_print, x_corr;
    std::map<std::string,FTYPE> nyquist; //< Nyquist frequencies of potential mesh, analyses mesh and particle separation
};

/* CHAMELEON */
struct Chi_Opt {
    /* cmd args */
    FTYPE beta, n, phi;
};

/**
 * @class:	Sim_Param
 * @brief:	class storing simulation parameters
 */

class Sim_Param
{
public:
    // CONSTRUCTOR
    Sim_Param(int ac, const char* const av[]); //< load from command line arguments
    Sim_Param(std::string file_name); //< load from sim_param.json file

    // VARIABLES
    Box_Opt box_opt;
    Integ_Opt integ_opt;
    Out_Opt out_opt;
    Comp_App comp_app;
    Cosmo_Param cosmo;
    App_Opt app_opt;
    Run_Opt run_opt;
    Other_par other_par;
    Chi_Opt chi_opt;

	// METHODS
    void print_info(std::string out, std::string app) const;
	void print_info() const;
	const FTYPE x_0() const{return box_opt.box_size/box_opt.mesh_num;}
    const FTYPE x_0_pwr() const{return box_opt.box_size/box_opt.mesh_num_pwr;}
    bool simulate() { return run_opt.simulate(); }
};

void to_json(nlohmann::json&, const Box_Opt&);
void to_json(nlohmann::json&, const Integ_Opt&);
void to_json(nlohmann::json&, const App_Opt&);
void to_json(nlohmann::json&, const Run_Opt&);
void to_json(nlohmann::json&, const Out_Opt&);

void from_json(const nlohmann::json&, Box_Opt&);
void from_json(const nlohmann::json&, Integ_Opt&);
void from_json(const nlohmann::json&, App_Opt&);
void from_json(const nlohmann::json&, Run_Opt&);
void from_json(const nlohmann::json&, Out_Opt&);

/**
 * @class:	Data_Vec
 * @brief:	class containing data [x, y,...]
 */

template <typename T, unsigned N>
class Data_Vec
{
public:
    // CONSTRUCTORS
    Data_Vec() = default;
    Data_Vec(size_t size) { data.fill(std::vector<T>(size)); }

    // VARIABLES
    std::array<std::vector<T>, N> data;

    // ELEMENT ACCESS
    std::vector<T>& operator[](int i){ return data[i]; }
    const std::vector<T>& operator[](int i) const { return data[i]; }

    // CAPACITY
    size_t dim() const noexcept{ return data.size(); }
    size_t size() const noexcept{ return data[0].size(); }
    void resize (size_t n){
        for (auto &vec : data) vec.resize(n);
    }
    void resize (size_t n, FTYPE val){
        for (auto &vec : data) vec.resize(n, val);
    }
    void reserve(size_t n){
        for (auto &vec : data) vec.reserve(n);
    }
    void erase(unsigned index){
        for (auto &vec : data) vec.erase(vec.begin() + index);
    }
    // MODIFIERS
    void fill(const T& val){
        for (auto &vec : data) std::fill(vec.begin(), vec.end(), val);
    }
};

// ignore following in SWIG wrapper
#ifndef SWIG

#include "core_power.h"

/**
 * @class:	App_Var
 * @brief:	class containing core variables for approximations
 */

template <class T> 
class App_Var
{
public:
	// CONSTRUCTORS
    App_Var(const Sim_Param &sim, std::string app_str);

	// DESTRUCTOR
	~App_Var();
	
    // VARIABLES
    const Sim_Param &sim;

	int step, print_every;
    FTYPE b, b_out, db;
    FTYPE D_init, dDda_init;
    const std::string app_str, z_suffix_const, out_dir_app;
    
    // LARGE FIELDS
	std::vector<Mesh> app_field;
    std::vector<Mesh> power_aux;
    std::vector<T> particles;

    // OTHER VARIABLES
    Data_Vec<FTYPE, 2> corr_func_binned, pwr_spec_binned, pwr_spec_binned_0, vel_pwr_spec_binned_0;
    Interp_obj pwr_spec_input;
	FFTW_PLAN_TYPE p_F, p_B, p_F_pwr, p_B_pwr;
	Tracking track;
	std::vector<int> dens_binned;
	
	// METHODS
	FTYPE z() const{ return 1/b - 1;}
	FTYPE b_half() const { return b - db/2; }
	bool integrate() const { return (b <= b_out) && (db > 0);}
	bool printing() const { return print_every ? ((step % print_every) == 0) or (b == b_out) : print_every ; }
    void print_output();
    void upd_time();
    void print_mem() const;
    void print_info() const;	
	std::string z_suffix();
	bool is_init_pwr_spec_0, is_init_vel_pwr_spec_0; //< be careful about setting these to true
protected:	
    std::stringstream z_suffix_num;
    uint64_t memory_alloc; // only the largest chunks
};

/**
 * @class:	App_Var_AA
 * @brief:	class containing variables for adhesion approximation
 */
 
 class App_Var_AA: public App_Var<Particle_v>
{
public:
	// CONSTRUCTORS & DESTRUCTOR
	App_Var_AA(const Sim_Param &sim, std::string app_str);

	// VARIABLES
	Mesh expotential;
};

/**
 * @class LinkedList
 * @brief class handling linked lists
 */

class LinkedList
{
public:
	// CONSTRUCTORS & DESTRUCTOR
	LinkedList(unsigned par_num, int m, FTYPE hc);
	
	// VARIABLES
	unsigned par_num;
	FTYPE Hc;
	std::vector<int> LL;
	Mesh_base<int> HOC;
	
	// METHODS
	void get_linked_list(const std::vector<Particle_v>& particles);
};

/**
 * @class:	App_Var_FP_mod
 * @brief:	class containing variables for modified Frozen-potential approximation
 */

class App_Var_FP_mod: public App_Var<Particle_v>
{
public:
	// CONSTRUCTORS & DESTRUCTOR
	App_Var_FP_mod(const Sim_Param &sim, std::string app_str);
	
	// VARIABLES
    LinkedList linked_list;
    Interp_obj fs_interp;
};
#endif
