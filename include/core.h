/**
 * @file:	core.h
 * @brief:	class definitions
 */

#pragma once

#include "stdafx.h"
#include <ccl.h>

constexpr double PI = M_PI;
/**
 * @class:	Vec_3D<T>
 * @brief:	class handling basic 3D-vector functions
 */
 
template <typename T> int sgn(T val)
{
	return (T(0) < val) - (val < T(0));
}

template <typename T>
class Vec_3D
{
public:
	// CONSTRUCTORS
	Vec_3D(){};
	Vec_3D(T x, T y, T z):
	vec({x, y, z}) {};
    
    // VARIABLES
    std::array<T, 3> vec;

    // ELEMENT ACCESS
    inline T& operator[](int i){ return vec[i]; }
    inline const T& operator[](int i) const { return vec[i]; }

    // ITERATORS
    typedef typename std::array<T, 3>::iterator iterator;
    typedef typename std::array<T, 3>::const_iterator const_iterator;

    inline iterator begin() noexcept { return vec.begin(); }
    inline const_iterator cbegin() const noexcept { return vec.cbegin(); }
    inline iterator end() noexcept { return vec.end(); }
    inline const_iterator cend() const noexcept { return vec.end(); }
    
	// METHODS
	double norm() const;
	inline void fill(const T& value){ vec.fill(value); }
		
	// OPERATORS	
	Vec_3D<T>& operator+=(const Vec_3D<T>& rhs);
    Vec_3D<T>& operator-=(const Vec_3D<T>& rhs);
    Vec_3D<T>& operator+=(T rhs);
	Vec_3D<T>& operator-=(T rhs);
	Vec_3D<T>& operator*=(T rhs);
	Vec_3D<T>& operator/=(T rhs);
	template<class U>
    explicit operator Vec_3D<U>() const;
};

// NON-MEMBER FUNCTIONS
template <typename T> Vec_3D<T> operator+(Vec_3D<T> lhs, const Vec_3D<T>& rhs);
template <typename T> Vec_3D<T> operator-(Vec_3D<T> lhs, const Vec_3D<T>& rhs);
template <typename T> Vec_3D<T> operator*(Vec_3D<T> lhs, T rhs);
template <typename T> Vec_3D<T> operator*(T lhs, Vec_3D<T> rhs);
template <typename T> Vec_3D<T> operator+(Vec_3D<T> lhs, T rhs);
template <typename T> Vec_3D<T> operator+(T lhs, Vec_3D<T> rhs);
template <typename T> Vec_3D<T> operator-(Vec_3D<T> lhs, T rhs);
template <typename T> Vec_3D<T> operator-(T lhs, Vec_3D<T> rhs);
template <typename T> Vec_3D<T> operator/(Vec_3D<T> lhs, T rhs);

template <typename T> inline bool operator==(const Vec_3D<T>& lhs, const Vec_3D<T>& rhs){return lhs.vec == rhs.vec; }
template <typename T> inline bool operator!=(const Vec_3D<T>& lhs, const Vec_3D<T>& rhs){ return lhs.vec != rhs.vec; }
template <typename T> inline bool operator<(const Vec_3D<T>& lhs, const Vec_3D<T>& rhs){ return lhs.vec < rhs.vec; }
template <typename T> inline bool operator<=(const Vec_3D<T>& lhs, const Vec_3D<T>& rhs){ return lhs.vec <= rhs.vec; }
template <typename T> inline bool operator>(const Vec_3D<T>& lhs, const Vec_3D<T>& rhs){ return lhs.vec > rhs.vec; }
template <typename T> inline bool operator>=(const Vec_3D<T>& lhs, const Vec_3D<T>& rhs){ return lhs.vec >= rhs.vec; }

/**
 * @class:	Mesh_base<T>
 * @brief:	class handling basic mesh functions, the most important are creating and destroing the underlying data structure
 *			creates a mesh of N1*N2*N3 cells
 */

template <typename T>
class Mesh_base
{
public:
	// CONSTRUCTORS & DESTRUCTOR
	Mesh_base(unsigned n1, unsigned n2, unsigned n3);
    Mesh_base(const Mesh_base& that);
    Mesh_base(Mesh_base<T>&& that) noexcept;
    Mesh_base& operator=(Mesh_base that) &;
    template <class U> friend void swap(Mesh_base<U>& first, Mesh_base<U>& second);
	~Mesh_base();
	
	// VARIABLES
	unsigned N1, N2, N3, length; // acces dimensions and length of mesh
	
	// METHODS
	inline T* real() const { return data;} // acces data
	void assign(T val);
	
	// OPERATORS
	inline T &operator[](int i){ return data[i]; }
	inline const T &operator[](int i) const{ return data[i]; }
	
	inline T& operator()(int i, int j, int k){ return data[i*N2*N3+j*N3+k]; }
	inline const T& operator()(int i, int j, int k) const{ return data[i*N2*N3+j*N3+k]; }
	
	inline T& operator()(int i, int j){ return data[i*N3+j]; }
	inline const T& operator()(int i, int j) const{ return data[i*N3+j]; }
	
	T& operator()(Vec_3D<int> pos);
	const T& operator()(Vec_3D<int> pos) const;
	
	Mesh_base& operator+=(const T& rhs);
	Mesh_base& operator-=(const T& rhs){ return *this+=-rhs; }
	Mesh_base& operator*=(const T& rhs);
	Mesh_base& operator/=(const T& rhs);

protected:
	// VARIABLES
	T* data;
};

// template <typename T> void swap(Mesh_base<T>& first, Mesh_base<T>& second);

/**
 * @class:	Mesh
 * @brief:	creates a mesh of N*N*(N+2) cells
 */

class Mesh : public Mesh_base<double>
{
public:
	// CONSTRUCTORS & DESTRUCTOR
    Mesh(unsigned n);
    // default constructors / assignemnts belowe
    // Mesh(const Mesh& that);
    // Mesh(Mesh&& that) noexcept;
    // Mesh& operator=(Mesh that) &;
    // friend void swap(Mesh& first, Mesh& second);
	
	// VARIABLES
	unsigned N; // acces dimension of mesh
	
	// METHODS
    inline fftw_complex* complex() const { return reinterpret_cast<fftw_complex*>(data);}
    void reset_part(bool part);
    inline void reset_re() { reset_part(0); }
    inline void reset_im() { reset_part(1); }
    
	// OPERATORS
	using Mesh_base<double>::operator ();
	double& operator()(Vec_3D<int> pos);
	const double& operator()(Vec_3D<int> pos) const;
	
	Mesh& operator+=(const double& rhs);
	Mesh& operator-=(const double& rhs){ return *this+=-rhs; }
	Mesh& operator*=(const double& rhs);
	Mesh& operator/=(const double& rhs);
};

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
	Particle_x(double x, double y, double z):
	position(x,y,z) {};
	Particle_x(Vec_3D<double> position):
	position(position) {};
	
	// VARIABLES
	Vec_3D<double> position;
	
	// OPERATORS
	inline double &operator[](int i){ return position[i]; }
	inline const double& operator[](int i) const{ return position[i]; }
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
	Particle_v(double x, double y, double z, double vx, double vy, double vz):
		Particle_x(x,y,z), velocity(vx,vy,vz) {};
	Particle_v(Vec_3D<double> position, Vec_3D<double> velocity):
		Particle_x(position), velocity(velocity) {};
	
	// VARIABLES
	Vec_3D<double> velocity;

	// OPERATORS
	double &operator()(int i){ return velocity[i]; }
	const double& operator()(int i) const{ return velocity[i]; }
};

/**
 * @class:	Pow_Spec_Param
 * @brief:	class storing parameters for power spectrum
 */

enum e_power_spec { power_law_T = 0, power_law = 1, flat = 2, single = 3,
                    ccl_EH = 4, ccl_BBKS = 5};

class Pow_Spec_Param
{
public:
    // CONSTRUCTOR
    void init();

    // DESTRUCTOR
    ~Pow_Spec_Param();

	double A = 1, ns, k2_G, sigma8;
    e_power_spec pwr_type;
    int pwr_type_i;

    // COSMOLOGY (flat LCDM)
    double Omega_c, Omega_b, H0, h;
    inline double Omega_m() const { return Omega_c + Omega_b; }
        
    // CCL VARIABLES
    ccl_configuration config;
    ccl_parameters params;
    ccl_cosmology* cosmo;

protected:
	bool is_init = 0;
};

struct Range
{
    double lower, upper;
};

struct Wavenumber_param
{
    // k-range where to use (linear) interpolation and k-range in which print 'pwr_spec_extrap_*'
    Range k_interp, k_print;

    // Nyquist frequencies of potential mesh, analyses mesh and particle separation
    std::map<std::string,double> nyquist;

    // use Padé approximant R[m=0; n=pade_order-1], i.e. pade_order of DOF
    unsigned pade_order = 3;

    // store values of k where the Padé approximant was evaluated
    mutable std::vector<double> k_pade;
};

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
	inline const unsigned num_step() const{return par_pos.size();};
	template <class T>  void update_track_par(T* particles);
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

/**
 * @class:	Sim_Param
 * @brief:	class storing simulation parameters
 */

class Sim_Param
{
public:
    // CONSTRUCTOR
    Sim_Param(int ac, char* av[]);
    /* SIMULATION BOX*/
    unsigned par_num, mesh_num, mesh_num_pwr, Ng, Ng_pwr, print_every;
    double box_size;
    unsigned order = 1, bin_num = 200;
    unsigned long seed = 12345678;

    /* TIME */
	double z_in, z_out;
    double b_in, b_out, db;
    
	double nu;
	unsigned nt;
	double rs, Hc, a;
	unsigned M;
	std::string out_dir; ///< where to save output of the simulation
    Pow_Spec_Param power; ///< all information about our cosmology
    Wavenumber_param k_par; ///< all the important k-values we need

    /* CORRELATION FUNCTION */
    Range x_corr; ///< range in which compute the correlation function
    corr_int_type corr_int;

    bool comp_ZA, comp_FF, comp_FP, comp_AA, comp_FP_pp;
	
	// METHODS
    void print_info(std::string out, std::string app) const;
	void print_info() const;
	inline const double x_0() const{return (double)box_size/mesh_num;}
    inline const double x_0_pwr() const{return (double)box_size/mesh_num_pwr;}
    
protected:
	bool is_init = 0;
};

/**
 * @class:	Data_x_y
 * @brief:	class containing data [x, y(x)]
 */

template <typename T>
class Data_x_y
{
public:
	// CONSTRUCTORS
    Data_x_y(){};
    Data_x_y(size_t size) : x(size), y(size) {;}
    
    // VARIABLES
    std::vector<T> x;
    std::vector<T> y;

    // CAPACITY
    size_t size() const noexcept{ return x.size(); }
    void resize (size_t n){
        x.resize(n);
        y.resize(n);
    }
    void resize (size_t n, double val){
        x.resize(n, val);
        y.resize(n, val);
    }
    void erase(unsigned index){
        x.erase(x.begin() + index);
        y.erase(y.begin() + index);
    }

    // MODIFIERS
    void fill(const T& val){
        std::fill(x.begin(), x.end(), val);
        std::fill(y.begin(), y.end(), val);
    }
};

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

	int err, step, print_every;
	double b, b_init, b_out, db;
    const std::string app_str, z_suffix_const, out_dir_app;
    
    // LARGE FIELDS
	std::vector<Mesh> app_field;
    std::vector<Mesh> power_aux;
    T* particles;

	Data_x_y<double> pwr_spec_binned, pwr_spec_binned_0, vel_pwr_spec_binned_0, corr_func_binned;
	fftw_plan p_F, p_B, p_F_pwr, p_B_pwr;
	Tracking track;
	std::vector<int> dens_binned;
	
	// METHODS
	inline double z() const{ return 1./b - 1.;}
	inline double b_half() const { return b - db/2.; }
	inline bool integrate() const { return (b <= b_out) && (db > 0);}
	inline bool printing() const { return ((step % print_every) == 0) or (b == b_out); }
    void print();
    void upd_time();
    void print_mem() const;
    void print_info() const;
	
	std::string z_suffix();
	
protected:	
    std::stringstream z_suffix_num;
    bool is_init_pwr_spec_0, is_init_vel_pwr_spec_0;
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
	LinkedList(unsigned par_num, int m, double hc);
	
	// VARIABLES
	unsigned par_num;
	double Hc;
	std::vector<int> LL;
	Mesh_base<int> HOC;
	
	// METHODS
	void get_linked_list(Particle_v* particles);
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
};