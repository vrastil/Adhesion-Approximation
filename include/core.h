/**
 * @file:	core.h
 * @brief:	class definitions
 */

#pragma once

#include "stdafx.h"
#include <fftw3.h>

typedef std::array<double, 2> double_2;

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
	x(x), y(y), z(z) {};
	
	// VARIABLES
	T x, y, z;
	
	// METHODS
	inline double norm(){ return sqrt(x*x+y*y+z*z); }
	void assign(T x_, T y_, T z_);
		
	// OPERATORS
	T& operator[](int i);
	const T& operator[](int i) const;
	Vec_3D<T>& operator+=(const Vec_3D<T>& rhs);
	Vec_3D<T>& operator-=(const Vec_3D<T>& rhs);
	Vec_3D<T>& operator*=(T rhs);
	Vec_3D<T>& operator/=(T rhs);
	template<class U>
	explicit operator Vec_3D<U>() const;
};

template <typename T> Vec_3D<T> operator+(Vec_3D<T> lhs, const Vec_3D<T>& rhs);
template <typename T> Vec_3D<T> operator-(Vec_3D<T> lhs, const Vec_3D<T>& rhs);
template <typename T> Vec_3D<T> operator*(Vec_3D<T> lhs, T rhs);
template <typename T> Vec_3D<T> operator*(T lhs, Vec_3D<T> rhs);
template <typename T> Vec_3D<T> operator/(Vec_3D<T> lhs, T rhs);

/**
 * @class:	Mesh_base<T>
 * @brief:	class handling basic mesh functions, the most important are creating and destroing the underlying data structure
 *			creates a mesh of N1*N2*N3 cells
 */

template <typename T>
class Mesh_base
{
protected:
	// VARIABLES
	T* data;
	
public:
	// CONSTRUCTORS & DESTRUCTOR
	Mesh_base(int n1, int n2, int n3);
	Mesh_base(const Mesh_base& that);
	~Mesh_base();
	
	// VARIABLES
	int N1, N2, N3, length; // acces dimensions and length of mesh
	
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
	
//	friend void swap(Mesh_base& first, Mesh_base& second);
	Mesh_base& operator=(const Mesh_base& other);
};

template <typename T> void swap(Mesh_base<T>& first, Mesh_base<T>& second);

/**
 * @class:	Mesh
 * @brief:	creates a mesh of N*N*(N+2) cells
 */

class Mesh : public Mesh_base<double>
{
public:
	// CONSTRUCTORS & DESTRUCTOR
	Mesh(int n);
	Mesh(const Mesh& that);
	
	// VARIABLES
	int N; // acces dimension of mesh
	
	// METHODS
	inline fftw_complex* complex() const { return reinterpret_cast<fftw_complex*>(data);}
	
	// OPERATORS
	using Mesh_base<double>::operator ();
	double& operator()(Vec_3D<int> pos);
	const double& operator()(Vec_3D<int> pos) const;
	
	Mesh& operator+=(const double& rhs);
	Mesh& operator-=(const double& rhs){ return *this+=-rhs; }
	Mesh& operator*=(const double& rhs);
	Mesh& operator/=(const double& rhs);
	
	friend void swap(Mesh& first, Mesh& second);
	Mesh& operator=(const Mesh& other);
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
	position(position.x,position.y,position.z) {};
	
	// VARIABLES
	Vec_3D<double> position;
	
	// OPERATORS
	double &operator[](int i){ return position[i]; }
	const double& operator[](int i) const{ return position[i]; }
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
		Particle_x(position), velocity(velocity.x, velocity.y, velocity.z) {};
	
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

enum e_power_spec { power_law_T = 0, power_law = 1, flat = 2, single = 3};

struct Pow_Spec_Param
{
	double A = 1, ns, k2_G, s8;
	e_power_spec pwr_type;
};

/**
 * @class:	Tracking
 * @brief:	class storing info about tracked particles
 */

class Tracking
{
public:
	// CONSTRUCTORS
	Tracking(int sqr_num_track_par, int par_num_per_dim);
	
	// VARIABLES
	int sqr_num_track_par, num_track_par; // square root of number of tracking particles
	std::vector<int> par_ids;
	std::vector<std::vector<Particle_x>> par_pos;
	
	// METHODS
	inline const int num_step() const{return par_pos.size();};
	void update_track_par(Particle_x* particles);
	void update_track_par(Particle_v* particles);
};

/**
 * @class:	Sim_Param
 * @brief:	class storing simulation parameters
 */
 
class Sim_Param
{
public:
	// VARIABLES
	int par_num, mesh_num, Ng, box_size;
	int order = 1, bin_num = 100;
	double k_min, k_max;
	unsigned long seed = 12345678;
	double z_in, z_out;
	double b_in, b_out, db;
	double nu;
	int nt;
	double rs, Hc, a;
	int M;
	std::string out_dir;
	Pow_Spec_Param power;
	bool comp_ZA, comp_FF, comp_FP, comp_AA, comp_FP_pp;
	
	// METHODS
	int init(int ac, char* av[]);
    void print_info(std::string out, std::string app) const;
	void print_info() const;
	inline const double x_0() const{return box_size/mesh_num;}
	
protected:
	bool is_init = 0;
};

/**
 * @class:	App_Var_base
 * @brief:	class containing core variables for approximations
 */
 
class App_Var_base
{
public:
	// CONSTRUCTORS
	App_Var_base(const Sim_Param &sim, std::string app_str);
	
	// DESTRUCTOR
	~App_Var_base();
	
	// VARIABLES
	int err, step, print_every;
	double b, b_out, db;
	const std::string z_suffix_const;
	std::vector<Mesh> app_field;
	Mesh power_aux;
	std::vector<double_2> pwr_spec_binned, pwr_spec_binned_0;
	fftw_plan p_F, p_B;
	Tracking track;
	std::vector<int> dens_binned;
	
	// METHODS
	inline double z(){ return 1./b - 1.;}
	inline double b_half() {return b - db/2.; }
	inline bool integrate(){return (b <= b_out) && (db > 0);}
	inline bool printing(){ return ((step % print_every) == 0) or (b == b_out); }
	void print_x(const Sim_Param &sim, std::string out_dir_app, Particle_x* particles);
	void print_v(const Sim_Param &sim, std::string out_dir_app, Particle_v* particles);
	void upd_time();
	
	std::string z_suffix();
	
protected:	
	std::stringstream z_suffix_num;
};

/**
 * @class:	App_Var
 * @brief:	class containing variables for approximations with particle positions only
 */
 
 class App_Var: public App_Var_base
{
public:
	// CONSTRUCTORS & DESTRUCTOR
	App_Var(const Sim_Param &sim, std::string app_str);
	~App_Var();
	
	// VARIABLES
	Particle_x* particles;
	
	// METHODS
	inline void print(const Sim_Param &sim, std::string out_dir_app) {print_x(sim, out_dir_app, particles);}
};

/**
 * @class:	App_Var_v
 * @brief:	class containing variables for approximations with particle velocities
 */
 
 class App_Var_v: public App_Var_base
{
public:
	// CONSTRUCTORS & DESTRUCTOR
	App_Var_v(const Sim_Param &sim, std::string app_str);
	~App_Var_v();
	
	// VARIABLES
	Particle_v* particles;
	
	// METHODS
	inline void print(const Sim_Param &sim, std::string out_dir_app) {print_v(sim, out_dir_app, particles);}
};

/**
 * @class:	App_Var_AA
 * @brief:	class containing variables for adhesion approximation
 */
 
 class App_Var_AA: public App_Var_base
{
public:
	// CONSTRUCTORS & DESTRUCTOR
	App_Var_AA(const Sim_Param &sim, std::string app_str);
	~App_Var_AA();
	
	// VARIABLES
	Particle_x* particles;
	Mesh expotential;
	
	// METHODS
	inline void print(const Sim_Param &sim, std::string out_dir_app) {print_x(sim, out_dir_app, particles);}
};

/**
 * @class LinkedList
 * @brief class handling linked lists
 */

class LinkedList
{
public:
	// CONSTRUCTORS & DESTRUCTOR
	LinkedList(int par_num, int m, double hc);
	
	// VARIABLES
	int par_num;
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
 
 class App_Var_FP_mod: public App_Var_base
{
public:
	// CONSTRUCTORS & DESTRUCTOR
	App_Var_FP_mod(const Sim_Param &sim, std::string app_str);
	~App_Var_FP_mod();
	
	// VARIABLES
	Particle_v* particles;
	LinkedList linked_list;
	
	// METHODS
	inline void print(const Sim_Param &sim, std::string out_dir_app) {print_v(sim, out_dir_app, particles);}
};

#include "core.hpp"