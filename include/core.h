#pragma once

#include "stdafx.h"
#include <fftw3.h>

/**
 * @class:	Vec_3D<T>
 * @brief:	class handling basic 3D-vector functions
 */

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
		
	// OPERATORS
	T& operator[](int i)
	{
		switch(i)
		{
			case 0 : return x;
			case 1 : return y;
			case 2 : return z;
			default:
			{
				printf("Invalid acces in class Vec_3D. Invalid postion '%d'.\n", i);
				if (i < 0) return x;
				else return z;
			}
		}
	}
	const T& operator[](int i) const
	{
		switch(i)
		{
			case 0 : return x;
			case 1 : return y;
			case 2 : return z;
			default:
			{
				printf("Invalid acces in class Vec_3D. Invalid postion '%d'.\n", i);
				if (i < 0) return x;
				else return z;
			}
		}
	}
	Vec_3D<T>& operator+=(const Vec_3D<T>& rhs)
	{
		x+=rhs.x;
		y+=rhs.y;
		z+=rhs.z;
		return *this;
	}
	friend Vec_3D<T> operator+(Vec_3D<T> lhs, const Vec_3D<T>& rhs)
	{
		lhs += rhs;
		return lhs;
	}
	Vec_3D<T>& operator*=(T rhs)
	{
		x*=rhs;
		y*=rhs;
		z*=rhs;
		return *this;
	}
	friend Vec_3D<T> operator*(Vec_3D<T> lhs, T rhs)
	{
		lhs *= rhs;
		return lhs;
	}
	
	template<class U>
	explicit operator Vec_3D<U>() const
	{
		Vec_3D<U> lhs;
		lhs.x = static_cast<U>(this->x);
		lhs.y = static_cast<U>(this->y);
		lhs.z = static_cast<U>(this->z);
		return lhs;
	}
};

/**
 * @class:	Mesh
 * @brief:	class handling basic mesh functions, the most important are creating and destroing the underlying data structure
 *			creates a mesh of N*N*(N+2) cells
 */

class Mesh
{
private:
	// VARIABLES
	double* data;
	
public:
	// CONSTRUCTORS & DESTRUCTOR
	Mesh(int n);
	Mesh(const Mesh& that);
	~Mesh();
	
	// VARIABLES
	int N, length; // acces dimensions and length of mesh
	
	// METHODS
	inline double* real() const { return data;} // acces data
	inline fftw_complex* complex() const { return reinterpret_cast<fftw_complex*>(data);}
	
	// OPERATORS
	inline double &operator[](int i){ return data[i]; }
	inline const double &operator[](int i) const{ return data[i]; }
	
	inline double& operator()(int i, int j, int k){ return data[i*N*(N+2)+j*(N+2)+k]; }
	inline const double& operator()(int i, int j, int k) const{ return data[i*N*(N+2)+j*(N+2)+k]; }
	
	double& operator()(Vec_3D<int> pos);
	const double& operator()(Vec_3D<int> pos) const;
	
	Mesh& operator+=(const double& rhs);
	Mesh& operator-=(const double& rhs){ return *this+=-rhs; }
	Mesh& operator*=(const double& rhs);
	Mesh& operator/=(const double& rhs);
	
	friend void swap(Mesh& first, Mesh& second);
	Mesh& operator=(Mesh& that);
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
	Tracking(int num_track_par, int par_num);
	
	// VARIABLES
	int num_track_par; // square root of number of tracking particles
	std::vector<int> par_ids;
	std::vector<std::vector<Particle_x>> par_pos;
	
	// METHODS
	int num_step(){return par_pos.size();};
	void update_track_par(Particle_x* particles);
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
	double b_in, b_out;
	double nu;
	int nt;
	std::string out_dir;
//	std::string out_dir_app;
	Pow_Spec_Param power;
	
	// METHODS
	int init(int ac, char* av[]);
	void print_info();
	
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
	int err = 0, step = 0, print_every = 1;
	double b, b_out, db;
	const std::string z_suffix_const;
	std::vector<Mesh> app_field;
	Mesh power_aux;
	std::vector<fftw_complex> pwr_spec_binned, pwr_spec_binned_0;
	fftw_plan p_F, p_B;
	Tracking track;
	
	// METHODS
	double z(){ return 1./b - 1.;}
	bool integrate(){return (b <= b_out) && (db > 0);}
	bool printing(){ return ((step % print_every) == 0) or (b == b_out); }
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
};