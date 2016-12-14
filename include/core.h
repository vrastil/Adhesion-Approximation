#pragma once

#include "stdafx.h"
#include <fftw3.h>

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
	~Mesh();
	
	// VARIABLES
	const int N, length; // acces dimensions and length of mesh
	
	// METHODS
	inline double* real() const { return data;} // acces data
	inline fftw_complex* complex() const { return reinterpret_cast<fftw_complex*>(data);}
	
	// OPERATORS
	inline double &operator[](int i){ return data[i]; }
	inline const double &operator[](int i) const{ return data[i]; }
	
	inline double& operator()(int i, int j, int k){ return data[i*N*(N+2)+j*(N+2)+k]; }
	inline const double& operator()(int i, int j, int k) const{ return data[i*N*(N+2)+j*(N+2)+k]; }
	
	Mesh& operator+=(const double& rhs);
	Mesh& operator-=(const double& rhs){ return *this+=-rhs; }
	Mesh& operator*=(const double& rhs);
	Mesh& operator/=(const double& rhs);
};

/**
 * @class:	Vec_3D
 * @brief:	class handling basic 3D-vector functions
 */

class Vec_3D
{
public:
	// CONSTRUCTORS & DESTRUCTOR
	Vec_3D();
	Vec_3D(double x, double y, double z);
	
	// VARIABLES
	double x, y, z;
		
	// OPERATORS
	double& operator[](int i);
	const double& operator[](int i) const;
};

/**
 * @class:	Particle_x
 * @brief:	class handling particles (position only)
 * @acces:	operator [] to get position coordinates
 */

class Particle_x
{
public:
	// CONSTRUCTORS & DESTRUCTOR
	Particle_x(double x, double y, double z):
		position(x, y, z) {};
	
	// VARIABLES
	Vec_3D position;
	
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
	// CONSTRUCTORS & DESTRUCTOR
	Particle_v();
	Particle_v(double x, double y, double z, double vx, double vy, double vz):
		Particle_x(x, y, z),
		velocity(vx, vy, vz) {};
	
	// VARIABLES
	Vec_3D velocity;

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
	std::string out_dir_app;
	Pow_Spec_Param power;
	
	// METHODS
	int init(int ac, char* av[]);
	void print_info();
	
protected:
	bool is_init = 0;
};