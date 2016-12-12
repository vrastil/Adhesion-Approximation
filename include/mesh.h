/*
Header for the mesh.cpp file
*/

#pragma once

#include "stdafx.h"
#include <fftw3.h>
#include "threadpool.hpp"

enum e_power_spec { power_law_T = 0, power_law = 1, flat = 2, single = 3};

// *******************
// * CLASS :: c_Mesh *
// *******************

template <class c_Type>
class c_Mesh
{
public:
	// CONSTRUCTORS
	c_Mesh(int n1, int n2, int n3);
	c_Mesh(int n); // mesh n*n*n
	c_Mesh(int n1, int n2, int n3, c_Type init_val);
	c_Mesh(int n, c_Type init_val);
	
	// VARIABLES
	const int N1, N2, N3, len; // acces dimensions and length of mesh
	std::vector<c_Type> p_data;
	
	// METHODS
	c_Type* data(){return p_data.data();}; // acces data
	c_Type get_val(int x, int y, int z) const {return p_data[index(x, y, z)];}
	int index(int x, int y, int z) const { return x*N2*N3+y*N3+z; }
	int index(int y[3]) const { return index(y[0], y[1], y[2]); }
	int index(int xy, int z) const { return xy*N3+z; }
	// get coordinates
	int get_x(int index) const { return index / (N2*N3); }
	int get_y(int index) const { return (index / N3) % N2; }
	int get_z(int index) const { return index % N3; }
	// get coordinates assuming complex numbers (half)
	int get_x_c(int index) const { return index / (N2*(N3/2)); }
	int get_y_c(int index) const { return (index / (N3/2)) % N2; }
	int get_z_c(int index) const{ return index % (N3/2); }
	int get_kx(int index) const;
	int get_ky(int index) const;
	int get_kz(int index) const;
	void get_k_vec(int index, int* k_vec) const;
	int get_k_sq(int index) const;
	void assign_to(const std::vector<double> &position, const double value, int order);
	void assign_from(const std::vector<double> &position, double* value, int order) const;
	
	// OPERATORS
	c_Type &operator[](int i){ return p_data[i]; }
	const c_Type &operator[](int i) const{ return p_data[i]; }
	c_Mesh<c_Type>& operator+=(const c_Type& rhs);
	c_Mesh<c_Type>& operator-=(const c_Type& rhs){ return *this+=-rhs; };
	c_Mesh<c_Type>& operator*=(const c_Type& rhs);
	c_Mesh<c_Type>& operator/=(const c_Type& rhs);
};

// *********************
// * CLASS :: c_Part_x *
// *********************

class c_Part_x
{
public:
	// CONSTRUCTORS
	c_Part_x(double x, double y, double z);
	
	// VARIABLES
	std::vector<double> position;
	
	// METHODS
	void get_per(int per);
	void drift(std::vector<double> velocity, double db);
	
	// OPERATORS
	double &operator[](int i){ return position[i]; }
};

// *********************
// * CLASS :: c_Part_v *
// *********************

class c_Part_v: public c_Part_x
{
public:
	// CONSTRUCTORS
	c_Part_v(double x, double y, double z);
	c_Part_v(double x, double y, double z, double vx, double vy, double vz);
	
	// VARIABLES
	std::vector<double> velocity;
	
	// METHODS
	void kick(std::vector<double> force, double db);
};

// *****************************
// * CLASS :: c_Pow_Spec_Param *
// *****************************

class c_Pow_Spec_Param
{
public:
	// VARIABLES
	double A = 1, ns, k2_G, s8;
	e_power_spec pwr_type;
	double i_pwr_type;
	
	// METHODS
	void eval_pwr();
	void norm_pwr();
	double lin_pow_spec(double k) const;
};

// ************************
// * CLASS :: c_Sim_Param *
// ************************

class c_Sim_Param
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
	c_Pow_Spec_Param power;
	
	// METHODS
	int init(int ac, char* av[]);
	void print_info();
	
protected:
	bool is_init = 0;
};

// *******************
// * CLASS :: c_Pool *
// *******************

class c_Pool
{
public:
	// CONSTRUCTORS
	c_Pool(int num_thread);
	
	// VARIABLES
	int num_thread;
	boost::threadpool::pool th_pool;
	
	// METHODS
	void add_task(int i_from, int i_to, std::function<void(int, int)> func);
};

// ***********************
// * CLASS :: c_Tracking *
// ***********************

class c_Tracking
{
public:
	// CONSTRUCTORS
	c_Tracking(int num_track_par, int par_num);
	
	// VARIABLES
	int num_track_par; // square root of number of tracking particles
	std::vector<int> par_ids;
	std::vector<std::vector<c_Part_x>> par_pos;
	
	// METHODS
	int num_step(){return par_pos.size();};
	
};

// ***************************
// * CLASS :: c_App_Var_base *
// ***************************

class c_App_Var_base
{
public:
	// CONSTRUCTORS
	c_App_Var_base(const c_Sim_Param &sim, std::string app_str);
	
	// DESTRUCTOR
	~c_App_Var_base();
	
	// VARIABLES
	int err, step = 0, print_every = 1;
	double b, b_out, db;
	const std::string z_suffix_const;
	std::vector< c_Mesh<double> > app_field;
	c_Mesh<double> power_aux;
	std::vector<fftw_complex> pwr_spec_binned, pwr_spec_binned_0;
	fftw_plan p_F, p_B;
	c_Pool pool;
	c_Tracking track;
	
	// METHODS
	double z(){ return 1./b - 1.;}
	bool integrate(){return (b <= b_out) && (db > 0);}
	bool printing(){ return ((step % print_every) == 0) or (b == b_out); }
	void upd_time();
	
	std::string z_suffix();
	
protected:	
	std::stringstream z_suffix_num;
};

// **********************
// * CLASS :: c_App_Var *
// **********************

class c_App_Var: public c_App_Var_base
{
public:
	// CONSTRUCTORS
	c_App_Var(const c_Sim_Param &sim, std::string app_str);
	
	// VARIABLES
	std::vector<c_Part_x> particles;
	
	// METHODS
	void set_unpert_pos(const c_Sim_Param &sim);
};

// ************************
// * CLASS :: c_App_Var_v *
// ************************

class c_App_Var_v: public c_App_Var_base
{
public:
	// CONSTRUCTORS
	c_App_Var_v(const c_Sim_Param &sim, std::string app_str);
	
	// VARIABLES
	std::vector<c_Part_v> particles;
	
	// METHODS
	void set_unpert_pos_w_vel(const c_Sim_Param &sim, const std::vector< c_Mesh<double> > &vel_field);
};

// ********************
// * VARIOUS ROUTINES *
// ********************

void fftw_execute_dft_r2c(const fftw_plan &p_F, c_Mesh<double>* rho, c_Pool* pool);
void fftw_execute_dft_c2r(const fftw_plan &p_B, c_Mesh<double>* rho, c_Pool* pool);
void fftw_execute_dft_r2c_triple(const fftw_plan &p_F, std::vector<c_Mesh<double>>* rho, c_Pool* pool);
void fftw_execute_dft_c2r_triple(const fftw_plan &p_B, std::vector<c_Mesh<double>>* rho, c_Pool* pool);
double wgh_sch(const std::vector<double> &x, const std::vector<int> &y, int mesh_num, int order);
// int get_per(int vec, int per); TODO, now static, have to get rid of grid_fce.cpp
namespace c11 {
void work_dir_over(std::string);
}
// *****************
// * MAIN ROUTINES *
// *****************

void get_rho_from_par(const std::vector<c_Part_v> &particles, c_Mesh<double>* rho, const int order, c_Pool* pool);
void gen_rho_dist_k(const c_Sim_Param &sim, c_Mesh<double>* rho, const fftw_plan &p_F, c_Pool* pool);
void pwr_spec_k(const c_Sim_Param &sim, const c_Mesh<double> &rho_k, c_Mesh<double>* power_aux, c_Pool* pool);
void pwr_spec(const c_Sim_Param &sim, c_Mesh<double>* rho, c_Mesh<double>* power_aux, const fftw_plan &p_F, c_Pool* pool);
void gen_pow_spec_binned(const c_Sim_Param &sim, const c_Mesh<double> &power_aux, std::vector<fftw_complex>* pwr_spec_binned);
void gen_pot_k(c_Mesh<double>* rho_k, c_Pool* pool);
void gen_displ_k(std::vector<c_Mesh<double>>* vel_field, c_Mesh<double>* pot_k, c_Pool* pool);

void print_par_pos_cut_small(const std::vector<c_Part_v> &particles, int mesh_num, int L, std::string out_dir, std::string suffix);
void print_pow_spec(const std::vector<fftw_complex> &pwr_spec_binned, std::string out_dir, std::string suffix);


// **********************
// * TEMPLATE FUNCTIONS *
// **********************

template <class c_Type>
void c_Mesh<c_Type>::assign_to(const std::vector<double> &position, const double value, int order)
{
	std::vector<int> y(3), z(3);
	for (int i = 0; i < 3; i++)
	{
		z[i] = (int)(position[i] - 0.5*(order - 1));
	}
	for (y[0] = z[0]; y[0] < z[0] + 1 + order; y[0]++)
	{
		for (y[1] = z[1]; y[1] < z[1] + 1 + order; y[1]++){
		
			for (y[2] = z[2]; y[2] < z[2] + 1 + order; y[2]++)
			{
				p_data[index(y.data())] += value * wgh_sch(position, y, N1, order);
			}
		}
	}
}

template <class c_Type>
void c_Mesh<c_Type>::assign_from(const std::vector<double> &position, double* value, int order) const
{
	std::vector<int> y(3), z(3);
	for (int i = 0; i < 3; i++)
	{
		z[i] = (int)(position[i] - 0.5*(order - 1));
	}
	for (y[0] = z[0]; y[0] < z[0] + 1 + order; y[0]++)
	{
		for (y[1] = z[1]; y[1] < z[1] + 1 + order; y[1]++){
		
			for (y[2] = z[2]; y[2] < z[2] + 1 + order; y[2]++)
			{
				*value += p_data[index(y.data())] * wgh_sch(position, y, N1, order);
			}
		}
	}
}