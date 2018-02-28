/**
 * @file:	mesh.hpp
 * @brief:	mesh data structure definitions
 */

#pragma once
#include "stdafx.h"
#include "core_power.h"
#include "templates/class_data_vec.hpp"
#include "templates/class_mesh.hpp"
#include "templates/class_particles.hpp"
#include "templates/class_vec_3d.hpp"

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
	const int sqr_num_track_par, num_track_par; // square root of number of tracking particles
	std::vector<int> par_ids;
	std::vector<std::vector<Particle_x<FTYPE>>> par_pos;
	
	// METHODS
	const unsigned num_step() const{ return par_pos.size(); };
	template <class T>  void update_track_par(const std::vector<T>& particles);
    void print_track_par(const Sim_Param &sim, std::string out_dir, std::string suffix);
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

private:
    void print_position();
    void print_density();
    void get_binned_power_spec();
    void print_power_spec();
    void print_extrap_pwr(const Extrap_Pk<FTYPE, 2>& P_k);
    void print_corr(const Extrap_Pk<FTYPE, 2>& P_k);
    void print_vel_pwr();
};

/**
 * @class:	App_Var_AA
 * @brief:	class containing variables for adhesion approximation
 */
 
 class App_Var_AA: public App_Var<Particle_v<FTYPE>>
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
	void get_linked_list(const std::vector<Particle_v<FTYPE>>& particles);
};

/**
 * @class:	App_Var_FP_mod
 * @brief:	class containing variables for modified Frozen-potential approximation
 */

class App_Var_FP_mod: public App_Var<Particle_v<FTYPE>>
{
public:
	// CONSTRUCTORS & DESTRUCTOR
	App_Var_FP_mod(const Sim_Param &sim, std::string app_str);
	
	// VARIABLES
    LinkedList linked_list;
    Interp_obj fs_interp;
};