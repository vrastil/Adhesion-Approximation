/**
 * @file:	app_var.hpp
 * @brief:	classes handling approximations data
 */

#pragma once
#include "stdafx.h"
#include "precision.hpp"
#include "templates/class_data_vec.hpp"
#include "templates/class_mesh.hpp"
class Sim_Param; //< params.hpp

/**
 * @class:	App_Var
 * @brief:	class containing core variables and methods for approximations
 * @temp:   <T> is a type of particle, implemented Particle_x and Particle_v
 */

template <class T> 
class App_Var
{
public:
	// CONSTRUCTORS & DESTRUCTOR
    App_Var(const Sim_Param &sim, const std::string& app_short, const std::string& app_long);
	~App_Var();

    // RUN THE SIMULATION
    void run_simulation();
	
protected:
    // VARIABLES
    const Sim_Param &sim;
    uint64_t memory_alloc; // only the largest chunks, NEED to increase in derived classes appropriately
    
    // LARGE FIELDS
	std::vector<Mesh> app_field;
    std::vector<Mesh> power_aux;
    std::vector<T> particles;

    // OTHER FIELDS
    Data_Vec<FTYPE_t, 2> corr_func_binned, pwr_spec_binned, pwr_spec_binned_0, vel_pwr_spec_binned_0;
	FFTW_PLAN_TYPE p_F, p_B, p_F_pwr, p_B_pwr;
	std::vector<int> dens_binned;
	
	// METHODS
    FTYPE_t a();
	FTYPE_t a_half();
    FTYPE_t da();
    std::string get_out_dir() const;
    std::string get_z_suffix() const;

    // METHODS THAT NEED TO BE OVERRIDEN IN DERIVED CLASSES
    virtual void print_output(); //< save info about simulation state
private:
    virtual void pot_corr(); //< CIC correction by default
    virtual void upd_pos() = 0;

    // IMPLEMENTATION
    class Impl;
    friend class Impl;
    const std::unique_ptr<Impl> m_impl;
};