/**
 * @brief modified frozen-potential approximation interface
 * 
 * @file mod_frozen_potential.hpp
 * @author Michal Vrastil
 * @date 2018-07-11
 */
#pragma once

#include "stdafx.h"
#include "app_var.hpp"
#include "precision.hpp"
#include "class_particles.hpp"

/********************//**
 * FORWARD DECLARATIONS *
 ************************/

class Sim_Param;

/**************//**
 * PUBLIC METHODS *
 ******************/

/**************//**
 * PUBLIC CLASSES *
 ******************/

/**
 * @class:	App_Var_FP_mod
 * @brief:	class containing variables and methods for modified Frozen-potential approximation
 * @ingroup APP
 */
class App_Var_FP_mod: public App_Var<Particle_v<FTYPE_t>>
{
public:
	// CONSTRUCTORS & DESTRUCTOR
	App_Var_FP_mod(const Sim_Param &sim);
    ~App_Var_FP_mod();

private:
    // IMPLEMENTATION
    class FP_ppImpl;
    const std::unique_ptr<FP_ppImpl> m_impl;

    // force interpolation corrections, long range potential for S2-shaped particles
    void pot_corr(std::vector<Mesh>& vel_field, Mesh& pot_k) override;

    // Leapfrog method for modified frozen-potential
    void upd_pos() override;
};