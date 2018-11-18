/**
 * @brief ahesion approximation interface
 * 
 * @file adhesion.hpp
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
 * @class:	App_Var_AA
 * @brief:	class containing variables and methods for adhesion approximation
 * @ingroup APP
 */
class App_Var_AA: public App_Var<Particle_v<FTYPE_t>>
{
public:
	// CONSTRUCTORS & DESTRUCTOR
	App_Var_AA(const Sim_Param &sim);
    ~App_Var_AA();

private:
    // IMPLEMENTATION
    class AAImpl;
    const std::unique_ptr<AAImpl> m_impl;

    // initialize potentials for adhesion
    void pot_corr() override;

    // Leapfrog method for adhesion
    void upd_pos() override;
};