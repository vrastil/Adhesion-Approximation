/**
 * @brief Frozen-potential approximation interface
 * 
 * @file frozen_potential.hpp
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
 * @class:	App_Var_FP
 * @brief:	class containing variables and methods for Frozen-potential approximation
 * @ingroup APP
 */
class App_Var_FP: public App_Var<Particle_v<FTYPE_t>>
{
public:
	// CONSTRUCTORS & DESTRUCTOR
	App_Var_FP(const Sim_Param &sim);

private:
    // Leapfrog method for frozen-potential
    void upd_pos() override;
};