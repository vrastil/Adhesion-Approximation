/**
 * @brief frozen-flow approximation interface
 * 
 * @file frozen_flow.hpp
 * @author Michal Vrastil
 * @date 2018-07-11
 */

#pragma once

#include "stdafx.h"
#include "app_var.hpp"
#include "precision.hpp"
#include "templates/class_particles.hpp"

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
 * @class:	App_Var_FF
 * @brief:	class containing variables and methods for Frozen-flow approximation
 * @ingroup APP
 */
class App_Var_FF: public App_Var<Particle_v<FTYPE_t>>
{
public:
	// CONSTRUCTORS & DESTRUCTOR
	App_Var_FF(const Sim_Param &sim);

private:
    // Leapfrog method for frozen-flow
    void upd_pos() override;
};