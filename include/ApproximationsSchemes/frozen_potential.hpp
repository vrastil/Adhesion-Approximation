#pragma once

#include "stdafx.h"
#include "app_var.hpp"
#include "precision.hpp"
#include "templates/class_particles.hpp"

class Sim_Param; //< declaration / definition in 'params.x'

/**
 * @class:	App_Var_FF
 * @brief:	class containing variables and methods for Frozen-flow approximation
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