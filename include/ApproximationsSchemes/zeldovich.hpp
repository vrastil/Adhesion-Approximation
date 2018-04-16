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

class App_Var_ZA: public App_Var<Particle_v<FTYPE_t>>
{
public:
	// CONSTRUCTORS & DESTRUCTOR
	App_Var_ZA(const Sim_Param &sim);

private:
    // no CIC correction for ZA
    void pot_corr() override;

    // ZA with velocitites
    void upd_pos() override;
};