#pragma once

#include "stdafx.h"
#include "app_var.hpp"
#include "precision.hpp"
#include "templates/class_particles.hpp"

class Sim_Param; //< declaration / definition in 'params.x'

/**
 * @class:	App_Var_FP_mod
 * @brief:	class containing variables and methods for modified Frozen-potential approximation
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
    void pot_corr() override;

    // Leapfrog method for modified frozen-potential
    void upd_pos() override;
};