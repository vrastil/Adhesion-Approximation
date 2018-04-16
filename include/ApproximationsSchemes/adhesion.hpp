#pragma once

#include "stdafx.h"
#include "app_var.hpp"
#include "precision.hpp"
#include "templates/class_particles.hpp"

class Sim_Param; //< declaration / definition in 'params.x'

/**
 * @class:	App_Var_AA
 * @brief:	class containing variables and methods for adhesion approximation
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