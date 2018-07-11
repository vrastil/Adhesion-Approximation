/**
 * @brief chameleon model of gravity interface
 * 
 * @file chameleon.hpp
 * @author Michal Vrastil
 * @date 2018-07-08
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
 * @class:	App_Var_Chi
 * @brief:	class containing variables and methods for chameleon gravity
 * @ingroup APP
 */
class App_Var_Chi: public App_Var<Particle_v<FTYPE_t>>
{
public:
	// CONSTRUCTORS & DESTRUCTOR
	App_Var_Chi(const Sim_Param &sim);
    ~App_Var_Chi();
private:
    // IMPLEMENTATION
    class ChiImpl;
    const std::unique_ptr<ChiImpl> m_impl;

    // Leapfrog method for chameleon gravity (frozen-potential)
    void upd_pos() override;

    // Print additional information about chameleon field
    void print_output() override;
};