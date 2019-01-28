/**
 * @brief Zel`dovich approximation interface
 * 
 * @file zeldovich.hpp
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
 * @class:	App_Var_ZA
 * @brief:	class containing variables and methods for Zel`dovich approximation
 * @ingroup APP
 */
class App_Var_ZA: public App_Var<Particle_v<FTYPE_t>>
{
public:
	// CONSTRUCTORS & DESTRUCTOR
	App_Var_ZA(const Sim_Param &sim);

protected:
    // for TZA
    App_Var_ZA(const Sim_Param &sim, const std::string& app_short, const std::string& app_long);

private:
    // no CIC correction for ZA
    void pot_corr() override;

    // ZA with velocitites
    void upd_pos() override;
};

class App_Var_TZA: public App_Var_ZA
{
public:
	// CONSTRUCTORS & DESTRUCTOR
	App_Var_TZA(const Sim_Param &sim);

    // truncation of initial power spectrum for TZA
    void update_cosmo() override;
};