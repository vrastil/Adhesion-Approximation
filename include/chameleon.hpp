#pragma once

#include "stdafx.h"
#include "params.hpp"
#include "app_var.hpp"

/**
 * @class:	App_Var_chi
 * @brief:	class containing variables for chameleon gravity
 */
 
class App_Var_chi: public App_Var<Particle_v<FTYPE_t>>
{
public:
	// CONSTRUCTORS & DESTRUCTOR
	App_Var_chi(const Sim_Param &sim, std::string app_str);
    ~App_Var_chi();

    // METHODS
    void print_output();
    void solve(FTYPE_t a);

private:
    // IMPLEMENTATION
    class ChiImpl;
    std::unique_ptr<ChiImpl> m_impl;
};