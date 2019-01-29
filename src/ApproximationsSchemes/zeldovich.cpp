/**
 * @brief Zel`dovich approximation implementation
 * 
 * @file zeldovich.cpp
 * @author Michal Vrastil
 * @date 2018-07-11
 */
#include "zeldovich.hpp"
#include "core_app.h"
#include "core_mesh.h"
#include "params.hpp"

App_Var_ZA::App_Var_ZA(const Sim_Param &sim):
    App_Var_ZA(sim, "ZA", "Zel`dovich approximation") {}

App_Var_ZA::App_Var_ZA(const Sim_Param &sim, const std::string& app_short, const std::string& app_long):
    App_Var<Particle_v<FTYPE_t>>(sim, app_short, app_long) {}

void App_Var_ZA::upd_pos()
{// ZA with velocitites
    set_pert_pos(sim, a(), particles, app_field);
}

void App_Var_ZA::pot_corr()
{
    return;
}

App_Var_TZA::App_Var_TZA(const Sim_Param &sim):
    App_Var_ZA(sim, "TZA", "Truncated Zel`dovich approximation") {}

void App_Var_TZA::update_cosmo(Cosmo_Param& cosmo)
{
    cosmo.truncated_pk = true;
}