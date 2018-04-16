#include "ApproximationsSchemes/zeldovich.hpp"
#include "core_app.h"
#include "core_mesh.h"
#include "params.hpp"

App_Var_ZA::App_Var_ZA(const Sim_Param &sim):
    App_Var<Particle_v<FTYPE_t>>(sim, "ZA", "Zel`dovich approximation") {}

void App_Var_ZA::upd_pos()
{// ZA with velocitites
    set_pert_pos(sim, a(), particles, app_field);
}

void App_Var_ZA::pot_corr()
{
    return;
}