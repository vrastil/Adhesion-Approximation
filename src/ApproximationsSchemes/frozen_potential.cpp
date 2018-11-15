/**
 * @brief frozen-potential approximation implementation
 * 
 * @file frozen_potential.cpp
 * @author Michal Vrastil
 * @date 2018-07-11
 */

#include "frozen_potential.hpp"
#include "core_app.h"
#include "core_mesh.h"
#include "integration.hpp"
#include "params.hpp"

App_Var_FP::App_Var_FP(const Sim_Param &sim):
    App_Var<Particle_v<FTYPE_t>>(sim, "FP", "Frozen-potential approximation") {}

void App_Var_FP::upd_pos()
{// Leapfrog method for frozen-potential
    auto kick_step = [&](){ kick_step_w_momentum(sim.cosmo, a_half(), da(), particles, app_field); };
    stream_kick_stream(da(), particles, kick_step, sim.box_opt.mesh_num);
}