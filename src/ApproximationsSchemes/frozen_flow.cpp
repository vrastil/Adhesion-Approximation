/**
 * @brief frozen-flow approximation implementation
 * 
 * @file frozen_flow.cpp
 * @author Michal Vrastil
 * @date 2018-07-11
 */

#include "frozen_flow.hpp"
#include "core_app.h"
#include "core_mesh.h"
#include "integration.hpp"
#include "params.hpp"

App_Var_FF::App_Var_FF(const Sim_Param &sim):
    App_Var<Particle_v<FTYPE_t>>(sim, "FF", "Frozen-flow approximation") {}

void App_Var_FF::upd_pos()
{// Leapfrog method for frozen-flow
    auto kick_step = [&](){ kick_step_no_momentum(sim.cosmo, a_half(), particles, app_field); };
    stream_kick_stream(da(), particles, kick_step, sim.box_opt.mesh_num);
}