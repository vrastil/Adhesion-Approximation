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

App_Var_FP::App_Var_FP(const Sim_Param &sim, const std::string& app_short, const std::string& app_long):
    App_Var<Particle_v<FTYPE_t>>(sim, app_short, app_long) {}

void App_Var_FP::upd_pos()
{// Leapfrog method for frozen-potential
    auto kick_step = [&](){ kick_step_w_momentum(sim.cosmo, a_half(), da(), particles, app_field); };
    stream_kick_stream(da(), particles, kick_step, sim.box_opt.mesh_num);
}

App_Var_PM::App_Var_PM(const Sim_Param &sim):
    App_Var_FP(sim, "PM", "Particle-mesh approximation") {}

void App_Var_PM::upd_pos()
{// Leapfrog method for frozen-potential
    auto kick_step = [&]()
    {
        // Get discrete density from particles
        get_rho_from_par(particles, app_field[0], sim);

        // get rho_k
        fftw_execute_dft_r2c(p_F, app_field[0]);

        // get potential
        gen_pot_k(app_field[0], app_field[0]);

        // CIC correction of potential, get forces in real space
        pot_corr(app_field, app_field[0]);

        // apply forces
        kick_step_w_momentum_pm(sim.cosmo, a_half(), da(), particles, app_field);
    };
    stream_kick_stream(da(), particles, kick_step, sim.box_opt.mesh_num);
}