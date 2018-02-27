#pragma once
#include "params.hpp"

// approximations in standard gravity

void zel_app(const Sim_Param &sim);
void frozen_flow(const Sim_Param &sim);
void frozen_potential(const Sim_Param &sim);
void mod_frozen_potential(const Sim_Param &sim);
void adhesion_approximation(const Sim_Param &sim);
void adhesion_approximation(const Sim_Param &sim);

// modified gravities

void chameleon_gravity(const Sim_Param &sim);