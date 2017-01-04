/**
 * @file:	core_app.h
 * @brief:	common functions for all types of approximations
 */
 
#include "stdafx.h"
#include "core.h"

void set_unpert_pos(const Sim_Param &sim, Particle_x* particles);
void set_unpert_pos_w_vel(const Sim_Param &sim, Particle_v* particles, const std::vector< Mesh> &vel_field);