/**
 * @brief functions for integration of particle trajectories
 * 
 * @file integration.hpp
 * @author Michal Vrastil
 * @date 2018-07-11
 */

#pragma once
#include "stdafx.h"
#include <functional>
#include "precision.hpp"
#include "class_mesh.hpp"
#include "class_particles.hpp"

class Cosmo_Param;

void stream_step(const FTYPE_t da, std::vector<Particle_v<FTYPE_t>>& particles);
void stream_kick_stream(const FTYPE_t da, std::vector<Particle_v<FTYPE_t>>& particles, std::function<void()> kick_step, size_t per);
void kick_step_no_momentum(const Cosmo_Param &cosmo, const FTYPE_t a, std::vector<Particle_v<FTYPE_t>>& particles, const std::vector< Mesh> &vel_field);
void kick_step_w_momentum(const Cosmo_Param &cosmo, const FTYPE_t a, const FTYPE_t da, std::vector<Particle_v<FTYPE_t>>& particles, const std::vector< Mesh> &force_field);
