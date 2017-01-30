/**
 * @file:	core_app.h
 * @brief:	common functions for all types of approximations
 */
 
#pragma once

#include "stdafx.h"
#include "core.h"
#include <fftw3.h>

void set_unpert_pos(const Sim_Param &sim, Particle_x* particles);
void set_unpert_pos_w_vel(const Sim_Param &sim, Particle_v* particles, const std::vector< Mesh> &vel_field);
void set_pert_pos(const Sim_Param &sim, double db, Particle_x* particles, const std::vector< Mesh> &vel_field);

void upd_pos_first_order(const Sim_Param &sim, double db, Particle_x* particles, const std::vector< Mesh> &vel_field);
// void upd_pos_first_order(const Sim_Param &sim, double db, Particle_v* particles, const std::vector< Mesh> &vel_field);
void upd_pos_second_order(const Sim_Param &sim, double db, double b, Particle_v* particles, const std::vector< Mesh> &force_field);

void gen_rho_dist_k(const Sim_Param &sim, Mesh* rho, const fftw_plan &p_F);
void gen_pot_k(Mesh* rho_k);
void gen_displ_k(std::vector<Mesh>* vel_field, const Mesh& pot_k);

void get_rho_from_par(Particle_x* particles, Mesh* rho, const Sim_Param &sim);
void get_rho_from_par(Particle_v* particles, Mesh* rho, const Sim_Param &sim);
void pwr_spec_k(const Sim_Param &sim, const Mesh &rho_k, Mesh* power_aux);
void gen_pow_spec_binned(const Sim_Param &sim, const Mesh &power_aux, std::vector<double_2>* pwr_spec_binned);
void gen_dens_binned(const Mesh& rho, std::vector<int> &dens_binned, const Sim_Param &sim);
