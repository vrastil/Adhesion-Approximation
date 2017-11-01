/**
 * @file:	core_app.h
 * @brief:	common functions for all types of approximations
 */
 
#pragma once

#include "stdafx.h"
#include "core.h"
#include "core_power.h"

void set_unpert_pos(const Sim_Param &sim, Particle_x* particles);
void set_unpert_pos_w_vel(const Sim_Param &sim, Particle_v* particles, const std::vector< Mesh> &vel_field);
void set_pert_pos(const Sim_Param &sim, const double db, Particle_x* particles, const std::vector< Mesh> &vel_field);
void set_pert_pos_w_vel(const Sim_Param &sim, const double db, Particle_v* particles, const std::vector< Mesh> &vel_field);

void upd_pos_first_order(const Sim_Param &sim, const double da, const double a, Particle_v* particles, const std::vector< Mesh> &vel_field);
void upd_pos_second_order(const Sim_Param &sim, const double db, const double b, Particle_v* particles, const std::vector< Mesh> &force_field);
void upd_pos_second_order_w_short_force(const Sim_Param &sim, LinkedList* linked_list, const double db, const double b, 
	                                    Particle_v* particles, const std::vector< Mesh> &force_field);

void gen_rho_dist_k(const Sim_Param &sim, Mesh* rho, const fftw_plan &p_F);
void gen_pot_k(const Mesh& rho_k, Mesh* pot_k);
void gen_pot_k(Mesh* rho_k);
void gen_displ_k(std::vector<Mesh>* vel_field, const Mesh& pot_k);
void gen_displ_k_cic(std::vector<Mesh>* vel_field, const Mesh& pot_k);
void gen_displ_k_S2(std::vector<Mesh>* vel_field, const Mesh& pot_k, const double a);

template <class T>
void get_rho_from_par(T* particles, Mesh* rho, const Sim_Param &sim);
int get_vel_from_par(Particle_v* particles, std::vector<Mesh>* vel_field, const Sim_Param &sim);
int get_vel_from_par(Particle_x* particles, std::vector<Mesh>* vel_field, const Sim_Param &sim);

void pwr_spec_k(const Sim_Param &sim, const Mesh &rho_k, Mesh* power_aux);
void vel_pwr_spec_k(const Sim_Param &sim, const std::vector<Mesh> &vel_field, Mesh* power_aux);
void gen_pow_spec_binned(const Sim_Param &sim, const Mesh &power_aux, Data_x_y<double>* pwr_spec_binned);
void gen_pow_spec_binned_from_extrap(const Sim_Param &sim, const Extrap_Pk &P_k, Data_x_y<double>* pwr_spec_binned);
void gen_corr_func_binned(const Sim_Param &sim, const Mesh &power_aux, Data_x_y<double>* corr_func_binned);
template<class T>
void gen_corr_func_binned_pp(const Sim_Param &sim, T* particles, Data_x_y<double>* corr_func_binned);
void gen_dens_binned(const Mesh& rho, std::vector<int> &dens_binned, const Sim_Param &sim);

void force_short(const Sim_Param &sim, const double D, const LinkedList& linked_list, Particle_v *particles,
				 const Vec_3D<double> position, Vec_3D<double>* force);