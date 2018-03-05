/**
 * @file:	core_app.h
 * @brief:	common functions for all types of approximations
 */
 
#pragma once
#include "stdafx.h"
#include "app_var.hpp"
#include "core_power.h"
#include "chameleon.hpp"
#include "params.hpp"
#include "precision.hpp"

void set_unpert_pos(const Sim_Param &sim, std::vector<Particle_x<FTYPE_t>>& particles);
void set_unpert_pos_w_vel(const Sim_Param &sim, std::vector<Particle_v<FTYPE_t>>& particles, const std::vector< Mesh> &vel_field);
void set_pert_pos(const Sim_Param &sim, const FTYPE_t db, std::vector<Particle_x<FTYPE_t>>& particles, const std::vector< Mesh> &vel_field);
void set_pert_pos_w_vel(const Sim_Param &sim, const FTYPE_t db, std::vector<Particle_v<FTYPE_t>>& particles, const std::vector< Mesh> &vel_field);

void upd_pos_first_order(const Sim_Param &sim, const FTYPE_t da, const FTYPE_t a, std::vector<Particle_v<FTYPE_t>>& particles, const std::vector< Mesh> &vel_field);
void upd_pos_second_order(const Sim_Param &sim, const FTYPE_t da, const FTYPE_t a, std::vector<Particle_v<FTYPE_t>>& particles, const std::vector< Mesh> &force_field);
void upd_pos_second_order_w_pp(const Sim_Param &sim, const FTYPE_t da, const FTYPE_t a, std::vector<Particle_v<FTYPE_t>>& particles, const std::vector< Mesh> &force_field,
                               LinkedList& linked_list, Interp_obj& fs_interp);
void upd_pos_second_order_w_chi(const Sim_Param &sim, const FTYPE_t da, const FTYPE_t a, std::vector<Particle_v<FTYPE_t>>& particles, const std::vector< Mesh> &force_field,
                                ChiSolver<CHI_PREC_t>& sol);

void gen_rho_dist_k(const Sim_Param &sim, Mesh& rho, const FFTW_PLAN_TYPE &p_F);
void gen_pot_k(const Mesh& rho_k, Mesh& pot_k);
void gen_pot_k(Mesh& rho_k);
void gen_displ_k(std::vector<Mesh>& vel_field, const Mesh& pot_k);
void gen_displ_k_cic(std::vector<Mesh>& vel_field, const Mesh& pot_k);
void gen_displ_k_S2(std::vector<Mesh>& vel_field, const Mesh& pot_k, const FTYPE_t a);

template <class T>
void get_rho_from_par(const std::vector<T>& particles, Mesh& rho, const Sim_Param &sim);
int get_vel_from_par(const std::vector<Particle_v<FTYPE_t>>& particles, std::vector<Mesh>& vel_field, const Sim_Param &sim);
int get_vel_from_par(const std::vector<Particle_x<FTYPE_t>>& particles, std::vector<Mesh>& vel_field, const Sim_Param &sim);

void pwr_spec_k(const Mesh &rho_k, Mesh& power_aux);
void pwr_spec_k_init(const Mesh &rho_k, Mesh& power_aux);
void vel_pwr_spec_k(const std::vector<Mesh> &vel_field, Mesh& power_aux);
void gen_pow_spec_binned(const Sim_Param &sim, const Mesh &power_aux, Data_Vec<FTYPE_t, 2>& pwr_spec_binned);
void gen_pow_spec_binned_init(const Sim_Param &sim, const Mesh &power_aux, const unsigned half_length, Data_Vec<FTYPE_t, 2>& pwr_spec_binned);
template<class P, typename T, unsigned N> // P = everything callable P_k(k), T = float-type, N = number
void gen_pow_spec_binned_from_extrap(const Sim_Param &sim, const P &P_k, Data_Vec<T, N>& pwr_spec_binned);
void gen_dens_binned(const Mesh& rho, std::vector<int> &dens_binned, const Sim_Param &sim);

FTYPE_t force_ref(const FTYPE_t r, const FTYPE_t a);
FTYPE_t force_tot(const FTYPE_t r, const FTYPE_t e2);
void force_short(const Sim_Param &sim, const FTYPE_t D, const LinkedList& linked_list, const std::vector<Particle_v<FTYPE_t>>& particles,
				 const Vec_3D<FTYPE_t> position, Vec_3D<FTYPE_t>& force, Interp_obj* fs_interp);