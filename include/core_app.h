/**
 * @file:	core_app.h
 * @brief:	common functions for all types of approximations
 */
 
#pragma once
#include "core.h"
#include "core_power.h"

void set_unpert_pos(const Sim_Param &sim, Particle_x* particles);
void set_unpert_pos_w_vel(const Sim_Param &sim, Particle_v* particles, const std::vector< Mesh> &vel_field);
void set_pert_pos(const Sim_Param &sim, const FTYPE db, Particle_x* particles, const std::vector< Mesh> &vel_field);
void set_pert_pos_w_vel(const Sim_Param &sim, const FTYPE db, Particle_v* particles, const std::vector< Mesh> &vel_field);

void upd_pos_first_order(const Sim_Param &sim, const FTYPE da, const FTYPE a, Particle_v* particles, const std::vector< Mesh> &vel_field);
void upd_pos_second_order(const Sim_Param &sim, const FTYPE da, const FTYPE a, Particle_v* particles, const std::vector< Mesh> &force_field);
void upd_pos_second_order_w_pp(const Sim_Param &sim, const FTYPE da, const FTYPE a, Particle_v* particles, const std::vector< Mesh> &force_field,
                               LinkedList* linked_list, Interp_obj* fs_interp);

void gen_rho_dist_k(const Sim_Param &sim, Mesh* rho, const fftw_plan &p_F);
void gen_pot_k(const Mesh& rho_k, Mesh* pot_k);
void gen_pot_k(Mesh* rho_k);
void gen_displ_k(std::vector<Mesh>* vel_field, const Mesh& pot_k);
void gen_displ_k_cic(std::vector<Mesh>* vel_field, const Mesh& pot_k);
void gen_displ_k_S2(std::vector<Mesh>* vel_field, const Mesh& pot_k, const FTYPE a);

template <class T>
void get_rho_from_par(T* particles, Mesh* rho, const Sim_Param &sim);
int get_vel_from_par(Particle_v* particles, std::vector<Mesh>* vel_field, const Sim_Param &sim);
int get_vel_from_par(Particle_x* particles, std::vector<Mesh>* vel_field, const Sim_Param &sim);

void pwr_spec_k(const Mesh &rho_k, Mesh* power_aux);
void pwr_spec_k_init(const Mesh &rho_k, Mesh* power_aux);
void vel_pwr_spec_k(const std::vector<Mesh> &vel_field, Mesh* power_aux);
void gen_pow_spec_binned(const Sim_Param &sim, const Mesh &power_aux, Data_Vec<FTYPE, 2>* pwr_spec_binned);
void gen_pow_spec_binned_init(const Sim_Param &sim, const Mesh &power_aux, const unsigned half_length, Data_Vec<FTYPE, 2>* pwr_spec_binned);
void gen_pow_spec_binned_from_extrap(const Sim_Param &sim, const Extrap_Pk &P_k, Data_Vec<FTYPE, 2>* pwr_spec_binned);
void gen_dens_binned(const Mesh& rho, std::vector<int> &dens_binned, const Sim_Param &sim);

FTYPE force_ref(const FTYPE r, const FTYPE a);
FTYPE force_tot(const FTYPE r, const FTYPE e2);
void force_short(const Sim_Param &sim, const FTYPE D, const LinkedList& linked_list, Particle_v *particles,
				 const Vec_3D<FTYPE> position, Vec_3D<FTYPE>* force, Interp_obj* fs_interp);