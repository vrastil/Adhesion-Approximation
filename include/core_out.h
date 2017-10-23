/**
 * @file:	core_out.h
 * @brief:	functions handling output of the program
 */

#pragma once

#include "stdafx.h"
#include "core.h"

std::string currentDateTime();
std::string std_out_dir(std::string pre_subdir, const Sim_Param &sim);

void create_dir(std::string out_dir);
void work_dir_over(std::string out_dir);

void print_pow_spec(const Data_x_y<double> &pwr_spec_binned, std::string out_dir, std::string suffix);
void print_vel_pow_spec(const Data_x_y<double> &pwr_spec_binned, std::string out_dir, std::string suffix);
void print_corr_func(const Data_x_y<double> &pwr_spec_binned, std::string out_dir, std::string suffix);
void print_pow_spec_diff(const Data_x_y<double> &pwr_spec_binned, const Data_x_y<double> &pwr_spec_binned_0, 
    double b, std::string out_dir, std::string suffix);
void print_vel_pow_spec_diff(const Data_x_y<double> &pwr_spec_binned, const Data_x_y<double> &pwr_spec_binned_0, 
    double b, std::string out_dir, std::string suffix);
template <class T>
void print_par_pos_cut_small(T* particles, const Sim_Param &sim, std::string out_dir, std::string suffix);
void print_track_par(const Tracking& track, const Sim_Param &sim, std::string out_dir, std::string suffix);
void print_rho_map(const Mesh& rho, const Sim_Param &sim, std::string out_dir, std::string suffix);
void print_projected_rho(const Mesh& delta, const Sim_Param &sim, std::string out_dir, std::string suffix);
void print_suppression(const Data_x_y<double> &supp, const Sim_Param &sim, std::string out_dir);
void print_dens_bin(const std::vector<int> &dens_binned, int mesh_num, std::string out_dir, std::string suffix);