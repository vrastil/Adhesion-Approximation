/**
 * @file:	core_out.h
 * @brief:	functions handling output of the program
 */

#pragma once

#include "stdafx.h"
#include "core.h"
#include <fftw3.h>

void work_dir_over(std::string out_dir);

void print_pow_spec(const std::vector<fftw_complex> &pwr_spec_binned, std::string out_dir, std::string suffix);
void print_pow_spec_diff(const std::vector<fftw_complex> &pwr_spec_binned, const std::vector<fftw_complex> &pwr_spec_binned_0, 
	double b, std::string out_dir, std::string suffix);
void print_par_pos_cut_small(Particle_x* particles, const Sim_Param &sim, std::string out_dir, std::string suffix);
void print_track_par(const Tracking& track, const Sim_Param &sim, std::string out_dir, std::string suffix);
void print_rho_map(const Mesh& rho, const Sim_Param &sim, std::string out_dir, std::string suffix);
void print_suppression(const std::vector<double_2> &supp, const Sim_Param &sim, std::string out_dir);
void print_dens_bin(const std::vector<int> &dens_binned, int mesh_num, std::string out_dir, std::string suffix);