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
void print_par_pos_cut_small(Particle_x* particles, const Sim_Param &sim, std::string out_dir, std::string suffix);
void print_track_par(const Tracking& track, const Sim_Param &sim, std::string out_dir, std::string suffix);