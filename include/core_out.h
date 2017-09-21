/**
 * @file:	core_out.h
 * @brief:	functions handling output of the program
 */

#pragma once

#include "stdafx.h"
#include "core.h"
#include <fftw3.h>

std::string currentDateTime();
std::string std_out_dir(std::string pre_subdir, const Sim_Param &sim);

void create_dir(std::string out_dir);
void work_dir_over(std::string out_dir);

void print_pow_spec(const std::vector<double_2> &pwr_spec_binned, std::string out_dir, std::string suffix);
void print_pow_spec_diff(const std::vector<double_2> &pwr_spec_binned, const std::vector<double_2> &pwr_spec_binned_0, 
    double b, std::string out_dir, std::string suffix);
template <class T>
void print_par_pos_cut_small(T* particles, const Sim_Param &sim, std::string out_dir, std::string suffix);
void print_track_par(const Tracking& track, const Sim_Param &sim, std::string out_dir, std::string suffix);
void print_rho_map(const Mesh& rho, const Sim_Param &sim, std::string out_dir, std::string suffix);
void print_projected_rho(const Mesh& delta, const Sim_Param &sim, std::string out_dir, std::string suffix);
void print_suppression(const std::vector<double_2> &supp, const Sim_Param &sim, std::string out_dir);
void print_dens_bin(const std::vector<int> &dens_binned, int mesh_num, std::string out_dir, std::string suffix);


/**
 * @brief:	template class functions definitions
 */

 template <class T>
 void print_par_pos_cut_small(T* particles, const Sim_Param &sim, std::string out_dir, std::string suffix)
 {
	out_dir += "par_cut/";
	std::string file_name = out_dir + "par_cut" + suffix + ".dat";
	FILE* pFile;
	pFile = fopen(file_name.c_str(), "w");
	if (pFile == NULL)
	{
		printf("Error while opening %s\n", file_name.c_str());
		perror("Error");
		return;
	}
	
	std::cout << "Writing small cut through the box of particles into file " << out_dir + "par_cut" + suffix + ".dat\n";
	fprintf (pFile, "# This file contains positions of particles in units [Mpc/h].\n");
    double x, y, z, dx;
    double x_0 = sim.x_0();
	for(int i=0; i < sim.par_num; i++)
	{
		x = particles[i].position.x;
		y = particles[i].position.y;
		z = particles[i].position.z;			
		dx = abs(y - sim.mesh_num/2.);
		if ((dx < 0.5) && (x < sim.mesh_num/4.) && (z < sim.mesh_num/4.))
		{
			// cut (L/4 x L/4 x 0.5)
			fprintf (pFile, "%f\t%f\t%f\n", x*x_0 , z*x_0, y*x_0);
		}
	}
	fclose (pFile);
 }