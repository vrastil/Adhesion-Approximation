/**
 * @brief functions handling output of the program
 * 
 * @file core_out.h
 * @author Michal Vrastil
 * @date 2018-07-11
 */

#pragma once
#include "stdafx.h"
#include <fstream>
#include "params.hpp"
#include "core_power.h"
#include "templates/class_data_vec.hpp"
#include "templates/class_mesh.hpp"

/**
 * class  Ofstream handles opening and closing files, has 16MB buffer for output 
 */

class Ofstream : public std::ofstream
{
public:
    Ofstream(std::string file_name);
    char* buf;
    ~Ofstream();
};

/**
 * class  Ifstream handles opening and closing files
 */

class Ifstream : public std::ifstream
{
public:
    Ifstream(std::string file_name);
    ~Ifstream();
};

std::string currentDateTime();
std::string std_out_dir(const std::string& pre_subdir, const Sim_Param &sim);

void create_dir(const std::string& out_dir);
void remove_dir(const std::string &out_dir);
void remove_all_files(const std::string &out_dir);

void print_pow_spec(const Data_Vec<FTYPE_t, 2> &pwr_spec_binned, std::string out_dir, std::string suffix);
void print_vel_pow_spec(const Data_Vec<FTYPE_t, 2> &pwr_spec_binned, std::string out_dir, std::string suffix);
void print_corr_func(const Data_Vec<FTYPE_t, 2> &pwr_spec_binned, std::string out_dir, std::string suffix);
void print_pow_spec_diff(const Data_Vec<FTYPE_t, 2> &pwr_spec_binned, const Data_Vec<FTYPE_t, 2> &pwr_spec_binned_0,
    FTYPE_t growth, std::string out_dir, std::string suffix);
void print_pow_spec_diff(const Data_Vec<FTYPE_t, 2> &pwr_spec_binned, const Interp_obj &pwr_spec_input,
    FTYPE_t growth, std::string out_dir, std::string suffix);
void print_pow_spec_diff(const Data_Vec<FTYPE_t, 2> &pwr_spec_binned, const Data_Vec<FTYPE_t, 2> &pwr_spec_binned_0,
    const Interp_obj &pwr_spec_input, FTYPE_t growth_now, FTYPE_t growth_init, std::string out_dir, std::string suffix);
void print_vel_pow_spec_diff(const Data_Vec<FTYPE_t, 2> &pwr_spec_binned, const Data_Vec<FTYPE_t, 2> &pwr_spec_binned_0,
    FTYPE_t b, std::string out_dir, std::string suffix);

template <class T>
void print_par_pos_cut_small(const std::vector<T>& particles, const Sim_Param &sim, std::string out_dir, std::string suffix);
void print_rho_map(const Mesh& rho, const Sim_Param &sim, std::string out_dir, std::string suffix);
void print_projected_rho(const Mesh& delta, const Sim_Param &sim, std::string out_dir, std::string suffix);
void print_dens_bin(const std::vector<int> &dens_binned, std::string out_dir, std::string suffix);