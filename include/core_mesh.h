/**
 * @file:	core_fce.h
 * @brief:	basic functions to work with mesh
 */

#pragma once
#include "stdafx.h"
#include "params.hpp"
#include "precision.hpp"
#include "templates/class_mesh.hpp"
#include "templates/class_particles.hpp"
#include "templates/class_vec_3d.hpp"

void get_k_vec(int N, int index, int* k_vec);
void get_k_vec(int N, int index, Vec_3D<int> &k_vec);
FTYPE_t get_k_sq(int N, int index);

template<typename T> void get_per(Vec_3D<T> &position, int per);
template<typename T> void get_per(Vec_3D<T> &position, int perx, int pery, int perz);
void get_per(std::vector<Particle_v<FTYPE_t>>& particles, const unsigned par_num, const int per);

FTYPE_t get_distance(const Vec_3D<FTYPE_t> &x_1, const Vec_3D<FTYPE_t> &x_2, int per);
Vec_3D<FTYPE_t> get_sgn_distance(const Vec_3D<FTYPE_t> &x_from, const Vec_3D<FTYPE_t> &x_to, int per);

void assign_to(Mesh& field, const Vec_3D<FTYPE_t> &position, const FTYPE_t value);
void assign_to(std::vector<Mesh>& field, const Vec_3D<FTYPE_t> &position, const Vec_3D<FTYPE_t>& value);
void assign_from(const Mesh &field, const Vec_3D<FTYPE_t> &position, FTYPE_t& value);
void assign_from(const std::vector<Mesh> &field, const Vec_3D<FTYPE_t> &position, Vec_3D<FTYPE_t>& value);

void fftw_execute_dft_r2c(const FFTW_PLAN_TYPE &p_F, Mesh& rho);
void fftw_execute_dft_c2r(const FFTW_PLAN_TYPE &p_B, Mesh& rho);
void fftw_execute_dft_r2c_triple(const FFTW_PLAN_TYPE &p_F, std::vector<Mesh>& rho);
void fftw_execute_dft_c2r_triple(const FFTW_PLAN_TYPE &p_B, std::vector<Mesh>& rho);

template<int points>
class IT
{
public:
    IT(const Vec_3D<FTYPE_t> &pos); // ctor for assignment scheme
    IT(const Vec_3D<FTYPE_t> &pos, FTYPE_t Hc); // ctor for chaining mesh

    unsigned counter;
    Vec_3D<int> vec;

    bool iter();
};

