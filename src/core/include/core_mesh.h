/**
 * @brief basic functions to work with mesh
 * 
 * @file core_mesh.h
 * @author Michal Vrastil
 * @date 2018-06-24
 */

#pragma once
#include "stdafx.h"
#include "params.hpp"
#include "precision.hpp"
#include "class_mesh.hpp"
#include "class_particles.hpp"
#include "class_vec_3d.hpp"

void get_k_vec(size_t N, size_t index, int* k_vec);
void get_k_vec(size_t N, size_t index, Vec_3D<int> &k_vec);
FTYPE_t get_k_sq(size_t N, size_t index);

template<typename T> void get_per(Vec_3D<T> &position, size_t per);
template<typename T> void get_per(Vec_3D<T> &position, size_t perx, size_t pery, size_t perz);
void get_per(std::vector<Particle_v<FTYPE_t>>& particles, const size_t per);

FTYPE_t get_distance(const Vec_3D<FTYPE_t> &x_1, const Vec_3D<FTYPE_t> &x_2, size_t per);
Vec_3D<FTYPE_t> get_sgn_distance(const Vec_3D<FTYPE_t> &x_from, const Vec_3D<FTYPE_t> &x_to, size_t per);

void assign_to(Mesh& field, const Vec_3D<FTYPE_t> &position, const FTYPE_t value);
void assign_to(std::vector<Mesh>& field, const Vec_3D<FTYPE_t> &position, const Vec_3D<FTYPE_t>& value);
void assign_from(const Mesh &field, const Vec_3D<FTYPE_t> &position, FTYPE_t& value, FTYPE_t mod = 1);
void assign_from(const std::vector<Mesh> &field, const Vec_3D<FTYPE_t> &position, Vec_3D<FTYPE_t>& value, FTYPE_t mod = 1);

/**
 * @brief compute forward (real to complex) FFT on mesh (inplace)
 * 
 * @param p_F plan for forward transformation
 * @param rho mesh upon which the transformation is performed
 */
void fftw_execute_dft_r2c(const FFTW_PLAN_TYPE &p_F, Mesh& rho);

/**
 * @brief compute backward (complex to real) FFT on mesh (inplace)
 * 
 * @param p_B plan for backward transformation
 * @param rho mesh upon which the transformation is performed
 */
void fftw_execute_dft_c2r(const FFTW_PLAN_TYPE &p_B, Mesh& rho);

/**
 * @brief compute three forward (real to complex) FFTs on vector of meshes (inplace)
 * 
 * @param p_F plan for forward transformations
 * @param rho vector of meshes upon which the transformations are performed
 */
void fftw_execute_dft_r2c_triple(const FFTW_PLAN_TYPE &p_F, std::vector<Mesh>& rho);

/**
 * @brief compute three backward (complex to real) FFTs on vector of meshes (inplace)
 * 
 * @param p_B plan for backward transformations
 * @param rho vector of meshes upon which the transformations are performed
 */
void fftw_execute_dft_c2r_triple(const FFTW_PLAN_TYPE &p_B, std::vector<Mesh>& rho);

template<unsigned int points>
class IT
{
public:
    IT(const Vec_3D<FTYPE_t> &pos); // ctor for assignment scheme
    IT(const Vec_3D<FTYPE_t> &pos, FTYPE_t Hc); // ctor for chaining mesh

    size_t counter;
    Vec_3D<size_t> vec;

    bool iter();
};

