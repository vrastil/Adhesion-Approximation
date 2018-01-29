/**
 * @file:	core_fce.h
 * @brief:	basic functions to work with mesh
 */

#pragma once
#include "core.h"

void get_k_vec(int N, int index, int* k_vec);
void get_k_vec(int N, int index, Vec_3D<int> &k_vec);
FTYPE get_k_sq(int N, int index);

template<typename T> void get_per(Vec_3D<T> &position, int per);
template<typename T> void get_per(Vec_3D<T> &position, int perx, int pery, int perz);
void get_per(std::vector<Particle_v>& particles, const unsigned par_num, const int per);

FTYPE get_distance(const Vec_3D<FTYPE> &x_1, const Vec_3D<FTYPE> &x_2, int per);
Vec_3D<FTYPE> get_sgn_distance(const Vec_3D<FTYPE> &x_from, const Vec_3D<FTYPE> &x_to, int per);

void assign_to(Mesh& field, const Vec_3D<FTYPE> &position, const FTYPE value);
void assign_to(std::vector<Mesh>& field, const Vec_3D<FTYPE> &position, const Vec_3D<FTYPE>& value);
void assign_from(const Mesh &field, const Vec_3D<FTYPE> &position, FTYPE& value);
void assign_from(const std::vector<Mesh> &field, const Vec_3D<FTYPE> &position, Vec_3D<FTYPE>& value);

void fftw_execute_dft_r2c(const FFTW_PLAN_TYPE &p_F, Mesh& rho);
void fftw_execute_dft_c2r(const FFTW_PLAN_TYPE &p_B, Mesh& rho);
void fftw_execute_dft_r2c_triple(const FFTW_PLAN_TYPE &p_F, std::vector<Mesh>& rho);
void fftw_execute_dft_c2r_triple(const FFTW_PLAN_TYPE &p_B, std::vector<Mesh>& rho);

template<typename T> T mean(const std::vector<T>& data)
{
    T tmp(0);
	
	#pragma omp parallel for reduction(+:tmp)
	for (auto it = data.begin(); it < data.end(); ++it) tmp += *it;
	
	return tmp / data.size();
}

inline FTYPE mean(const Mesh& data){ return mean(data.data); }

template<typename T> T std_dev(const std::vector<T>& data, T mean)
{
    T tmp(0);
	
	#pragma omp parallel for reduction(+:tmp)
	for (auto it = data.begin(); it < data.end(); ++it) tmp += pow2(*it-mean);
	
	return sqrt(tmp / data.size());
}

inline FTYPE std_dev(const Mesh& data, FTYPE mean){ return std_dev(data.data, mean); }

template<typename T> T min(const std::vector<T>& data)
{
    return *std::min_element(data.begin(), data.end());
}

inline FTYPE min(const Mesh& data){ return min(data.data); }

template<typename T> T max(const std::vector<T>& data)
{
    return *std::max_element(data.begin(), data.end());
}

inline FTYPE max(const Mesh& data){ return max(data.data); }

template<int points>
class IT
{
public:
    IT(const Vec_3D<FTYPE> &pos); // ctor for assignment scheme
    IT(const Vec_3D<FTYPE> &pos, FTYPE Hc); // ctor for chaining mesh

    unsigned counter;
    Vec_3D<int> vec;

    bool iter();
};

