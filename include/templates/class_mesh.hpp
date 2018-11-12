/**
 * @brief define container Mesh
 * 
 * @file class_mesh.hpp
 * @author Michal Vrastil
 * @date 2018-06-24
 */

#pragma once
#include "stdafx.h"
#include "../precision.hpp"
#include "class_vec_3d.hpp"

// from "core_mesh.hpp"
template<typename T> void get_per(Vec_3D<T> &position, size_t per);
template<typename T> void get_per(Vec_3D<T> &position, size_t perx, size_t pery, size_t perz);

/**
 * @class:	Mesh_base
 * @brief:	class handling basic mesh functions, the most important are creating and destroing the underlying data structure
 *			creates a mesh of N1*N2*N3 cells
 */
template <typename T>
class Mesh_base
{
public:
	// CONSTRUCTOR
	Mesh_base(size_t n1, size_t n2, size_t n3):
    N1(n1), N2(n2), N3(n3), length(n1*n2*n3), data(length) {}
	
	// VARIABLES
	size_t N1, N2, N3, length; // acces dimensions and length of mesh
    std::vector<T> data; // data stored on the mesh
	
	// METHODS
    T* real() { return data.data();} // acces data through pointer
    const T* real() const { return data.data();} // acces data through const pointer
	void assign(T val)
    {
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++) this->data[i]=val;
    }
	
	// OPERATORS
	T &operator[](size_t i){ return data[i]; }
	const T &operator[](size_t i) const{ return data[i]; }
	
	T& operator()(size_t i, size_t j, size_t k){ return data[i*N2*N3+j*N3+k]; }
	const T& operator()(size_t i, size_t j, size_t k) const{ return data[i*N2*N3+j*N3+k]; }
	
	T& operator()(size_t i, size_t j){ return data[i*N3+j]; }
	const T& operator()(size_t i, size_t j) const{ return data[i*N3+j]; }

    template<typename U> T& operator()(Vec_3D<U> pos)
    {
        get_per(pos, N1, N2, N3);
        return data[size_t(pos[0])*N2*N3+size_t(pos[1])*N3+size_t(pos[2])]; 
    }

	template<typename U> const T& operator()(Vec_3D<U> pos) const
    {
        get_per(pos, N1, N2, N3);
        return data[size_t(pos[0])*N2*N3+size_t(pos[1])*N3+size_t(pos[2])]; 
    }
	
	Mesh_base& operator+=(const T& rhs)
    {
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++) this->data[i]+=rhs; 
        return *this;
    }

	Mesh_base& operator-=(const T& rhs){ return *this+=-rhs; }

	Mesh_base& operator*=(const T& rhs)
    {
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++) this->data[i]*=rhs; 
        return *this;
    }

	Mesh_base& operator/=(const T& rhs)
    {
        #pragma omp parallel for
        for (size_t i = 0; i < length; i++) this->data[i]/=rhs; 
        return *this;
    }

};

/**
 * @class:	Mesh
 * @brief:	creates a mesh of N*N*(N+2) cells
 */
class Mesh : public Mesh_base<FTYPE_t>
{
public:
	// CONSTRUCTORS & DESTRUCTOR
    Mesh(size_t n): Mesh_base(n, n, n+2), N(n) {}
	
	// VARIABLES
	size_t N; // acces dimension of mesh
	
	// METHODS

    /**
     * @brief get fftw_complex pointer to data
     * 
     * @return FFTW_COMPLEX_TYPE* 
     */
    FFTW_COMPLEX_TYPE* complex() { return reinterpret_cast<FFTW_COMPLEX_TYPE*>(data.data());}

    /**
     * @brief get const fftw_complex pointer to data
     * 
     * @return FFTW_COMPLEX_TYPE* 
     */
    const FFTW_COMPLEX_TYPE* complex() const { return reinterpret_cast<const FFTW_COMPLEX_TYPE*>(data.data());}

    void reset_part(bool part)
    {/* nullify real (part = 0) or complex (part = 1) part of a field */
        #pragma omp parallel for
        for (size_t i = part; i < this->length; i+=2){
            data[i] = 0;
        }
    }

    void reset_re() { reset_part(0); }
    void reset_im() { reset_part(1); }
    
	// OPERATORS
	using Mesh_base<FTYPE_t>::operator ();

    template<typename U> FTYPE_t& operator()(Vec_3D<U> pos)
    {
        get_per(pos, N);
        return data[size_t(pos[0])*N2*N3+size_t(pos[1])*N3+size_t(pos[2])]; 
    }

	template<typename U> const FTYPE_t& operator()(Vec_3D<U> pos) const
    {
        get_per(pos, N);
        return data[size_t(pos[0])*N2*N3+size_t(pos[1])*N3+size_t(pos[2])]; 
    }
};