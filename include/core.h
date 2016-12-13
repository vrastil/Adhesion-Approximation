#pragma once

#include "stdafx.h"
#include <fftw3.h>

/**
 * @class Mesh
 * @author Michal Vrastil
 * @date 12/12/16
 * @file core.h
 * @brief class handling basic mesh functions, the most important are creating and destroing the underlying data structure
 */

class Mesh
{
private:
	// VARIABLES
	double* data;
	
public:
	// CONSTRUCTORS & DESTRUCTOR
	Mesh(int n);
	~Mesh();
	
	// VARIABLES
	const int N, length; // acces dimensions and length of mesh
	
	// METHODS
	inline double* real() const { return data;} // acces data
	inline fftw_complex* complex() const { return reinterpret_cast<fftw_complex*>(data);}
	
	void get_k_vec(int index, int* k_vec) const;
	int get_k_sq(int index) const;
	
	// OPERATORS
	inline double &operator[](int i){ return data[i]; }
	inline const double &operator[](int i) const{ return data[i]; }
	
	inline double &operator()(int i, int j, int k){ return data[i*N*(N+2)+j*(N+2)+k]; }
	inline const double &operator()(int i, int j, int k) const{ return data[i*N*(N+2)+j*(N+2)+k]; }
	
	Mesh& operator+=(const double& rhs);
	Mesh& operator-=(const double& rhs){ return *this+=-rhs; };
	Mesh& operator*=(const double& rhs);
	Mesh& operator/=(const double& rhs);
};