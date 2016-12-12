#pragma once

#include "stdafx.h"
#include <fftw3.h>

/**
 * @class Mesh
 * @author Michal Vrastil
 * @date 12/12/16
 * @file core.h
 * @brief class handling basic mesh functions, the most important are creating and destroing underlying data structure
 */

class Mesh
{
public:
	// CONSTRUCTORS
	c_Mesh(int n);
	
	// VARIABLES
	const int mesh_num, length; // acces dimensions and length of mesh
	double* data;
	
	// METHODS
	double* data(){return p_data.data();}; // acces data
	double get_val(int x, int y, int z) const {return p_data[index(x, y, z)];}
	int index(int x, int y, int z) const { return x*N2*N3+y*N3+z; }
	int index(int y[3]) const { return index(y[0], y[1], y[2]); }
	int index(int xy, int z) const { return xy*N3+z; }
	// get coordinates
	int get_x(int index) const { return index / (N2*N3); }
	int get_y(int index) const { return (index / N3) % N2; }
	int get_z(int index) const { return index % N3; }
	// get coordinates assuming complex numbers (half)
	int get_x_c(int index) const { return index / (N2*(N3/2)); }
	int get_y_c(int index) const { return (index / (N3/2)) % N2; }
	int get_z_c(int index) const{ return index % (N3/2); }
	int get_kx(int index) const;
	int get_ky(int index) const;
	int get_kz(int index) const;
	void get_k_vec(int index, int* k_vec) const;
	int get_k_sq(int index) const;
	void assign_to(const std::vector<double> &position, const double value, int order);
	void assign_from(const std::vector<double> &position, double* value, int order) const;
	
	// OPERATORS
	double &operator[](int i){ return p_data[i]; }
	const double &operator[](int i) const{ return p_data[i]; }
	Mesh& operator+=(const double& rhs);
	Mesh& operator-=(const double& rhs){ return *this+=-rhs; };
	Mesh& operator*=(const double& rhs);
	Mesh& operator/=(const double& rhs);
};