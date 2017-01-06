/**
 * @file:	core_fce.h
 * @brief:	basic functions to work with mesh
 */
 
#include "stdafx.h"
#include "core.h"

void get_k_vec(int N, int index, int* k_vec);
int get_k_sq(int N, int index);

void assign_to(Mesh* field, const Vec_3D<double> &position, const double value, int order);
void assign_from(const Mesh &field, const Vec_3D<double> &position, double* value, int order);