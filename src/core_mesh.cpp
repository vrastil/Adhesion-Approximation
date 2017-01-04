
#include "stdafx.h"
#include "core.h"

using namespace std;
const double PI = acos(-1.);

void get_k_vec(int N, int index, int* k_vec)
{
	k_vec[0] = index / ((N/2 + 1)*N);
	k_vec[1] = (index / (N/2 + 1)) % N;
	k_vec[2] = index % (N/2 + 1);
}

int get_k_sq(int N, int index)
{
	int k_vec[3];
	double tmp = 0;
	get_k_vec(N, index, k_vec);
	for (int i =0; i<3; i++) tmp += ((k_vec[i]<N/2.) ? pow(k_vec[i],2) : pow(k_vec[i] - N, 2.));
	return tmp;
}

int get_per(int vec, int per)
{
	if (vec > per) vec = vec % per;
	else if (vec < 0) vec = vec % per + per;
	return vec;
}

double wgh_sch(const Vec_3D<double> &x, const Vec_3D<int> &y, int mesh_num, const int order)
{
	// The weighting scheme used to assign values to the mesh points or vice versa
	// Return value of assigment function on mesh point y from particle in x
	double dx, w = 1;
	int y_per;

	switch (order){
	case 0: {	// NGP: Nearest grid point
				for (int i = 0; i < 3; i++)
				{
					y_per = get_per(y[i], mesh_num);
					if ((int)x[i] != y_per) w *= 0;
				}
				return w;
	}
	case 1: {	// CIC: Cloud in cells
				for (int i = 0; i < 3; i++)
				{
					y_per = get_per(y[i], mesh_num);					
					dx = fmin(fmin(abs(x[i] - y_per), x[i] + mesh_num - y_per), y_per + mesh_num - x[i]);
					if (dx > 1) w *= 0;
					else w *= 1 - dx;
				}
				return w;
	}
	case 2: {	// TSC: Triangular shaped clouds
				for (int i = 0; i < 3; i++)
				{
					y_per = get_per(y[i], mesh_num);
					dx = fmin(fmin(abs(x[i] - y_per), x[i] + mesh_num - y_per), y_per + mesh_num - x[i]);
					if (dx > 1.5) w *= 0;
					else if (dx > 0.5) w *= (1.5 - dx)*(1.5 - dx) / 2.0;
					else w *= 3 / 4.0 - dx*dx;
				}
				return w;
	}
	}
	return 0;
}

void assign_to(Mesh* field, const Vec_3D<double> &position, const double value, int order)
{
	Vec_3D<int> y, z;
	for (int i = 0; i < 3; i++) z[i] = (int)(position[i] - 0.5*(order - 1));
	for (y[0] = z[0]; y[0] < z[0] + 1 + order; y[0]++)
	{
		for (y[1] = z[1]; y[1] < z[1] + 1 + order; y[1]++){
		
			for (y[2] = z[2]; y[2] < z[2] + 1 + order; y[2]++)
			{
				(*field)(y.x, y.y, y.z) += value * wgh_sch(position, y, field->N, order);
			}
		}
	}
}

void assign_from(const Mesh &field, const Vec_3D<double> &position, double* value, int order)
{
	Vec_3D<int> y, z;
	for (int i = 0; i < 3; i++) z[i] = (int)(position[i] - 0.5*(order - 1));
	for (y[0] = z[0]; y[0] < z[0] + 1 + order; y[0]++)
	{
		for (y[1] = z[1]; y[1] < z[1] + 1 + order; y[1]++){
		
			for (y[2] = z[2]; y[2] < z[2] + 1 + order; y[2]++)
			{
				*value += field(y.x, y.y, y.z) * wgh_sch(position, y, field.N, order);
			}
		}
	}
}