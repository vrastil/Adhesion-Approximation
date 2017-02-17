
#include "stdafx.h"
#include "core.h"
#include <fftw3.h>

using namespace std;
const double PI = acos(-1.);

void get_k_vec(int N, int index, int* k_vec)
{
	k_vec[0] = index / ((N/2 + 1)*N);
	k_vec[1] = (index / (N/2 + 1)) % N;
	k_vec[2] = index % (N/2 + 1);
	
	for (int i =0; i<3; i++) k_vec[i] = ((k_vec[i]<=N/2.) ? k_vec[i] : k_vec[i] - N);
}

void get_k_vec(int N, int index, Vec_3D<int> &k_vec)
{
	k_vec[0] = index / ((N/2 + 1)*N);
	k_vec[1] = (index / (N/2 + 1)) % N;
	k_vec[2] = index % (N/2 + 1);
	
	for (int i =0; i<3; i++) k_vec[i] = ((k_vec[i]<=N/2.) ? k_vec[i] : k_vec[i] - N);
}

int get_k_sq(int N, int index)
{
	int k_vec[3];
	double tmp = 0;
	get_k_vec(N, index, k_vec);
	for (int i =0; i<3; i++) tmp += pow(k_vec[i],2);
	return tmp;
}

int get_per(int vec, int per)
{
	if (vec > per) vec = vec % per;
	else if (vec < 0) vec = vec % per + per;
	return vec;
}

void get_per(Vec_3D<double> &position, int per)
{
	for (int i = 0; i < 3; i++)
	{
		if (position[i] >= per) position[i] = fmod(position[i], per);
		else if (position[i] < 0) position[i] = fmod(position[i], per) + per;
	}
}

void get_per(Vec_3D<int> &position, int per)
{
	for (int i = 0; i < 3; i++)
	{
		if (position[i] >= per) position[i] = position[i] % per;
		else if (position[i] < 0) position[i] = position[i] % per + per;
	}
}

void get_per(Vec_3D<int> &position, const Vec_3D<int> &per)
{
	for (int i = 0; i < 3; i++)
	{
		if (position[i] >= per[i]) position[i] = position[i] % per[i];
		else if (position[i] < 0) position[i] = position[i] % per[i] + per[i];
	}
}

double wgh_sch(const Vec_3D<double> &x, Vec_3D<int> y, int mesh_num, const int order)
{
	// The weighting scheme used to assign values to the mesh points or vice versa
	// Return value of assigment function on mesh point y from particle in x
	double dx, w = 1;
	get_per(y, mesh_num);
	
	switch (order){
	case 0: {	// NGP: Nearest grid point
				for (int i = 0; i < 3; i++)
				{
					if ((int)x[i] != y[i]) w *= 0;
				}
				return w;
	}
	case 1: {	// CIC: Cloud in cells
				for (int i = 0; i < 3; i++)
				{
					dx = fmin(fmin(abs(x[i] - y[i]), x[i] + mesh_num - y[i]), y[i] + mesh_num - x[i]);
					if (dx > 1) w *= 0;
					else w *= 1 - dx;
				}
				return w;
	}
	case 2: {	// TSC: Triangular shaped clouds
				for (int i = 0; i < 3; i++)
				{
					dx = fmin(fmin(abs(x[i] - y[i]), x[i] + mesh_num - y[i]), y[i] + mesh_num - x[i]);
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
				#pragma omp atomic
				(*field)(y) += value * wgh_sch(position, y, field->N, order);
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
				#pragma omp atomic
				*value += field(y) * wgh_sch(position, y, field.N, order);
			}
		}
	}
}

void assign_from(const vector< Mesh> &field, const Vec_3D<double> &position, Vec_3D<double>* value, int order)
{
	for (int i = 0; i < 3; i++) assign_from(field[i], position, &((*value)[i]), order);
}

static inline void normalize_FFT_FORWARD(Mesh& rho)
{
	rho /= pow(rho.N, 1.5);
}

static inline void normalize_FFT_BACKWARD(Mesh& rho)
{
	rho /= pow(rho.N, 1.5);
}

void fftw_execute_dft_r2c(const fftw_plan &p_F, Mesh& rho)
{
	fftw_execute_dft_r2c(p_F, rho.real(), rho.complex());
	normalize_FFT_FORWARD(rho);
}

void fftw_execute_dft_c2r(const fftw_plan &p_B, Mesh& rho)
{
	fftw_execute_dft_c2r(p_B, rho.complex(), rho.real());
	normalize_FFT_BACKWARD(rho);
}

void fftw_execute_dft_r2c_triple(const fftw_plan &p_F, vector<Mesh>& rho)
{
	for (int i = 0; i < 3; i++) fftw_execute_dft_r2c(p_F, rho[i]);
}

void fftw_execute_dft_c2r_triple(const fftw_plan &p_B, vector<Mesh>& rho)
{
	for (int i = 0; i < 3; i++) fftw_execute_dft_c2r(p_B, rho[i]);
}

double mean(double* p_data, int len)
{
	double tmp = 0;
	
	#pragma omp parallel for reduction(+:tmp)
	for (int i = 0; i < len; i++) tmp += p_data[i];
	
	return tmp / len;
}

double std_dev(double* p_data, int len, double t_mean)
{
	double tmp = 0;
	
	#pragma omp parallel for reduction(+:tmp)
	for (int i = 0; i < len; i++) tmp += pow(p_data[i]-t_mean, 2.);
	
	return sqrt(tmp / len);
}

