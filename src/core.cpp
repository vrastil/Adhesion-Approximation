
#include "stdafx.h"
#include "core.h"
#include "core_cmd.h"

using namespace std;
const double PI = acos(-1.);

/**
 * @class:	Mesh
 * @brief:	class handling basic mesh functions, the most important are creating and destroing the underlying data structure
 *			creates a mesh of N*N*(N+2) cells
 */

Mesh::Mesh(int n):N(n), length(n*n*(n+2))
{
	data = new double[length];
}

Mesh::~Mesh()
{
	delete[] data;
}

Mesh& Mesh::operator+=(const double& rhs)
{
	#pragma omp parallel for
		for (int i = 0; i < length; i++) this->data[i]+=rhs;
		
	return *this;
}

Mesh& Mesh::operator*=(const double& rhs)
{
	#pragma omp parallel for
		for (int i = 0; i < length; i++) this->data[i]*=rhs;
		
	return *this;
}

Mesh& Mesh::operator/=(const double& rhs)
{
	#pragma omp parallel for
		for (int i = 0; i < length; i++) this->data[i]/=rhs;
		
	return *this;
}

/**
 * @class:	Vec_3D
 * @brief:	class handling basic 3D-vector functions
 */
 
double& Vec_3D::operator[](int i)
{
	switch(i)
	{
		case 0 : return x;
		case 1 : return y;
		case 2 : return z;
		default:
		{
			printf("Invalid acces in class Vec_3D. Invalid postion '%d'.\n", i);
			if (i < 0) return x;
			else return z;
		}
	}
}

const double& Vec_3D::operator[](int i) const
{
	switch(i)
	{
		case 0 : return x;
		case 1 : return y;
		case 2 : return z;
		default:
		{
			printf("Invalid acces in class Vec_3D. Invalid postion '%d'.\n", i);
			if (i < 0) return x;
			else return z;
		}
	}
}

/**
 * @class:	Sim_Param
 * @brief:	class storing simulation parameters
 */
 
 int Sim_Param::init(int ac, char* av[])
{
	int err = handle_cmd_line(ac, av, this);
	if (err) {is_init = 0; return err;}
	else {
		is_init = 1;
		par_num = pow(mesh_num / Ng, 3);
		power.k2_G *= power.k2_G;
		power.eval_pwr();
		b_in = 1./(z_in + 1);
		b_out = 1./(z_out + 1);
		k_min = 2.*PI/box_size;
		k_max = 2.*PI*mesh_num/box_size;
		return err;
	}
}

void Sim_Param::print_info()
{
	if (is_init) 
	{
		printf("\n");
		printf("Num_par:\t%i\n", Ng);
		printf("Num_mesh:\t%i^3\n", mesh_num);
		printf("Box size:\t%i Mpc/h\n", box_size);
		printf("Starting redshift:\t%G\n", z_in);
		printf("The primordial power spectrum 'P(k)=A*k^ns' has amplitude A = %G and spectral index ns = %G.\n", power.A, power.ns);
		if (power.k2_G == 0) printf("Smoothing length was not set.\n");
		else printf("Smoothing wavenumber is %G h/Mpc.\n", sqrt(power.k2_G));
		printf("'viscozity' for adhesion approximation is %G px^2.\n", nu);
		printf("The program will try to use %i threads.\n", nt);
		cout << "Output will be written to folder '"<< out_dir << "'\n";
		printf("\n");
	}
	else printf("WARNING! Simulation parameters are not initialized!\n");
}