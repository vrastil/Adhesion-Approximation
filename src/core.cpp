
#include "stdafx.h"
#include "core.h"

using namespace std;
const double PI = acos(-1.);

Mesh::Mesh(int n):N(n), length(n*n*(n+2))
{
	data = new double[length];
}

Mesh::~Mesh()
{
	delete[] data;
}

void Mesh::get_k_vec(int index, int* k_vec) const
{
	k_vec[0] = index / ((N/2 + 1)*N);
	k_vec[1] = (index / (N/2 + 1)) % N;
	k_vec[2] = index % (N/2 + 1);
}

int Mesh::get_k_sq(int index) const
{
	int k_vec[3];
	double tmp = 0;
	get_k_vec(index, k_vec);
	for (int i =0; i<3; i++) tmp += ((k_vec[i]<N/2.) ? pow(k_vec[i],2) : pow(k_vec[i] - N, 2.));
	return tmp;
}

Mesh& Mesh::operator+=(const double& rhs)
{
	for (int i = 0; i < length; i++) this->data[i]+=rhs;
	return *this;
}