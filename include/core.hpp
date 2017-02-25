/**
 * @file:	core.hpp
 * @brief:	template class functions definitions
 */
 
#include "core_mesh.h"

/**
 * @class:	Vec_3D<T>
 * @brief:	class handling basic 3D-vector functions
 */

template <typename T>
void Vec_3D<T>::assign(T x_, T y_, T z_)
{
	x = x_;
	y = y_;
	z = z_;
}

template <typename T>
T& Vec_3D<T>::operator[](int i)
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

template <typename T>
const T& Vec_3D<T>::operator[](int i) const
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

template <typename T>
Vec_3D<T>& Vec_3D<T>::operator+=(const Vec_3D<T>& rhs)
{
	x+=rhs.x;
	y+=rhs.y;
	z+=rhs.z;
	return *this;
}

template <typename T>
Vec_3D<T> operator+(Vec_3D<T> lhs, const Vec_3D<T>& rhs)
{
	lhs += rhs;
	return lhs;
}

template <typename T>
Vec_3D<T>& Vec_3D<T>::operator-=(const Vec_3D<T>& rhs)
{
	x-=rhs.x;
	y-=rhs.y;
	z-=rhs.z;
	return *this;
}

template <typename T>
Vec_3D<T> operator-(Vec_3D<T> lhs, const Vec_3D<T>& rhs)
{
	lhs -= rhs;
	return lhs;
}

template <typename T>
Vec_3D<T>& Vec_3D<T>::operator*=(T rhs)
{
	x*=rhs;
	y*=rhs;
	z*=rhs;
	return *this;
}

template <typename T>
Vec_3D<T> operator*(Vec_3D<T> lhs, T rhs)
{
	lhs *= rhs;
	return lhs;
}

template <typename T>
Vec_3D<T> operator*(T lhs, Vec_3D<T> rhs)
{
	rhs *= lhs;
	return rhs;
}

template <typename T>
Vec_3D<T>& Vec_3D<T>::operator/=(T rhs)
{
	x/=rhs;
	y/=rhs;
	z/=rhs;
	return *this;
}

template <typename T>
Vec_3D<T> operator/(Vec_3D<T> lhs, T rhs)
{
	lhs /= rhs;
	return lhs;
}

template <typename T>
template<class U>
Vec_3D<T>::operator Vec_3D<U>() const
{
	Vec_3D<U> lhs;
	lhs.x = static_cast<U>(this->x);
	lhs.y = static_cast<U>(this->y);
	lhs.z = static_cast<U>(this->z);
	return lhs;
}

/**
 * @class:	Mesh_base<T>
 * @brief:	class handling basic mesh functions, the most important are creating and destroing the underlying data structure
 *			creates a mesh of N1*N2*N3 cells
 */

template <typename T>
Mesh_base<T>::Mesh_base(int n1, int n2, int n3):N1(n1), N2(n2), N3(n3), length(n1*n2*n3)
{
	data = new T[length];
//	printf("Normal base ctor %p, N1 = %i, N2 = %i, N3 = %i\n", this, N1, N2, N3); 
}

template <typename T>
Mesh_base<T>::Mesh_base(const Mesh_base<T>& that): N1(that.N1), N2(that.N2), N3(that.N3), length(that.length)
{
	data = new T[length];
	
	#pragma omp parallel for
	for (int i = 0; i < length; i++) data[i] = that.data[i];
//	printf("Copy base ctor %p\n", this);
}

template <typename T>
void swap(Mesh_base<T>& first, Mesh_base<T>& second)
{
	std::swap(first.length, second.length);
	std::swap(first.N1, second.N1);
	std::swap(first.N2, second.N2);
	std::swap(first.N3, second.N3);
	std::swap(first.data, second.data);
}

template <typename T>
Mesh_base<T>& Mesh_base<T>::operator=(const Mesh_base<T>& other)
{
//	printf("Copy base assignemnt %p\n", this);
	Mesh_base<T> temp(other);
	swap(*this, temp);
    return *this;
}

template <typename T>
Mesh_base<T>::~Mesh_base<T>()
{
	delete[] data;
//	printf("dtor base %p\n", this);
}

template <typename T>
T& Mesh_base<T>::operator()(Vec_3D<int> pos)
{
	get_per(pos, N1, N2, N3);
	return data[pos.x*N2*N3+pos.y*N3+pos.z]; 
}

template <typename T>
const T& Mesh_base<T>::operator()(Vec_3D<int> pos) const
{
	get_per(pos, N1, N2, N3);
	return data[pos.x*N2*N3+pos.y*N3+pos.z];
}

template <typename T>
Mesh_base<T>& Mesh_base<T>::operator+=(const T& rhs)
{
	#pragma omp parallel for
		for (int i = 0; i < length; i++) this->data[i]+=rhs;
		
	return *this;
}

template <typename T>
Mesh_base<T>& Mesh_base<T>::operator*=(const T& rhs)
{
	#pragma omp parallel for
		for (int i = 0; i < length; i++) this->data[i]*=rhs;
		
	return *this;
}

template <typename T>
Mesh_base<T>& Mesh_base<T>::operator/=(const T& rhs)
{
	#pragma omp parallel for
		for (int i = 0; i < length; i++) this->data[i]/=rhs;
		
	return *this;
}

template <typename T>
void Mesh_base<T>::assign(T val)
{
	#pragma omp parallel for
		for (int i = 0; i < length; i++) this->data[i]=val;
}