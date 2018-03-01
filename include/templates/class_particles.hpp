/**
 * @file:	particles.hpp
 * @brief:	particle definitions  (with and without velocity)
 */

#pragma once
#include "stdafx.h"
#include "class_vec_3d.hpp"

/**
 * @class:	Particle_x
 * @brief:	class handling particles (position only)
 * @acces:	operator [] to get position coordinates
 */

template<typename T>
class Particle_x
{
public:
	// CONSTRUCTORS
	Particle_x(){};
    template<typename U>
	Particle_x(Vec_3D<U> position):
	position(position) {};
	
	// VARIABLES
	Vec_3D<T> position;
	
	// OPERATORS
	T &operator[](int i){ return position[i]; }
	const T& operator[](int i) const{ return position[i]; }
};

/**
 * @class:	Particle_v
 * @brief:	class handling particles (with velocitites)
 * @acces:	operator [] to get position coordinates
 * 			operator () to get velocity coordinates
 */

template<typename T>
class Particle_v : public Particle_x<T>
{
public:
	// CONSTRUCTORS
	Particle_v(){};
    template<typename U, typename V>
	Particle_v(Vec_3D<U> position, Vec_3D<V> velocity):
		Particle_x<T>(position), velocity(velocity) {};
	
	// VARIABLES
	Vec_3D<T> velocity;

	// OPERATORS
	T &operator()(int i){ return velocity[i]; }
	const T& operator()(int i) const{ return velocity[i]; }
};