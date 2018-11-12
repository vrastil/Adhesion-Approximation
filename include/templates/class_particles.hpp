/**
 * @brief define container Particle (with and without velocity)
 * 
 * @file class_particles.hpp
 * @author Michal Vrastil
 * @date 2018-06-24
 */

#pragma once
#include "stdafx.h"
#include "class_vec_3d.hpp"

/**
 * @class:	Particle_x
 * @brief class handling particles (position only)
 * 
 * @tparam T coordinate type (float, double or long double)
 * 
 * operator [] to get position coordinates
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
	
    /**
     * @brief get position coordinates
     * 
     * @param i axis
     * @return T& position
     */
	T &operator[](unsigned int i){ return position[i]; }
	const T& operator[](int i) const{ return position[i]; }
};

/**
 * @class:	Particle_v
 * @brief class handling particles (position only)
 * 
 * @tparam T coordinate type (float, double or long double)
 * 
 * operator [] to get position coordinates
 * operator () to get velocity coordinates
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