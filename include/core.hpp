/**
 * @file:	core.hpp
 * @brief:	template class functions definitions
 */
 
#include "core_mesh.h"
#include "core_out.h"
#include "core_app.h"

/**
 * @class:	Vec_3D<T>
 * @brief:	class handling basic 3D-vector functions
 */

template <typename T>
double Vec_3D<T>::norm() const
{
    T tmp(0);
    for (const T val : vec)
    {
        tmp += val*val;
    }
    return sqrt(tmp);
}

template <typename T>
Vec_3D<T>& Vec_3D<T>::operator+=(const Vec_3D<T>& rhs)
{
    for(unsigned i = 0; i < 3; ++i)
    {
        vec[i] += rhs[i];
    }
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
    for(unsigned i = 0; i < 3; ++i)
    {
        vec[i] -= rhs[i];
    }
	return *this;
}

template <typename T>
Vec_3D<T> operator-(Vec_3D<T> lhs, const Vec_3D<T>& rhs)
{
	lhs -= rhs;
	return lhs;
}

template <typename T>
Vec_3D<T>& Vec_3D<T>::operator+=(T rhs)
{
    for(T& val : vec)
    {
        val += rhs;
    }
	return *this;
}

template <typename T>
Vec_3D<T>& Vec_3D<T>::operator-=(T rhs)
{
    for(T& val : vec)
    {
        val -= rhs;
    }
	return *this;
}

template <typename T>
Vec_3D<T>& Vec_3D<T>::operator*=(T rhs)
{
    for(T& val : vec)
    {
        val *= rhs;
    }
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
Vec_3D<T> operator+(Vec_3D<T> lhs, T rhs)
{
	lhs += rhs;
	return lhs;
}

template <typename T>
Vec_3D<T> operator+(T lhs, Vec_3D<T> rhs)
{
	rhs += lhs;
	return rhs;
}

template <typename T>
Vec_3D<T> operator-(Vec_3D<T> lhs, T rhs)
{
	lhs -= rhs;
	return lhs;
}

template <typename T>
Vec_3D<T> operator-(T lhs, Vec_3D<T> rhs)
{
	rhs -= lhs;
	return rhs;
}

template <typename T>
Vec_3D<T>& Vec_3D<T>::operator/=(T rhs)
{
    for(T& val : vec)
    {
        val /= rhs;
    }
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
    for(unsigned i = 0; i < 3; ++i)
    {
        lhs[i] = static_cast<U>((*this)[i]);
    }
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
	return data[pos[0]*N2*N3+pos[1]*N3+pos[2]]; 
}

template <typename T>
const T& Mesh_base<T>::operator()(Vec_3D<int> pos) const
{
	get_per(pos, N1, N2, N3);
	return data[pos[0]*N2*N3+pos[1]*N3+pos[2]];
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


/**
 * @class:	Tracking
 * @brief:	class storing info about tracked particles
 */

template <class T>
void Tracking::update_track_par(T* particles)
{
	std::vector<Particle_x> par_pos_step;
	par_pos_step.reserve(num_track_par);
	for (int i=0; i<num_track_par; i++){
		par_pos_step.push_back(particles[par_ids[i]]);
	}
	par_pos.push_back(par_pos_step);
}

/**
 * @class:	App_Var_base
 * @brief:	class containing variables for approximations
 */

template <class T> 
void App_Var_base::print(const Sim_Param &sim, std::string out_dir_app, T* particles)
{
	/* Printing positions */
	print_par_pos_cut_small(particles, sim, out_dir_app, z_suffix());
	print_track_par(track, sim, out_dir_app, z_suffix());

	/* Printing density */
    get_rho_from_par(particles, &power_aux, sim);
    gen_dens_binned(power_aux, dens_binned, sim);    
	print_rho_map(power_aux, sim, out_dir_app, z_suffix());
	print_dens_bin(dens_binned, sim.mesh_num, out_dir_app, z_suffix());

	/* Printing power spectrum */
	fftw_execute_dft_r2c(p_F_pwr, power_aux);
	pwr_spec_k(sim, power_aux, &power_aux);
	gen_pow_spec_binned(sim, power_aux, &pwr_spec_binned);
    print_pow_spec(pwr_spec_binned, out_dir_app, z_suffix());
    if (!is_init_pwr_spec_0){
        pwr_spec_binned_0 = pwr_spec_binned;
        b_init = b;
        is_init_pwr_spec_0 = true;
    }
	print_pow_spec_diff(pwr_spec_binned, pwr_spec_binned_0, b / b_init, out_dir_app, z_suffix());
}