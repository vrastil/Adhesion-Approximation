/**
 * @brief basic functions to work with mesh
 * 
 * @file core_mesh.cpp
 * @author Michal Vrastil
 * @date 2018-07-11
 */

#include "core_mesh.h"

#ifndef ORDER
#define ORDER 1
#endif

template <typename T> static int sgn(T val)
{
	return (T(0) < val) - (val < T(0));
}

void get_k_vec(size_t N, size_t index, int* k_vec)
{
	k_vec[0] = index / ((N/2 + 1)*N);
	k_vec[1] = (index / (N/2 + 1)) % N;
	k_vec[2] = index % (N/2 + 1);
	
	for (unsigned int i = 0; i<2; i++) if (k_vec[i] > N/2) k_vec[i] -= N; // k_vec[2] is ALWAYS less or equal than N/2 (real FFTW)
}

void get_k_vec(size_t N, size_t index, Vec_3D<int> &k_vec)
{
	k_vec[0] = index / ((N/2 + 1)*N);
	k_vec[1] = (index / (N/2 + 1)) % N;
	k_vec[2] = index % (N/2 + 1);
	
	for (unsigned int i = 0; i<2; i++) if (k_vec[i] > N/2) k_vec[i] -= N; // k_vec[2] is ALWAYS less or equal than N/2 (real FFTW)
}

FTYPE_t get_k_sq(size_t N, size_t index)
{
	int k_vec[3];
	FTYPE_t tmp = 0;
	get_k_vec(N, index, k_vec);
	for (unsigned int i = 0; i<3; i++) tmp += pow2(k_vec[i]);
	return tmp;
}

template<typename T>
static typename std::enable_if<std::is_integral<T>::value, T>::type get_per(T vec, size_t per)
{
    static_assert(std::is_integral<T>::value, "Integral required.");
    if ((vec >= per) || (vec < 0) ){
        vec %= per;
        return (vec < 0) ? vec + per : vec;
    }
    else return vec;
}

template<typename T>
static typename std::enable_if<!std::is_integral<T>::value, T>::type get_per(T vec, size_t per)
{
    static_assert(!std::is_integral<T>::value, "Non-integral required.");
    return ((vec >= per) || (vec < 0) ) ? vec - per * floor( vec / per ) : vec;
}

template<typename T>
void get_per(Vec_3D<T> &position, size_t per)
{
    for (T& pos : position) pos = get_per(pos, per);
}

template<typename T>
void get_per(Vec_3D<T> &position, size_t perx, size_t pery, size_t perz)
{
    position[0] = get_per(position[0], perx);
    position[1] = get_per(position[1], pery);
    position[2] = get_per(position[2], perz);
}

void get_per(std::vector<Particle_v<FTYPE_t>>& particles, const size_t per)
{
    const size_t Np = particles.size();
    #pragma omp parallel for
    for (size_t i = 0; i < Np; i++)
    {
        get_per(particles[i].position, per);
    }
}

template<typename T>
inline T get_distance_1D(const T x_1, const T x_2, const size_t per)
{	
	T d = std::abs(x_2 - x_1);
	if (d <= per / 2) return d; // most probable, first condition
	else if (d <= per) return per - d;
	else d = fmod(d, T(per)); // fmod unlikely to evaluate, speed up code
	if (d <= per / 2) return d;
	else return per - d;
}

template<typename T>
inline T get_distance_1D(const T x_1, const int x_2, const size_t per){ return get_distance_1D(x_1, (T)x_2, per); }

template<typename T>
inline T get_sgn_distance_1D(T x_from, T x_to, size_t per)
{	// return signed (oriented) distance, e.g. x_to - x_from (no periodicity)
	// periodicity makes this a little bit trickier
	T d = x_from - x_to;
	if (std::abs(d) <= per / 2) return d; // most probable, first condition
	else if (std::abs(d) <= per) return d-sgn(d)*per;
	else d = fmod(d, T(per)); // fmod unlikely to evaluate, speed up code
	if (std::abs(d) <= per / 2) return d;
	else return d-sgn(d)*per;
}

FTYPE_t get_distance(const Vec_3D<FTYPE_t> &x_1, const Vec_3D<FTYPE_t> &x_2, int per)
{
	FTYPE_t dst = 0;
	for (unsigned int i = 0; i < 3; i++) dst+= pow2(get_distance_1D(x_1[i], x_2[i], per));
	return sqrt(dst);
}

Vec_3D<FTYPE_t> get_sgn_distance(const Vec_3D<FTYPE_t> &x_from, const Vec_3D<FTYPE_t> &x_to, size_t per)
{
	Vec_3D<FTYPE_t> dst;
	for (unsigned int i = 0; i < 3; i++) dst[i] = get_sgn_distance_1D(x_from[i], x_to[i], per);
	return dst;
}

template<size_t N> static FTYPE_t wgh_sch(const Vec_3D<FTYPE_t> &x, const Vec_3D<int>& y, size_t mesh_num);
// The weighting scheme used to assign values to the mesh points or vice versa
// Return value of assigment function on mesh point y from particle in x

// template<> FTYPE_t wgh_sch<0>(const Vec_3D<FTYPE_t> &x, Vec_3D<int> y, size_t mesh_num)
// { // NGP: Nearest grid point
//     get_per(y, mesh_num);
//     for (unsigned int i = 0; i < 3; i++) if ((int)x[i] != y[i]) return 0;
//     return 1;
// }

template<> FTYPE_t wgh_sch<1>(const Vec_3D<FTYPE_t> &x, const Vec_3D<int>& y, size_t mesh_num)
{ // CIC: Cloud in cells
    FTYPE_t dx, w = 1;
    for (unsigned int i = 0; i < 3; i++)
    {
        dx = get_distance_1D(x[i], y[i], mesh_num);
        if (dx > 1) return 0;
        else w *= 1 - dx;
    }
    return w;
}

template<> FTYPE_t wgh_sch<2>(const Vec_3D<FTYPE_t> &x, const Vec_3D<int>& y, size_t mesh_num)
{ // TSC: Triangular shaped clouds
    FTYPE_t dx, w = 1;
    for (unsigned int i = 0; i < 3; i++)
    {
        dx = get_distance_1D(x[i], y[i], mesh_num);
        if (dx > 1.5) return 0;
        else if (dx > 0.5) w *= pow2(1.5 - dx) / 2;
        else w *= 3 / FTYPE_t(4) - dx*dx;
    }
    return w;
}

/**
 * @class:	IT
 * @brief:	class for effective iteration of cube of mesh cells
 */

template<unsigned int points>
IT<points>::IT(const Vec_3D<FTYPE_t> &pos): counter(0)
{
    for(size_t i = 0; i < 3; i++){
        vec[i] = (int)(pos[i] - 0.5*(int(points) - 2));
    }
}

template<unsigned int points>
IT<points>::IT(const Vec_3D<FTYPE_t> &pos, FTYPE_t Hc): counter(0)
{
    for(size_t i = 0; i < 3; i++){
        vec[i] = (int)(pos[i]/Hc) - 1;
    }
}

template<unsigned int points>
bool IT<points>::iter(){
    if (++counter == (points*points*points)) return false;
    vec[2]++;
    if ((counter % points) == 0){
        vec[2] -= points;
        vec[1]++;
        if ((counter % (points*points)) == 0){
            vec[1] -= points;
            vec[0]++;
        }
    }
    return true;
}

void assign_to(Mesh& field, const Vec_3D<FTYPE_t> &position, const FTYPE_t value)
{
    IT<ORDER+1> it(position);
    do{
        #pragma omp atomic
        field(it.vec) += value * wgh_sch<ORDER>(position, it.vec, field.N);
    } while( it.iter() );
}

void assign_to(std::vector<Mesh>& field, const Vec_3D<FTYPE_t> &position, const Vec_3D<FTYPE_t>& value)
{
    IT<ORDER+1> it(position);
    FTYPE_t w;
    do{
        w = wgh_sch<ORDER>(position, it.vec, field[0].N); ///< reuse the same weigh for every field in std::vector
        for (size_t i = 0; i < 3; i++)
        {
            #pragma omp atomic
            field[i](it.vec) += value[i] * w;
        }
	} while( it.iter() );
}

void assign_from(const Mesh &field, const Vec_3D<FTYPE_t> &position, FTYPE_t& value, FTYPE_t mod)
{
	IT<ORDER+1> it(position);
    do{
        #pragma omp atomic
        value += field(it.vec) * mod * wgh_sch<ORDER>(position, it.vec, field.N);
	} while( it.iter() );
}

void assign_from(const std::vector<Mesh> &field, const Vec_3D<FTYPE_t> &position, Vec_3D<FTYPE_t>& value, FTYPE_t mod)
{
    IT<ORDER+1> it(position);
    FTYPE_t w;
    do{
        w = mod * wgh_sch<ORDER>(position, it.vec, field[0].N); ///< reuse the same weigh for every field in std::vector
        for (unsigned int i = 0; i < 3; i++)
        {
            #pragma omp atomic
            value[i] += field[i](it.vec) * w;
        }
	} while( it.iter() );
}

void fftw_execute_dft_r2c(const FFTW_PLAN_TYPE &p_F, Mesh& rho)
{
	FFTW_EXEC_R2C(p_F, rho.real(), rho.complex());
	rho /= pow((FTYPE_t)rho.N, 3); //< normalization
}

void fftw_execute_dft_c2r(const FFTW_PLAN_TYPE &p_B, Mesh& rho)
{
	FFTW_EXEC_C2R(p_B, rho.complex(), rho.real());
}

void fftw_execute_dft_r2c_triple(const FFTW_PLAN_TYPE &p_F, std::vector<Mesh>& rho)
{
	for (unsigned int i = 0; i < 3; i++) fftw_execute_dft_r2c(p_F, rho[i]);
}

void fftw_execute_dft_c2r_triple(const FFTW_PLAN_TYPE &p_B, std::vector<Mesh>& rho)
{
	for (unsigned int i = 0; i < 3; i++) fftw_execute_dft_c2r(p_B, rho[i]);
}

template void get_per(Vec_3D<int>&, size_t);
template void get_per(Vec_3D<int>&, size_t, size_t, size_t);
template void get_per(Vec_3D<size_t>&, size_t);
template void get_per(Vec_3D<size_t>&, size_t, size_t, size_t);
template void get_per(Vec_3D<FTYPE_t>&, size_t);
template void get_per(Vec_3D<FTYPE_t>&, size_t, size_t, size_t);

template class IT<ORDER+1>;
template class IT<3>;