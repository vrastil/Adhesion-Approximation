/**
 * @brief define container Vec_3D
 * 
 * @file class_vec_3d.hpp
 * @author Michal Vrastil
 * @date 2018-07-11
 */

#pragma once
#include <array>

/**
 * @class:	Vec_3D<T>
 * @brief:	class handling basic 3D-vector functions, definitions
 */

template <typename T>
class Vec_3D : public std::array<T, 3>
{
public:
	// CONSTRUCTORS
	Vec_3D(){};

    template<typename ...U> Vec_3D(U...init):
    std::array<T, 3>{init...} {}

    template<typename U> Vec_3D(const Vec_3D<U>& vec):
    std::array<T, 3>({T(std::get<0>(vec)), T(std::get<1>(vec)), T(std::get<2>(vec))}) {}

    // METHODS
    T norm2() const
    {
        T tmp(0);
        for (T val : *this) tmp += val*val;
        return tmp;
    }

	auto norm() const -> decltype(sqrt(norm2())) { return sqrt(norm2()); }
		
	// OPERATORS
    template<typename U>
	Vec_3D<T>& operator+=(const Vec_3D<U>& rhs)
    {
        for(size_t i = 0; i < 3; ++i) (*this)[i] += rhs[i];
        return *this;
    }
    
    template<typename U>
    Vec_3D<T>& operator-=(const Vec_3D<U>& rhs)
    {
        for(size_t i = 0; i < 3; ++i) (*this)[i] -= rhs[i];
        return *this;
    }

    template<typename U>
	Vec_3D<T>& operator*=(U rhs)
    {
        for(T& val : *this) val *= rhs;
        return *this;
    }

    template<typename U>
	Vec_3D<T>& operator/=(U rhs)
    {
        for(T& val : *this) val /= rhs;
        return *this;
    }
};

// NON-MEMBER FUNCTIONS
template <typename T, typename U>
Vec_3D<T> operator+(Vec_3D<T> lhs, const Vec_3D<U>& rhs)
{
    lhs += rhs;
    return lhs;
}

template <typename T, typename U>
Vec_3D<T> operator-(Vec_3D<T> lhs, const Vec_3D<U>& rhs)
{
    lhs -= rhs;
    return lhs;
}

template <typename T, typename U>
Vec_3D<T> operator*(Vec_3D<T> lhs, U rhs)
{
    lhs *= rhs;
    return lhs;
}

template <typename T, typename U>
Vec_3D<T> operator*(T lhs, const Vec_3D<U>& rhs)
{
    Vec_3D<T> tmp;
    for(size_t i = 0; i < 3; ++i) tmp[i] = lhs*rhs[i];
    return tmp;
}

template <typename T, typename U>
Vec_3D<T> operator/(Vec_3D<T> lhs, U rhs)
{
    for(T& val : lhs) val /= rhs;
    return lhs;
}