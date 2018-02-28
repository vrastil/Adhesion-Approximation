#pragma once
#include "stdafx.h"
#include <array>

/**
 * @class:	Data_Vec
 * @brief:	class containing data [x, y,...]
 */

template <typename T, unsigned N>
class Data_Vec
{
public:
    // CONSTRUCTORS
    Data_Vec() = default;
    Data_Vec(size_t size) { data.fill(std::vector<T>(size)); }

    // VARIABLES
    std::array<std::vector<T>, N> data;

    // ELEMENT ACCESS
    std::vector<T>& operator[](int i){ return data[i]; }
    const std::vector<T>& operator[](int i) const { return data[i]; }

    // CAPACITY
    size_t dim() const noexcept{ return data.size(); }
    size_t size() const noexcept{ return data[0].size(); }
    void resize (size_t n){
        for (auto &vec : data) vec.resize(n);
    }
    void resize (size_t n, T val){
        for (auto &vec : data) vec.resize(n, val);
    }
    void reserve(size_t n){
        for (auto &vec : data) vec.reserve(n);
    }
    void erase(unsigned index){
        for (auto &vec : data) vec.erase(vec.begin() + index);
    }
    // MODIFIERS
    void fill(T val){
        for (auto &vec : data) std::fill(vec.begin(), vec.end(), val);
    }
};