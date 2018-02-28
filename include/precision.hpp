/**
 * @file:	precision.hpp
 * @brief:	single / double / long double definitions
 */

#pragma once

#ifndef PRECISION
#define PRECISION 2 // default double precision
#endif

#if PRECISION == 1
typedef float FTYPE;
#define MAKE_FFTW_NAME(FUNC_NAME) fftwf_ ## FUNC_NAME
#elif PRECISION == 2
typedef double FTYPE;
#define MAKE_FFTW_NAME(FUNC_NAME) fftw_ ## FUNC_NAME
#elif PRECISION == 3
typedef long double FTYPE;
#define MAKE_FFTW_NAME(FUNC_NAME) fftwl_ ## FUNC_NAME
#endif

#define FFTW_PLAN_TYPE MAKE_FFTW_NAME(plan)
#define FFTW_DEST_PLAN MAKE_FFTW_NAME(destroy_plan)
#define FFTW_COMPLEX_TYPE MAKE_FFTW_NAME(complex)
#define FFTW_PLAN_R2C MAKE_FFTW_NAME(plan_dft_r2c_3d)
#define FFTW_PLAN_C2R MAKE_FFTW_NAME(plan_dft_c2r_3d)
#define FFTW_PLAN_OMP MAKE_FFTW_NAME(plan_with_nthreads)
#define FFTW_PLAN_OMP_INIT MAKE_FFTW_NAME(init_threads)
#define FFTW_PLAN_OMP_CLEAN MAKE_FFTW_NAME(cleanup_threads)
#define FFTW_EXEC_R2C MAKE_FFTW_NAME(execute_dft_r2c)
#define FFTW_EXEC_C2R MAKE_FFTW_NAME(execute_dft_c2r)

constexpr FTYPE PI = FTYPE(3.14159265358979323846); // 20 digits

inline float pow(float base, unsigned exp)
{
    float result = 1.f;
    while (exp)
    {
        if (exp & 1)
            result *= base;
        exp >>= 1;
        base *= base;
    }
    return result;
}

template <typename T> inline T pow2(T base){ return base*base; } //< most often used