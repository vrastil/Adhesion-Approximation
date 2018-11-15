# FastSim project
cmake_minimum_required(VERSION 3.2)
project (FastSim VERSION 0.1)

# add custom FindFFTW
list(APPEND CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/src/3party/findFFTW")

# check dependencies
find_package(GSL 2.1)
find_package(Boost 1.66)
find_package(FFTW 3.1)

# precision of the simulation
SET(PRECISION 2 CACHE STRING "Precision of the simulation")

# set compile flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} \
    -std=c++11 -pipe -MMD \
    -fopenmp -flto=jobserver -fPIC"
)

# set FFTW precision (single / double / long double)
if(${PRECISION} MATCHES 1)
    set(FFTW_LIB "-lfftw3f -lfftw3f_omp")
elseif(${PRECISION} MATCHES 2)
    set(FFTW_LIB "-lfftw3 -lfftw3_omp")
elseif(${PRECISION} MATCHES 3)
    set(FFTW_LIB "-lfftw3l -lfftw3l_omp")
else(${PRECISION} MATCHES 1)
    message(FATAL_ERROR "Invalid value of PRECISON (${PRECISION})")
endif(${PRECISION} MATCHES 1)

# build individual components
add_subdirectory (src/3party/)
add_subdirectory (src/ApproximationsSchemes)
add_subdirectory (src/core/)
add_subdirectory (src/data/)
add_subdirectory (src/io/)
add_subdirectory (src/main/)