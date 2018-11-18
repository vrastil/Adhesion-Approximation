# core functions for simulations
set(LIBRARY_NAME core)

# source files
set(SOURCE_FILES 
    ${CMAKE_CURRENT_SOURCE_DIR}/core_app.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/core_mesh.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/core_power.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/integration.cpp
)

# create static library
add_library(${LIBRARY_NAME} STATIC ${SOURCE_FILES})

# set compile flags
target_compile_definitions(${LIBRARY_NAME} PRIVATE NOISE_HALF)

# dependencies
target_link_libraries(${LIBRARY_NAME}
        main
        Approximations
        CBRNG
        gsl
        ${FFTW_LIB}
)

# publish headers for this library
target_include_directories(${LIBRARY_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")