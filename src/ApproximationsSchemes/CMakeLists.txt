# various approximation schemes
set(LIBRARY_NAME Approximations)

# source files
set(SOURCE_FILES 
    ${CMAKE_CURRENT_SOURCE_DIR}/app_var.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/adhesion.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/chameleon.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/frozen_flow.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/frozen_potential.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/mod_frozen_potential.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/zeldovich.cpp
)

# create static library
add_library(${LIBRARY_NAME} STATIC ${SOURCE_FILES})

# dependencies
target_link_libraries(${LIBRARY_NAME}
        main
        data
        ccl
        MultiGridSolver
        ${FFTW_LIB}
)

# publish headers for this library
target_include_directories(${LIBRARY_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")