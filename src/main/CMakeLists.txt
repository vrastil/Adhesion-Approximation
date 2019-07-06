# main interface library -- config.hpp, precison.hpp, stdafx.h
set(LIBRARY_NAME main)

# create interface
add_library(${LIBRARY_NAME} INTERFACE)

# publish headers for this library
target_include_directories(${LIBRARY_NAME} INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include")

# set compile flags
target_compile_definitions(${LIBRARY_NAME} INTERFACE PRECISION=${PRECISION})

# dependencies
target_link_libraries(${LIBRARY_NAME}
        INTERFACE boost_log
)

# config file
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/include/config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/include/config.h)

#############################
# main executable - FastSim #
#############################

set(EXECUTABLE FastSim)

# source files
set(SOURCE_FILES 
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
)

# create executable
add_executable(${EXECUTABLE} ${SOURCE_FILES})

# dependencies
target_link_libraries(${EXECUTABLE}
        Approximations
        ${FFTW_LIB}
        boost_log_setup
        boost_timer
        boost_thread
        -fuse-linker-plugin
)

# install Fast-Sim into binaries
install(TARGETS FastSim
        RUNTIME DESTINATION bin
)

##########################
# test executable - test #
##########################

set(EXECUTABLE test_all)

# source files
set(SOURCE_FILES 
    ${CMAKE_CURRENT_SOURCE_DIR}/test_main.cpp
)

# create executable
add_executable(${EXECUTABLE} EXCLUDE_FROM_ALL ${SOURCE_FILES})

# dependencies
target_link_libraries(${EXECUTABLE}
        main
        catch
        test_data
        test_core
)
