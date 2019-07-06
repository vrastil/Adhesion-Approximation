# library for various data types
set(LIBRARY_NAME data)

# source files
set(SOURCE_FILES 
    ${CMAKE_CURRENT_SOURCE_DIR}/params.cpp
)

# create static library
add_library(${LIBRARY_NAME} STATIC ${SOURCE_FILES})

# dependencies
target_link_libraries(${LIBRARY_NAME}
        main
        ccl
        io
        json
        ${FFTW_LIB}
)

# publish headers for this library
target_include_directories(${LIBRARY_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")

# unit tests
set(TEST_NAME test_data)

# source files
set(SOURCE_FILES 
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_vec_3d.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/tests/test_mesh.cpp
)

# create static library
add_library(${TEST_NAME} EXCLUDE_FROM_ALL OBJECT ${SOURCE_FILES})

# dependencies
target_link_libraries(${TEST_NAME}
        main
        catch
        ${LIBRARY_NAME}
)