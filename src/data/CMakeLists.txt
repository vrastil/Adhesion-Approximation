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