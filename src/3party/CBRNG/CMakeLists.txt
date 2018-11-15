# random number generator library
set(LIBRARY_NAME CBRNG)

# source files
set(SOURCE_FILES 
    ${CMAKE_CURRENT_SOURCE_DIR}/CBRNG_Random.cpp
)

# create static library
add_library(${LIBRARY_NAME} STATIC ${SOURCE_FILES})

# dependencies
target_link_libraries(${LIBRARY_NAME}
        main
        gsl
)

# publish headers for this library
target_include_directories(${LIBRARY_NAME} PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})