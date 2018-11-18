# input/output library
set(LIBRARY_NAME io)

# source files
set(SOURCE_FILES 
    ${CMAKE_CURRENT_SOURCE_DIR}/core_cmd.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/core_out.cpp
)

# create static library
add_library(${LIBRARY_NAME} STATIC ${SOURCE_FILES})

# dependencies
target_link_libraries(${LIBRARY_NAME}
        main
        data
        core
        boost_program_options
        boost_filesystem
        boost_system
)

# publish headers for this library
target_include_directories(${LIBRARY_NAME} PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include")