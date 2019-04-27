# check new policy (from cmake 3.13), use legacy behavior
if(POLICY CMP0078)
  cmake_policy(SET CMP0078 OLD)
else()
  set(UseSWIG_TARGET_NAME_PREFERENCE LEGACY)
endif()

# check new policy (from cmake 3.14), use legacy behavior
if(POLICY CMP0086)
  cmake_policy(SET CMP0086 OLD)
endif()

# FastSim python wrapper
set(LIBRARY_NAME fastsim)

# source files
set(SOURCE_FILE ${CMAKE_CURRENT_SOURCE_DIR}/all.i)

# ???
include(UseSWIG)
include_directories(BEFORE ${PYTHON_INCLUDE_DIRS})
SET(CMAKE_SWIG_FLAGS "")

# identify source file as C++
set_property(SOURCE ${SOURCE_FILE} PROPERTY CPLUSPLUS ON)

# get dependency include directories
set(DEP_LIBS core data ccl main)
foreach(LIB ${DEP_LIBS})
    get_property(TMP TARGET ${LIB} PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
    set(LIB_INC ${LIB_INC} -I${TMP})
endforeach(LIB)

# set include directories
set_property(SOURCE ${SOURCE_FILE} PROPERTY SWIG_FLAGS ${LIB_INC} "-DPRECISION=2")

# create wrapper
if(${CMAKE_VERSION} VERSION_LESS "3.8.0")
  swig_add_module(${LIBRARY_NAME} python ${SOURCE_FILE})
else()
  swig_add_library(${LIBRARY_NAME}
                   TYPE SHARED
                   LANGUAGE python
                   SOURCES ${SOURCE_FILE})
endif()

# dependencies
swig_link_libraries(${LIBRARY_NAME}
        core data ccl main ${PYTHON_LIBRARIES}
)

# install python module 
install(TARGETS ${SWIG_MODULE_${LIBRARY_NAME}_REAL_NAME}
        DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../simpy/
)
install(FILES ${CMAKE_BINARY_DIR}/swig/${LIBRARY_NAME}.py
        DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../simpy/)
