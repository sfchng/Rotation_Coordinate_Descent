cmake_minimum_required (VERSION 3.0)
project (ROTAVG)

set(default_build_type "Release")

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
    set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
        STRING "Choose the type of build." FORCE)
    # Set the possible values of build type for cmake-gui
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
        "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()


set (CMAKE_CXX_STANDARD 11)


# Make CMake aware of the cmake folder for local FindXXX scripts,
# append rather than set in case the user has passed their own
# additional paths via -D.
list(APPEND CMAKE_MODULE_PATH "${ROTAVG_SOURCE_DIR}/cmake")
include(UpdateCacheVariable)

# Eigen
#find_package (Eigen3 3.3 REQUIRED NO_MODULE)
SET( EIGEN3_INCLUDE_DIR "${ROTAVG_SOURCE_DIR}/3rdparty/Eigen" )
INCLUDE_DIRECTORIES ( "${EIGEN3_INCLUDE_DIR}" )
link_directories ( "${EIGEN3_INCLUDE_DIR}" )
message("${CMAKE_MODULE_PATH}")
message("${EIGEN3_INCLUDE_DIR}")
IF( NOT EIGEN3_INCLUDE_DIR )
    MESSAGE( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
else( EIGEN3_INCLUDE_DIR )
    message("-- Found Eigen3 ${EIGEN3_INCLUDE_DIR}, "
        "building with Eigen3.")
ENDIF()


# SuiteSparse
find_package(SuiteSparse)
if (SUITESPARSE_FOUND)
    # On Ubuntu the system install of SuiteSparse (v3.4.0) up to at least
    # Ubuntu 13.10 cannot be used to link shared libraries.
    if (BUILD_SHARED_LIBS AND
        SUITESPARSE_IS_BROKEN_SHARED_LINKING_UBUNTU_SYSTEM_VERSION)

        message(FATAL_ERROR "You are attempting to build Ceres as a shared "
        "library on Ubuntu using a system package install of SuiteSparse "
        "3.4.0. This package is broken and does not support the "
        "construction of shared libraries (you can still build Ceres as "
        " a static library).  If you wish to build a shared version of Ceres "
        "you should uninstall the system install of SuiteSparse "
        "(libsuitesparse-dev) and perform a source install of SuiteSparse "
        "(we recommend that you use the latest version), "
        "see http://ceres-solver.org/building.html for more information.")
    endif (BUILD_SHARED_LIBS AND
           SUITESPARSE_IS_BROKEN_SHARED_LINKING_UBUNTU_SYSTEM_VERSION)

    message("-- Found SuiteSparse ${SUITESPARSE_VERSION}, "
        "building with SuiteSparse.")
else (SUITESPARSE_FOUND)
    message(FATAL_ERROR "Unable to find SuiteSparse!")
endif (SUITESPARSE_FOUND)

link_directories( ${SUITESPARSE_LIBRARY_DIR_HINTS} )

add_executable( rcd main.cpp rcd.cpp utils.cpp)
add_executable( rcdl main_rcdl.cpp rcdl.cpp)

#message(" SuiteSparse INCLUDE ${SUITESPARSE_INCLUDE_DIRS}, "
#"building with SuiteSparse.")


target_include_directories( rcd PUBLIC
 ${EIGEN3_INCLUDE_DIR}
 ${SUITESPARSE_INCLUDE_DIRS}
)

target_include_directories( rcdl PUBLIC
 ${EIGEN3_INCLUDE_DIR}
 ${SUITESPARSE_INCLUDE_DIRS}
)


target_link_libraries (rcd ${EIGEN3_LIBRARIES} ${SUITESPARSE_LIBRARIES} umfpack)
target_link_libraries (rcdl ${EIGEN3_LIBRARIES} ${SUITESPARSE_LIBRARIES} umfpack)