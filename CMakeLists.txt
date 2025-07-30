cmake_minimum_required(VERSION 3.10.0)
project(QEDispInv VERSION 1.0.0 LANGUAGES CXX Fortran)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-non-virtual-dtor")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -fPIC -ffp-contract=fast")
set(CMAKE_Fortran_FLAGS "-Wall -g -O3")

option(ENABLE_DEBUG "Enable debug outputs" OFF)
option(USE_OPENMP "with openmp" ON)

# OpenMP
if (USE_OPENMP)
    find_package(OpenMP REQUIRED)
endif()

# Eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})

# fmt
add_subdirectory(${CMAKE_SOURCE_DIR}/external/fmt EXCLUDE_FROM_ALL)

# HighFive
set(HIGHFIVE_USE_BOOST OFF)
set(HIGHFIVE_BUILD_DOCS OFF)
set(HIGHFIVE_UNIT_TESTS OFF)
set(HIGHFIVE_EXAMPLES OFF)
add_subdirectory(${CMAKE_SOURCE_DIR}/external/highfive)

# toml11
add_subdirectory(${CMAKE_SOURCE_DIR}/external/toml11 EXCLUDE_FROM_ALL)

# LBFGSpp
add_subdirectory(${CMAKE_SOURCE_DIR}/external/LBFGSpp EXCLUDE_FROM_ALL)

include_directories(${CMAKE_SOURCE_DIR}/include)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

add_subdirectory(src)
add_subdirectory(fortran)
