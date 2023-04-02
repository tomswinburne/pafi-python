
cmake_minimum_required(VERSION 3.5)

set(DEBUG 0)

set(PREFIX "$ENV{PREFIX}") # example value
set(CMAKE_CXX_COMPILER /usr/local/bin/mpicxx) # example value
set(CMAKE_C_COMPILER /usr/local/bin/mpicc) # example value
set(CMAKE_INCLUDE_PATH ${PREFIX}/include) # location of copied lammps header directory (lammps/*.h)
set(CMAKE_LIBRARY_PATH ${PREFIX}/lib) # location of liblammps_mpi.a


enable_testing()

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules)

# Enable C++ and ensure C++11
enable_language(CXX)
set(CMAKE_CXX_STANDARD 11) # C++11...
set(CMAKE_CXX_STANDARD_REQUIRED ON) #...is required...
set(CMAKE_CXX_EXTENSIONS ON) #...without compiler extensions like gnu++11


project (pafi)

include_directories(core)
include_directories(core/rapidxml)

include_directories(${CMAKE_INCLUDE_PATH})
link_directories(${CMAKE_LIBRARY_PATH})


find_package(LAMMPS REQUIRED)
include_directories(${LAMMPS_INCLUDE_DIRS})


# include VERBOSE for debugging output
add_definitions(-std=c++11 -O3 -DLAMMPS_EXCEPTIONS -DLAMMPS_LIB_MPI)
if(DEBUG)
	add_definitions(-g -DVERBOSE)
endif(DEBUG)

file(GLOB src_files "core/*")

add_executable(pafi ${src_files} pafi-pafi/standard/pafi.cpp)

add_executable(pafi-lammps-preparation ${src_files} path-tools/pafi-lammps-preparation.cpp)
add_executable(pafi-path-test ${src_files} path-tools/pafi-path-test.cpp)

target_link_libraries(pafi lammps_mpi)
target_link_libraries(pafi-lammps-preparation lammps_mpi)
target_link_libraries(pafi-path-test lammps_mpi)

# custom pafi example
add_executable(pafi-custom ${src_files} pafi-pafi/custom/pafi.cpp)
target_link_libraries(pafi-custom lammps_mpi)
