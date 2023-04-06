# FOR LAMMPS CMAKE BUILD

# set installation location
set(CMAKE_INSTALL_PREFIX "$ENV{PREFIX}" CACHE BOOL "" FORCE)

#set(Python_FIND_FRAMEWORK LAST)

find_package(Python 3.11 EXACT REQUIRED)

# enforce c++11 standards
set(CCFLAGS -g -O3 -std=c++11)

# compile a binary and shared library
set(BUILD_SHARED_LIBS ON CACHE BOOL "" FORCE)

# allow error messages (very useful)
set(LAMMPS_EXCEPTIONS ON CACHE BOOL "" FORCE)

# minimal packages to run example (MANYBODY, EXTRA-FIX) and generate new pathway (REPLICA for "fix neb")
set(ALL_PACKAGES MANYBODY EXTRA-FIX REPLICA ML-SNAP)

foreach(PKG ${ALL_PACKAGES})
  set(PKG_${PKG} ON CACHE BOOL "" FORCE)
endforeach()
