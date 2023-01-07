cmake_minimum_required(VERSION 3.22.2)



set(CMAKE_CXX_FLAGS "-std=c++17 ${CMAKE_CXX_FLAGS}")
set(CMAKE_C_COMPILER "/home4/shiweil/software/gcc-12/bin/gcc")
set(CMAKE_CXX_COMPILER "/home4/shiweil/software/gcc-12/bin/g++")

# libcint and lapack 
link_libraries(m)
link_libraries(-lcint)
link_libraries(-llapacke)
link_libraries(-llapack)
link_libraries(-lcblas)
link_libraries(-lrefblas)
link_libraries(-lm)
link_libraries(-lgfortran)
include_directories("/home4/shiweil/clib/libcint/include")
link_directories("/home4/shiweil/clib/libcint/lib64")
include_directories("/home4/shiweil/software/lapack-3.10.1/CBLAS/include")
include_directories("/home4/shiweil/software/lapack-3.10.1/LAPACKE/include")
link_directories("/home4/shiweil/software/lapack-3.10.1")

# xmol
project(xmol)
add_executable(xmol molecule.cpp chemutils.cpp scf.cpp utils.cpp)
