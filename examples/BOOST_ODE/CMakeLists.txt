cmake_minimum_required(VERSION 3.21)
project(BOOST_ODE)

find_package(Matplot++ REQUIRED)
find_package(simforode REQUIRED)
find_package(Boost 1.83.0 REQUIRED)
find_package(OpenMP 4.5 REQUIRED)
find_package(MPI REQUIRED)

set(ODE_FILES
    odu_func.cpp
    sodu_func.cpp
    boost_bench_sodu.cpp
    boost_bench_sodu_omp.cpp
    boost_bench_sodu_mpi.cpp

    boost_matr_vec.cpp
    boost_matr_vec_omp.cpp
    boost_matr_vec_mpi.cpp
)

foreach(EXE ${ODE_FILES})
    string(REPLACE ".cpp" "" EXE_NAME ${EXE})
    add_executable(${EXE_NAME} ${EXE})
    target_link_libraries(${EXE_NAME} PUBLIC simforode::simforode Matplot++::matplot)
endforeach()
