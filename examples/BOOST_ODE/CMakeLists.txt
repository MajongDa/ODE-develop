cmake_minimum_required(VERSION 3.21)
project(BOOST_ODE)

find_package(Matplot++ REQUIRED)

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
    target_link_libraries(${EXE_NAME} PUBLIC simforode Matplot++::matplot)
endforeach()
