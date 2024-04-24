cmake_minimum_required(VERSION 3.21)
project(ODEtest)


set(ODE_X_FILES
    odu_func.cpp
    sodu_test.cpp
    multVecMatrTest.cpp
    bench_sodu.cpp
    bench_sodu_omp.cpp
    bench_sodu_mpi.cpp
)

foreach(EXE ${ODE_X_FILES})
    string(REPLACE ".cpp" "" EXE_NAME ${EXE})
    add_executable(${EXE_NAME} ${EXE})
    target_link_libraries(${EXE_NAME} PUBLIC simfor)
endforeach()
