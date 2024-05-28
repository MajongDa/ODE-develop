# cmake_minimum_required(VERSION 3.21)
# project(ODE)


# set(ODE_FILES
#     #boostless multiplication test
#     matrix_vector.cpp
#     matrix_vector_omp.cpp
#     matrix_vector_mpi.cpp
#     #boostless SDDE benchmark
#     bench_sodu.cpp
#     bench_sodu_omp.cpp
#     bench_sodu_mpi.cpp
# )


# add_library ( sodu sodu.cpp sodu.hpp )
# target_link_libraries ( sodu boost_mpi OpenMP::OpenMP_CXX MPI::MPI_CXX)


# foreach(EXE ${ODE_FILES})
#     string(REPLACE ".cpp" "" EXE_NAME ${EXE})
#     add_executable(${EXE_NAME} ${EXE})
#     target_link_libraries(${EXE_NAME} sodu boost_mpi MPI::MPI_CXX OpenMP::OpenMP_CXX)
# endforeach()
cmake_minimum_required(VERSION 3.21)
project(ODE)

find_package(Matplot++ REQUIRED)
find_package(simforode REQUIRED)
find_package(OpenMP 4.5 REQUIRED)
find_package(MPI REQUIRED)

include_directories(${Boost_INCLUDE_DIR})
add_compile_options(-O3)

set(ODE_FILES
    bench_sodu.cpp
    bench_sodu_omp.cpp
    bench_sodu_mpi.cpp
    matrix_vector.cpp
    matrix_vector_omp.cpp
    matrix_vector_mpi.cpp
)

foreach(EXE ${ODE_FILES})
    string(REPLACE ".cpp" "" EXE_NAME ${EXE})
    add_executable(${EXE_NAME} ${EXE})
    target_link_libraries(${EXE_NAME} PUBLIC simforode::simforode Matplot++::matplot MPI::MPI_CXX OpenMP::OpenMP_CXX)
endforeach()
