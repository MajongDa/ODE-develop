#ifndef SIMFOR_ODU_HPP
#define SIMFOR_ODU_HPP

#include <omp.h>
#include <boost/mpi.hpp>
#include "internal/types.hpp"
#include "simfor/multMatrVec.hpp"
#include <chrono>
#include <iostream>

namespace mpi = boost::mpi;

namespace simfor
{
//коэффаценты соотв. методов
float adams_bashford_koeff( int );
float adams_moulton_koeff( int );

//метотоды Адамса-Милна
vec AdMiln_function_solve ( float, float, float, int, float ( * ) ( float, float ), float );
matr AdMiln_function_solve_vector ( float, float, float, vec ( * ) ( float, vec), vec );
matr AdMiln_system_solve_matrix ( float, int, vec, matr );
matr AdMiln_system_solve_matrix_omp ( float, int, vec, matr );
matr AdMiln_system_solve_matrix_mpi ( float, int, vec, matr );

//метотоды Адамса-Башфорта
vec adams5_function_solve ( float, float, float, int, float ( * ) ( float, float ), float);
matr adams5_function_solve_vector ( float, float, float, int, vec ( * ) ( float, vec ), vec );
matr adams5_system_solve_matrix ( float, int, vec, matr );
matr adams5_system_solve_matrix_omp ( float, int, vec, matr );
matr adams5_system_solve_matrix_mpi ( float, int, vec, matr );

//метотоды Адамса-Моултона
vec AdMltn_function_solve ( float, float, float, int, float ( * ) ( float, float ), float);
matr AdMltn_function_solve_vector ( float, float, float, int, vec ( * ) ( float, vec), vec );
matr AdMltn_system_solve_matrix ( float, int, vec, matr );
matr AdMltn_system_solve_matrix_omp ( float, int, vec, matr );
matr AdMltn_system_solve_matrix_mpi ( float, int, vec, matr );

//eiler
vec eiler_function_solve ( float, float, float, int, float ( * ) ( float, float ), float );
matr eiler_system_solve_vector ( float, float, float, int, vec ( * ) ( float, vec ), vec );
matr eiler_system_solve_matrix ( float, int, vec, matr);
matr eiler_system_solve_matrix_omp ( float, int, vec, matr);
matr eiler_system_solve_matrix_mpi ( float, int, vec, matr);

//rk(Runge-Kutta)
vec rk_function_solve ( float, float, float, int, float ( * ) ( float, float ), float );
matr rk_system_function_solve_vector ( float, float, float, int, vec( * ) ( float, vec), vec );
matr rk_system_solve_matrix ( float, int, vec, matr);
matr rk_system_solve_matrix_omp ( float, int, vec, matr);
matr rk_system_solve_matrix_mpi ( float, int, vec, matr);
}

#endif
