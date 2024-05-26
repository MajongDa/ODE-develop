#ifndef SIMFOR_ODU_HPP
#define SIMFOR_ODU_HPP

#include <omp.h>
#include <boost/mpi.hpp>
#include "internal/types.hpp"
#include "multMatrVec.hpp"
#include <iostream>



namespace simfor
{
namespace mpi = boost::mpi;
//Coeffisients for Adams methods
double adams_bashford_koeff( int );
double adams_moulton_koeff( int );

//Adams-Milne method
vec AdMiln_function_solve ( double, double, double, int, double ( * ) ( double, double ), double );
matr AdMiln_function_solve_vector ( double, double, double, vec ( * ) ( double, vec), vec );

//метотоды Адамса-Башфорта
vec adams5_function_solve ( double, double, double, int, double ( * ) ( double, double ), double);
matr adams5_function_solve_vector ( double, double, double, int, vec ( * ) ( double, vec ), vec );
matr adams5_system_solve_matrix ( double, int, vec, matr );
matr adams5_system_solve_matrix_omp ( double, int, vec, matr );
matr adams5_system_solve_matrix_mpi ( double, int, vec, matr );

//Adams-Moulton method
vec AdMltn_function_solve ( double, double, double, int, double ( * ) ( double, double ), double);
matr AdMltn_function_solve_vector ( double, double, double, int, vec ( * ) ( double, vec), vec );
matr AdMltn_system_solve_matrix ( double, int, vec, matr );
matr AdMltn_system_solve_matrix_omp ( double, int, vec, matr );
matr AdMltn_system_solve_matrix_mpi ( double, int, vec, matr );

//Euler method
vec eiler_function_solve ( double, double, double, int, double ( * ) ( double, double ), double );
matr eiler_system_solve_vector ( double, double, double, int, vec ( * ) ( double, vec ), vec );
matr eiler_system_solve_matrix ( double, int, vec, matr);
matr eiler_system_solve_matrix_omp ( double, int, vec, matr);
matr eiler_system_solve_matrix_mpi ( double, int, vec, matr);

//Runge-Kutta classical 4 stage method
vec rk4_function_solve ( double, double, double, int, double ( * ) ( double, double ), double );
matr rk4_system_function_solve_vector ( double, double, double, int, vec( * ) ( double, vec), vec );
matr rk4_system_solve_matrix ( double, int, vec, matr);
matr rk4_system_solve_matrix_omp ( double, int, vec, matr);
matr rk4_system_solve_matrix_mpi ( double, int, vec, matr);

//Runge-Kutta classical 6 stage method
vec rk6_function_solve ( double, double, double, int, double ( * ) ( double, double ), double );
matr rk6_system_function_solve_vector ( double, double, double, int, vec( * ) ( double, vec), vec );
matr rk6_system_solve_matrix ( double, int, vec, matr);
matr rk6_system_solve_matrix_omp ( double, int, vec, matr);
matr rk6_system_solve_matrix_mpi ( double, int, vec, matr);

//Runge-Kutta Dormand-Prince method
matr rk54( double, double, double ( * ) ( double, double ), double , double, double, double);
matr rk54( double, double, vec( * ) ( double, vec), vec, double , double, double );

//Runge-Kutta Fahlberg method
matr rk45( double, double, double ( * ) ( double, double ), double , double, double, double);
matr rk45( double, double, vec( * ) ( double, vec), vec, double , double, double );
}

#endif
