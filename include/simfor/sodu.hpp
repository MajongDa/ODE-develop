
#ifndef SIMFOR_SODU_HPP
#define SIMFOR_SODU_HPP

#include <omp.h>
#include <iostream>
#include "multMatrixVector.hpp"


namespace simfor
{
//коэффаценты соотв. методов
float adams_bashford_koeff( int );
float adams_moulton_koeff( int );

//метотоды Адамса-Милна
matrix<float> AdMiln_system_solve_matrix ( float, int, vector<float>, matrix<float> );
matrix<float> AdMiln_system_solve_matrix_omp ( float, int, vector<float>, matrix<float> );
matrix<float> AdMiln_system_solve_matrix_mpi ( float, int, vector<float>, matrix<float> );

//метотоды Адамса-БашфортаW
matrix<float> adams5_system_solve_matrix ( float, int, vector<float>, matrix<float> );
matrix<float> adams5_system_solve_matrix_omp ( float, int, vector<float>, matrix<float> );
matrix<float> adams5_system_solve_matrix_mpi ( float, int, vector<float>, matrix<float> );

//метотоды Адамса-Моултона

matrix<float> AdMltn_system_solve_matrix ( float, int, vector<float>, matrix<float> );
matrix<float> AdMltn_system_solve_matrix_omp ( float, int, vector<float>, matrix<float> );
matrix<float> AdMltn_system_solve_matrix_mpi ( float, int, vector<float>, matrix<float> );

//eiler
matrix<float> eiler_system_solve_matrix ( float, int, vector<float>, matrix<float>);
matrix<float> eiler_system_solve_matrix_omp ( float, int, vector<float>, matrix<float>);
matrix<float> eiler_system_solve_matrix_mpi ( float, int, vector<float>, matrix<float>);

//rk(Runge-Kutta)
matrix<float> rk_system_solve_matrix ( float, int, vector<float>, matrix<float>);
matrix<float> rk_system_solve_matrix_omp ( float, int, vector<float>, matrix<float>);
matrix<float> rk_system_solve_matrix_mpi ( float, int, vector<float>, matrix<float>);
}

#endif
