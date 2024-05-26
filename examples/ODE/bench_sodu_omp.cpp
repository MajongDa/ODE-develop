#include <iostream>
#include "sodu.hpp"

float koef_matrix_settings ( unsigned j, unsigned i, unsigned n )
    {
    float value;
    if ( j==i ) value = -4.f;
    else if ( j>0 && ( ( j-1 ) == i ) ) value = 1.f;
    else if ( j< ( n-1 ) && ( ( j+1 ) == i ) ) value = 1.f;
    else value = 0.0f;
    return value;
    }

simfor::matrix<float> odu_matrix_create ( unsigned n )
    {
    simfor::matrix<float> M ( n, n );
    for ( int i = 0; i < n; i++ )
        for ( int j = 0; j < n; j++ )
            M ( i, j ) = koef_matrix_settings ( j,i,n );
    return M;
    }

float integral_function ( float x, int number )
    {
    return exp ( -1 * powf ( x - number, 2 ) );
    }

float integral_result ( float x, int number )
    {
    float integral = 0;
    float gauss_x[7] = {-1.f,
                        - sqrtf ( 5. / 11. + 2. / 11. * sqrtf ( 5. / 3. ) ),
                        - sqrtf ( 5. / 11. - 2. / 11. * sqrtf ( 5. / 3. ) ), 0.0,
                        sqrtf ( 5. / 11. - 2. / 11. * sqrtf ( 5. / 3. ) ),
                        sqrtf ( 5. / 11. + 2. / 11. * sqrtf ( 5. / 3. ) ), 1.0
                       },
                       gauss_weight[7] = {1. / 21.,
                                          ( 124 - 7 * sqrtf ( 15 ) ) / 350.0f,
                                          ( 124 + 7 * sqrtf ( 15 ) ) / 350.0f,
                                          256. / 525.,
                                          ( 124 + 7 * sqrtf ( 15 ) ) / 350.0f,
                                          ( 124 - 7 * sqrtf ( 15 ) ) / 350.0f,
                                          1. / 21.
                                         };

    for ( int i = 0; i < 7; ++i )
        integral +=  gauss_weight[i] *
                     integral_function ( x / 2. * gauss_x[i] + x / 2., number );
    integral = integral * x / 2.;

    return integral;
    }

simfor::vector<float> initial_cond ( float a, int n )
    {
    simfor::vector<float> y0 ( n );
    for ( unsigned i = 0; i < n; ++i )
        y0 ( i ) = integral_result ( a, i );
    return y0;
    }

int main ( int argc, char* argv [] )
    {
    int n = atoi ( argv [ 1 ] ), p = atoi ( argv [ 2 ] );
    float a=0, b=1, h;
    h = ( b - a ) / n;
    simfor::matrix<float> F = odu_matrix_create ( n ), Y;
    simfor::vector<float> x0 = initial_cond ( 1, n );

    omp_set_num_threads ( p );

    //burnin
    for ( int i = 0; i < 2; ++i )
        Y = simfor::eiler_system_solve_matrix_omp ( h, n, x0, F );

    std::cout << "Threads: " << p << "\n";

    double t;

    t = omp_get_wtime();
    Y = simfor::eiler_system_solve_matrix_omp ( h, n, x0, F );
    t = ( omp_get_wtime() - t );
    std::cout << "E_OMP " << t  << "\n";

    t = omp_get_wtime();
    Y = simfor::rk_system_solve_matrix_omp ( h, n, x0, F );
    t = ( omp_get_wtime() - t );
    std::cout << "RK_OMP " << t  << "\n";

    t = omp_get_wtime();
    Y = simfor::adams5_system_solve_matrix_omp ( h, n, x0, F );
    t = ( omp_get_wtime() - t );
    std::cout << "ADB_OMP " << t <<  "\n";

    // double t, avg_t = 0, min_t = 1000;
    // for ( int i = 0; i < 10; ++i )
    //     {
    //     t = omp_get_wtime();
    //     Y = simfor::eiler_system_solve_matrix_omp ( h, n, x0, F );
    //     avg_t += t = ( omp_get_wtime() - t );
    //     min_t = min_t > t ? t : min_t;
    //     }
    // std::cout << "E_OMP avg " << avg_t/10  << " min_time "<< min_t << "\n";
    // min_t = 1000; avg_t = 0;
    // for ( int i = 0; i < 10; ++i )
    //     {
    //     t = omp_get_wtime();
    //     Y = simfor::rk_system_solve_matrix_omp ( h, n, x0, F );
    //     avg_t += t = ( omp_get_wtime() - t );
    //     min_t = min_t > t ? t : min_t;
    //     }
    // std::cout << "RK_OMP avg " << avg_t/10  << " min_time "<< min_t << "\n";
    // min_t = 1000; avg_t = 0;
    // for ( int i = 0; i < 10; ++i )
    //     {
    //     t = omp_get_wtime();
    //     Y = simfor::adams5_system_solve_matrix_omp ( h, n, x0, F );
    //     avg_t += t = ( omp_get_wtime() - t );
    //     min_t = min_t > t ? t : min_t;
    //     }
    //     std::cout << "ADB_OMP avg " << avg_t/10 << " min_time "<< min_t <<  "\n";
    }

