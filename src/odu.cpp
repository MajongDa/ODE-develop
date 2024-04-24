#include "odu.hpp"

namespace simfor
{

float adams_bashford_koeff ( int i )
    {
    float a[5] = {1901. / 720., 1387. / 360., 109. / 30., 637. / 360., 251. / 720.};
    return a[i];
    }
float adams_moulton_koeff ( int i )
    {
    float a[5] = {251. / 720., 646. / 720., 264. / 720., 106. / 720., 19. / 720.};
    return a[i];
    }

vec AdMiln_function_solve ( float a, float b, float h, int n, float ( *f ) ( float, float ), float y0 )
    {
    float xp;
    vec y ( n + 1 );

    //start step
    subrange ( y, 0, 4 ) = rk_function_solve ( a, b, h, 3, f, y0 );

    for ( int i = 3; i < n; ++i )
        {
        xp = y ( i ) + h / 24 * (
                 -9 * f ( a + h * ( i - 3 ), y ( i - 3 ) ) +
                 37 * f ( a + h * ( i - 2 ), y ( i - 2 ) ) -
                 59 * f ( a + h * ( i - 1 ), y ( i - 1 ) ) +
                 55 * f ( a + h * ( i ), y ( i ) )
             );

        y ( i + 1 ) = y ( i ) + h / 24 * (
                          f ( a + h * ( i - 2 ), y ( i - 2 ) ) -
                          5 * f ( a + h * ( i - 1 ), y ( i - 1 ) ) +
                          19 * f ( a + h * ( i ), y ( i ) ) +
                          9 * f ( a + h * ( i + 1 ), xp )
                      );
        }
    return y;
    }
//Метод Адамса-Башфорда
vec adams5_function_solve ( float a, float b, float h, int n, float ( *f ) ( float, float ), float y0 )
    {
    vec y ( n + 1 );

    subrange ( y, 0, 5 ) = rk_function_solve ( a, b, h, 4, f, y0 );

    for ( int i = 4; i < n; ++i )
        {

        y ( i + 1 ) = y ( i ) + h * (
                          adams_bashford_koeff ( 0 ) * f ( a + h * ( i ),     y ( i ) ) -
                          adams_bashford_koeff ( 1 ) * f ( a + h * ( i - 1 ), y ( i - 1 ) ) +
                          adams_bashford_koeff ( 2 ) * f ( a + h * ( i - 2 ), y ( i - 2 ) ) -
                          adams_bashford_koeff ( 3 ) * f ( a + h * ( i - 3 ), y ( i - 3 ) ) +
                          adams_bashford_koeff ( 4 ) * f ( a + h * ( i - 4 ), y ( i - 4 ) )
                      );
        }
    return y;
    }

matr adams5_function_solve_vector ( float a, float b, float h, int n, vec ( *vf ) ( float, vec ), vec vs )
    {
    if ( n < 5 ) return matr();

    matr Y ( n + 1, vs.size() );

    subrange ( Y, 0, 5, 0, vs.size() ) =
    rk_system_function_solve_vector ( a, b, h, 4, vf, vs );

    for ( int i = 4; i < n; ++i )
        {
        row ( Y, i + 1 ) = row ( Y, i ) +
                           h * ( vf ( a + h * ( i ),     row ( Y, i ) )     * adams_bashford_koeff ( 0 ) -
                                 vf ( a + h * ( i - 1 ), row ( Y, i - 1 ) ) * adams_bashford_koeff ( 1 ) +
                                 vf ( a + h * ( i - 2 ), row ( Y, i - 2 ) ) * adams_bashford_koeff ( 2 ) -
                                 vf ( a + h * ( i - 3 ), row ( Y, i - 3 ) ) * adams_bashford_koeff ( 3 ) +
                                 vf ( a + h * ( i - 4 ), row ( Y, i - 4 ) ) * adams_bashford_koeff ( 4 )
                               );
        }
    return Y;
    }


matr adams5_system_solve_matrix ( float h, int n, vec y_0, matr F )
    {
    vec k1, k2, k3, k4;
    matr Y ( n + 1, y_0.size() );

    subrange ( Y, 0, 5, 0, y_0.size() ) =
        rk_system_solve_matrix ( h, 4, y_0, F );


//     row ( Y, 0 ) =  y_0;
//     for ( int i = 0; i < 4; ++i )
//         {
//         k1 = prod ( F, row ( Y, i ) );
//         k2 = prod ( F, row ( Y, i ) + h / 2. * k1 );
//         k3 = prod ( F, row ( Y, i ) + h / 2. * k2 );
//         k4 = prod ( F, row ( Y, i ) + h      * k3 );
//         row ( Y, i + 1 ) = row ( Y, i ) + h / 6.* ( k1 + 2 * k2 + 2 * k3 + k4 );
//         }

    for ( int i = 4; i < n; ++i )
        {
        row ( Y, i + 1 ) = row ( Y, i ) + h * (
                               adams_bashford_koeff ( 0 ) * prod ( F, row ( Y,  i ) ) -
                               adams_bashford_koeff ( 1 ) * prod ( F, row ( Y, i - 1 ) ) +
                               adams_bashford_koeff ( 2 ) * prod ( F, row ( Y, i - 2 ) ) -
                               adams_bashford_koeff ( 3 ) * prod ( F, row ( Y, i - 3 ) ) +
                               adams_bashford_koeff ( 4 ) * prod ( F, row ( Y, i - 4 ) )
                           );

        }


    return Y;
    }

//TODO add omp for scalar-vector mult
matr adams5_system_solve_matrix_omp ( float h, int n, vec y_0, matr F )
    {
    vec k1, k2, k3, k4;
    matr Y ( n + 1, y_0.size() );

    subrange ( Y, 0, 5, 0, y_0.size() ) =
        rk_system_solve_matrix_omp ( h, 4, y_0, F );


//     row ( Y, 0 ) =  y_0;
//
//     for ( int i = 0; i < 4; ++i )
//         {
//         k1 = multMatrVec_omp ( F, row ( Y, i ) );
//         k2 = multMatrVec_omp ( F, row ( Y, i ) + h / 2. * k1 );
//         k3 = multMatrVec_omp ( F, row ( Y, i ) + h / 2. * k2 );
//         k4 = multMatrVec_omp ( F, row ( Y, i ) + h      * k3 );
//         row ( Y, i + 1 ) = row ( Y, i ) + h / 6.* ( k1 + 2 * k2 + 2 * k3 + k4 );
//         }

    for ( int i = 4; i < n; ++i )
        {
        row ( Y, i + 1 ) = row ( Y, i ) + h * (
                                 adams_bashford_koeff ( 0 ) * multMatrVec_omp ( F, row ( Y,  i ) ) -
                                 adams_bashford_koeff ( 1 ) * multMatrVec_omp ( F, row ( Y, i - 1 ) ) +
                                 adams_bashford_koeff ( 2 ) * multMatrVec_omp ( F, row ( Y, i - 2 ) ) -
                                 adams_bashford_koeff ( 3 ) * multMatrVec_omp ( F, row ( Y, i - 3 ) ) +
                                 adams_bashford_koeff ( 4 ) * multMatrVec_omp ( F, row ( Y, i - 4 ) )
                             );
        }
    return Y;
    }

//TODO add omp for scalar-vector mult
matr adams5_system_solve_matrix_mpi ( float h, int n, vec y_0, matr F )
    {
    vec k1, k2, k3, k4;
    matr Y ( n + 1, y_0.size() );


    subrange ( Y, 0, 5, 0, y_0.size() ) =
        rk_system_solve_matrix_mpi ( h, 4, y_0, F );

//     for ( int i = 0; i < 4; ++i )
//         {
//         k1 = multMatrVec_mpi ( F, row ( Y, i ) );
//         k2 = multMatrVec_mpi ( F, row ( Y, i ) + h / 2. * k1 );
//         k3 = multMatrVec_mpi ( F, row ( Y, i ) + h / 2. * k2 );
//         k4 = multMatrVec_mpi ( F, row ( Y, i ) + h      * k3 );
//         row ( Y, i + 1 ) = row ( Y, i ) + h / 6.* ( k1 + 2 * k2 + 2 * k3 + k4 );
//         }

    for ( int i = 4; i < n; ++i )
        {
        row ( Y, i + 1 ) = row ( Y, i ) + h * (
                                 adams_bashford_koeff ( 0 ) * multMatrVec_mpi ( F, row ( Y,  i ) ) -
                                 adams_bashford_koeff ( 1 ) * multMatrVec_mpi ( F, row ( Y, i - 1 ) ) +
                                 adams_bashford_koeff ( 2 ) * multMatrVec_mpi ( F, row ( Y, i - 2 ) ) -
                                 adams_bashford_koeff ( 3 ) * multMatrVec_mpi ( F, row ( Y, i - 3 ) ) +
                                 adams_bashford_koeff ( 4 ) * multMatrVec_mpi ( F, row ( Y, i - 4 ) )
                             );
        }


    return Y;
    }

//Метод Адамса-Моултона
vec AdMltn_function_solve ( float a, float b, float h, int n, float ( *f ) ( float, float ), float y0 )
    {
    vec y ( n + 1 );
    float k1;

    subrange ( y, 0, 5 ) = rk_function_solve ( a, b, h, 4, f, y0 );

    for ( int i = 4; i < n; ++i )
        {

        k1 = y ( i ) + h * (
                 adams_bashford_koeff ( 0 ) * f ( a + h * ( i ),     y ( i ) ) -
                 adams_bashford_koeff ( 1 ) * f ( a + h * ( i - 1 ), y ( i - 1 ) ) +
                 adams_bashford_koeff ( 2 ) * f ( a + h * ( i - 2 ), y ( i - 2 ) ) -
                 adams_bashford_koeff ( 3 ) * f ( a + h * ( i - 3 ), y ( i - 3 ) ) +
                 adams_bashford_koeff ( 4 ) * f ( a + h * ( i - 4 ), y ( i - 4 ) )
             );
        k1 = y ( i ) + h * (
                          adams_moulton_koeff ( 0 ) * f ( a + h * ( i + 1 ),     k1 ) +
                          adams_moulton_koeff ( 1 ) * f ( a + h * ( i ), y ( i ) ) -
                          adams_moulton_koeff ( 2 ) * f ( a + h * ( i - 1 ), y ( i - 1 ) ) +
                          adams_moulton_koeff ( 3 ) * f ( a + h * ( i - 2 ), y ( i - 2 ) ) -
                          adams_moulton_koeff ( 4 ) * f ( a + h * ( i - 3 ), y ( i - 3 ) )
                      );
        y ( i + 1 ) = y ( i ) + h * (
                          adams_moulton_koeff ( 0 ) * f ( a + h * ( i + 1 ),     k1 ) +
                          adams_moulton_koeff ( 1 ) * f ( a + h * ( i ), y ( i ) ) -
                          adams_moulton_koeff ( 2 ) * f ( a + h * ( i - 1 ), y ( i - 1 ) ) +
                          adams_moulton_koeff ( 3 ) * f ( a + h * ( i - 2 ), y ( i - 2 ) ) -
                          adams_moulton_koeff ( 4 ) * f ( a + h * ( i - 3 ), y ( i - 3 ) )
                      );
        }
    return y;
    }

matr AdMltn_function_solve_vector ( float a, float b, float h, int n, vec ( *vf ) ( float, vec ), vec vs )
    {

    if ( n < 5 ) return matr();

    vec P;

    matr Y ( ( n + 1 ), vs.size() );

    subrange ( Y, 0, 5, 0, vs.size() ) =
    rk_system_function_solve_vector ( a, b, h, 4, vf, vs );

    for ( int i = 4; i < n; i++ )
        {
        row ( Y, i + 1 ) = row ( Y, i ) +
                           h * ( vf ( a + h * ( i ),     row ( Y, i ) )     * adams_bashford_koeff ( 0 ) -
                                 vf ( a + h * ( i - 1 ), row ( Y, i - 1 ) ) * adams_bashford_koeff ( 1 ) +
                                 vf ( a + h * ( i - 2 ), row ( Y, i - 2 ) ) * adams_bashford_koeff ( 2 ) -
                                 vf ( a + h * ( i - 3 ), row ( Y, i - 3 ) ) * adams_bashford_koeff ( 3 ) +
                                 vf ( a + h * ( i - 4 ), row ( Y, i - 4 ) ) * adams_bashford_koeff ( 4 )
                               );

        P = row ( Y, i ) + h * (
                 adams_bashford_koeff ( 0 ) * vf ( a + h * ( i ),     row ( Y, i ) ) -
                 adams_bashford_koeff ( 1 ) * vf ( a + h * ( i - 1 ), row ( Y, i - 1 ) ) +
                 adams_bashford_koeff ( 2 ) * vf ( a + h * ( i - 2 ), row ( Y, i - 2 ) ) -
                 adams_bashford_koeff ( 3 ) * vf ( a + h * ( i - 3 ), row ( Y, i - 3 ) ) +
                 adams_bashford_koeff ( 4 ) * vf ( a + h * ( i - 4 ), row ( Y, i - 4 ) )
             );
        row ( Y, i + 1 ) = row ( Y, i ) + h * (
                               adams_moulton_koeff ( 0 ) * vf ( a + h * ( i + 1 ),     P ) -
                               adams_moulton_koeff ( 1 ) * vf ( a + h * ( i ), row ( Y, i ) ) +
                               adams_moulton_koeff ( 2 ) * vf ( a + h * ( i - 1 ), row ( Y, i - 1 ) ) +
                               adams_moulton_koeff ( 3 ) * vf ( a + h * ( i - 2 ), row ( Y, i - 2 ) ) +
                               adams_moulton_koeff ( 4 ) * vf ( a + h * ( i - 3 ), row ( Y, i - 3 ) )
                           );
        }



    return Y;
    }


matr AdMltn_system_solve_matrix ( float h, int n, vec y_0, matr F )
    {
    vec k1, k2, k3, k4;
    matr Y ( n + 1, y_0.size() );

    row ( Y, 0 ) =  y_0;

    for ( int i = 0; i < 4; ++i )
        {
        k1 = prod ( F, row ( Y, i ) );
        k2 = prod ( F, row ( Y, i ) + h / 2. * k1 );
        k3 = prod ( F, row ( Y, i ) + h / 2. * k2 );
        k4 = prod ( F, row ( Y, i ) + h      * k3 );
        row ( Y, i + 1 ) = row ( Y, i ) + h / 6.* ( k1 + 2 * k2 + 2 * k3 + k4 );
        }

    for ( int i = 4; i < n; ++i )
        {
        k1 = row ( Y, i ) + h * (
                 adams_bashford_koeff ( 0 ) * prod ( F, row ( Y,  i ) ) -
                 adams_bashford_koeff ( 1 ) * prod ( F, row ( Y, i - 1 ) ) +
                 adams_bashford_koeff ( 2 ) * prod ( F, row ( Y, i - 2 ) ) -
                 adams_bashford_koeff ( 3 ) * prod ( F, row ( Y, i - 3 ) ) +
                 adams_bashford_koeff ( 4 ) * prod ( F, row ( Y, i - 4 ) )
             );

        row ( Y, i + 1 ) = row ( Y, i ) + h* (
                                 adams_moulton_koeff ( 0 ) * prod ( F, k1 ) -
                                 adams_moulton_koeff ( 1 ) * prod ( F, row ( Y, i ) ) +
                                 adams_moulton_koeff ( 2 ) * prod ( F, row ( Y, i - 1 ) ) -
                                 adams_moulton_koeff ( 3 ) * prod ( F, row ( Y, i - 2 ) ) +
                                 adams_moulton_koeff ( 4 ) * prod ( F, row ( Y, i - 3 ) )
                             );
        }

    return Y;
    }

matr AdMltn_system_solve_matrix_omp ( float h, int n, vec y_0, matr F )
    {
    vec k1, k2, k3, k4;
    matr Y ( n + 1, y_0.size() );

    row ( Y, 0 ) =  y_0;

    for ( int i = 0; i < 4; ++i )
        {
        k1 = multMatrVec_omp ( F, row ( Y, i ) );
        k2 = multMatrVec_omp ( F, row ( Y, i ) + h / 2. * k1 );
        k3 = multMatrVec_omp ( F, row ( Y, i ) + h / 2. * k2 );
        k4 = multMatrVec_omp ( F, row ( Y, i ) + h      * k3 );
        row ( Y, i + 1 ) = row ( Y, i ) + h / 6.* ( k1 + 2 * k2 + 2 * k3 + k4 );
        }

    for ( int i = 4; i < n; ++i )
        {
        k1 = row ( Y, i ) + h * (
                 adams_bashford_koeff ( 0 ) * multMatrVec_omp ( F, row ( Y,  i ) ) -
                 adams_bashford_koeff ( 1 ) * multMatrVec_omp ( F, row ( Y, i - 1 ) ) +
                 adams_bashford_koeff ( 2 ) * multMatrVec_omp ( F, row ( Y, i - 2 ) ) -
                 adams_bashford_koeff ( 3 ) * multMatrVec_omp ( F, row ( Y, i - 3 ) ) +
                 adams_bashford_koeff ( 4 ) * multMatrVec_omp ( F, row ( Y, i - 4 ) )
             );

        row ( Y, i + 1 ) = row ( Y, i ) + h* (
                                 adams_moulton_koeff ( 0 ) * multMatrVec_omp ( F, k1 ) -
                                 adams_moulton_koeff ( 1 ) * multMatrVec_omp ( F, row ( Y, i ) ) +
                                 adams_moulton_koeff ( 2 ) * multMatrVec_omp ( F, row ( Y, i - 1 ) ) -
                                 adams_moulton_koeff ( 3 ) * multMatrVec_omp ( F, row ( Y, i - 2 ) ) +
                                 adams_moulton_koeff ( 4 ) * multMatrVec_omp ( F, row ( Y, i - 3 ) )
                             );
        }


    return Y;
    }

matr AdMltn_system_solve_matrix_mpi ( float h, int n, vec y_0, matr F )
    {
    vec k1, k2, k3, k4;
    matr Y ( n + 1, y_0.size() );

    row ( Y, 0 ) =  y_0;

    for ( int i = 0; i < 4; ++i )
        {
        k1 = multMatrVec_mpi ( F, row ( Y, i ) );
        k2 = multMatrVec_mpi ( F, row ( Y, i ) + h / 2. * k1 );
        k3 = multMatrVec_mpi ( F, row ( Y, i ) + h / 2. * k2 );
        k4 = multMatrVec_mpi ( F, row ( Y, i ) + h      * k3 );
        row ( Y, i + 1 ) = row ( Y, i ) + h / 6.* ( k1 + 2 * k2 + 2 * k3 + k4 );
        }

    for ( int i = 4; i < n; ++i )
        {
        k1 = row ( Y, i ) + h * (
                 adams_bashford_koeff ( 0 ) * multMatrVec_mpi ( F, row ( Y,  i ) ) -
                 adams_bashford_koeff ( 1 ) * multMatrVec_mpi ( F, row ( Y, i - 1 ) ) +
                 adams_bashford_koeff ( 2 ) * multMatrVec_mpi ( F, row ( Y, i - 2 ) ) -
                 adams_bashford_koeff ( 3 ) * multMatrVec_mpi ( F, row ( Y, i - 3 ) ) +
                 adams_bashford_koeff ( 4 ) * multMatrVec_mpi ( F, row ( Y, i - 4 ) )
             );

        row ( Y, i + 1 ) = row ( Y, i ) + h* (
                                 adams_moulton_koeff ( 0 ) * multMatrVec_mpi ( F, k1 ) -
                                 adams_moulton_koeff ( 1 ) * multMatrVec_mpi ( F, row ( Y, i ) ) +
                                 adams_moulton_koeff ( 2 ) * multMatrVec_mpi ( F, row ( Y, i - 1 ) ) -
                                 adams_moulton_koeff ( 3 ) * multMatrVec_mpi ( F, row ( Y, i - 2 ) ) +
                                 adams_moulton_koeff ( 4 ) * multMatrVec_mpi ( F, row ( Y, i - 3 ) )
                             );
        }


    return Y;
    }


vec eiler_function_solve ( float a, float b, float h, int n, float ( *f ) ( float, float ), float y0 )
    {

    float tp, xp;
    vec y ( n + 1 );

    y ( 0 ) = y0;


    for ( int i = 0; i < n; ++i )
        {
        tp = a + h * ( i ) + h/2;
        xp = y ( i ) + h/2 * f ( a + h * ( i ), y ( i ) );
        y ( i + 1 ) = y ( i ) + h * f ( tp, xp );
        }
    return y;
    }

matr eiler_system_solve_vector ( float a, float b, float h, int n, vec ( *vf ) ( float, vec ), vec vs )
    {
    matr Y ( n + 1, vs.size() );

    row ( Y, 0 ) = vs;

    for ( int i = 0; i < n; ++i )
        {
        row ( Y, i + 1 ) = row ( Y, i ) +
                             h * vf ( a + h * ( i ) + h/2,
                                      row ( Y, i ) + h/2 * vf ( a + h * ( i ),
                                              row ( Y, i ) )
                                    );
        }
    return Y;
    }

//TODO усовершенствованный метод Эйлера
matr eiler_system_solve_matrix ( float h, int n, vec y_0, matr F )
    {
    matr Y ( n + 1, y_0.size() );
    row ( Y, 0 ) = y_0;

    for ( int i = 0; i < n; ++i )
        row ( Y, i+1 ) = row ( Y, i ) + h * prod ( F, row ( Y, i ) );
    return Y;
    }
//TODO omp VM_mult
matr eiler_system_solve_matrix_omp ( float h, int n, vec y_0, matr F )
    {
    matr Y ( n + 1, y_0.size() );
    row ( Y, 0 ) = y_0;

    for ( int i = 0; i < n; ++i )
        row ( Y, i+1 ) = row ( Y, i ) + h * multMatrVec_omp ( F, row ( Y, i ) );
    return Y;
    }
//TODO mpi VM_mult
matr eiler_system_solve_matrix_mpi ( float h, int n, vec y_0, matr F )
    {
    matr Y ( n + 1, y_0.size() );
    row ( Y, 0 ) = y_0;

    for ( int i = 0; i < n; ++i )
        row ( Y, i+1 ) = row ( Y, i ) + h * multMatrVec_mpi ( F, row ( Y, i ) );
    return Y;
    }



vec rk_function_solve ( float a, float b, float h, int n, float ( *f ) ( float, float ), float y0 )
    {
    float k1, k2, k3, k4, k5, k6;
    vec y ( n+1 );

    y ( 0 ) = y0;
    for ( int i = 0; i < n; ++i )
        {
//         k1 = f ( a + h * i, y ( i ) );
//         k2 = f ( a + h * i + h / 2., y ( i ) + h/2 * k1 );
//         k3 = f ( a + h * i + h / 2., y ( i ) + h/2 * k2 );
//         k4 = f ( a + h * i + h, y ( i ) + h*k3 );
//         y ( i + 1 ) = y ( i ) + h / 6. * ( k1 + 2*k2 + 2*k3 + k4 );
        k1 = f ( a + h * i,          y ( i ) );
        k2 = f ( a + h * i + h / 4., y ( i ) + h / 4. * k1 );
        k3 = f ( a + h * i + h / 4., y ( i ) + h / 8. * k1 + h / 8. * k2 );
        k4 = f ( a + h * i + h / 2., y ( i ) - h / 2. * k2 + h *      k3 );
        k5 = f ( a + h * i + h * 3 / 4., y ( i ) + h * 3 / 16. * k1 + h * 9 / 16. * k4 );
        k6 = f ( a + h * i + h,      y ( i ) - h * 3 / 7. * k1 +
                 h * 2 / 7. * k2 + h * 12 / 7. * k3 - h * 12 / 7. * k4 + h * 8 / 7. * k5 );

        y ( i + 1 ) = y ( i ) + h / 90. * (
                          7 * k1 + 32 * k3 + 12 * k4 + 32 * k5 + 7 * k6
                      );
        }
    return y;
    }

matr rk_system_function_solve_vector ( float a, float b, float h, int n, vec ( *vf ) ( float, vec ), vec vs )
    {

    vec k1, k2, k3, k4, k5, k6;
    matr Y ( n+1, vs.size() );


    row ( Y, 0 ) = vs;

    for ( int i = 0; i < n; ++i )
        {
//         k1 = vf ( a + h * i, row ( Y, i ) );
//         k2 = vf ( a + h * i + h / 2., row ( Y, i ) + h / 2. * k1 );
//         k3 = vf ( a + h * ( i ) + h / 2., row ( Y, i ) + h / 2. * k2 );
//         k4 = vf ( a + h * ( i ) + h, row ( Y, i ) + h * k3 );
//         row ( Y, i ) = row ( Y, i ) + h / 6. * ( k1 + 2*k2 + 2*k3 + k4 );

        k1 = vf ( a + h * i,          row ( Y, i ) );
        k2 = vf ( a + h * i + h / 4., row ( Y, i ) + h / 4. * k1 );
        k3 = vf ( a + h * i + h / 4., row ( Y, i ) + h / 8. * k1 + h / 8. * k2 );
        k4 = vf ( a + h * i + h / 2., row ( Y, i ) - h / 2. * k2 + h      * k3 );
        k5 = vf ( a + h * i + h * 3 / 4., row ( Y, i ) + h * 3 / 16. * k1 + h * 9 / 16. * k4 );
        k6 = vf ( a + h * i + h,      row ( Y, i ) - h * 3 / 7. * k1 +
                  h * 2 / 7. * k2 + h * 12 / 7. * k3 - h * 12 / 7. * k4 + h * 8 / 7. * k5 );

        row ( Y, i + 1 ) = row ( Y, i ) + h / 90. * (
                                 7 * k1 + 32 * k3 + 12 * k4 + 32 * k5 + 7 * k6
                             );
        }
    return Y;
    }

//TODO протестировать
matr rk_system_solve_matrix ( float h, int n, vec y_0, matr F )
    {
    matr Y ( n + 1, y_0.size() );
    vec k1, k2, k3, k4;
    row ( Y, 0 ) = y_0;
    for ( int i = 0; i < n; ++i )
        {
        k1 = prod ( F, row ( Y, i ) );
        k2 = prod ( F, row ( Y, i ) + h / 2. * k1 );
        k3 = prod ( F, row ( Y, i ) + h / 2. * k2 );
        k4 = prod ( F, row ( Y, i ) + h      * k3 );
        row ( Y, i + 1 ) = row ( Y, i ) + h / 6.* ( k1 + 2 * k2 + 2 * k3 + k4 );
        }

    return Y;
    }

//FIXME добавить omp умножение
matr rk_system_solve_matrix_omp ( float h, int n, vec y_0, matr F )
    {
    vec k1, k2, k3, k4;
    matr Y ( n + 1, y_0.size() );
    row ( Y, 0 ) = y_0;


    for ( int i = 0; i < n; ++i )
        {
        k1 = multMatrVec_omp ( F, row ( Y, i ) );
        k2 = multMatrVec_omp ( F, row ( Y, i ) + h / 2. * k1 );
        k3 = multMatrVec_omp ( F, row ( Y, i ) + h / 2. * k2 );
        k4 = multMatrVec_omp ( F, row ( Y, i ) + h      * k3 );
        row ( Y, i + 1 ) = row ( Y, i ) + h / 6.* ( k1 + 2 * k2 + 2 * k3 + k4 );
        }

    return Y;
    }

//TODO add/test mpi VM_prod
matr rk_system_solve_matrix_mpi ( float h, int n, vec y_0, matr F )
    {
    vec k1, k2, k3, k4;
    matr Y ( n + 1, y_0.size() );
    row ( Y, 0 ) = y_0;


    for ( int i = 0; i < n; ++i )
        {
        k1 = multMatrVec_mpi ( F, row ( Y, i ) );
        k2 = multMatrVec_mpi ( F, row ( Y, i ) + h / 2. * k1 );
        k3 = multMatrVec_mpi ( F, row ( Y, i ) + h / 2. * k2 );
        k4 = multMatrVec_mpi ( F, row ( Y, i ) + h      * k3 );
        row ( Y, i + 1 ) = row ( Y, i ) + h / 6.* ( k1 + 2 * k2 + 2 * k3 + k4 );
        }

    return Y;
    }

}
