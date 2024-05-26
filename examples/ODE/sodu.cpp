
#include "sodu.hpp"

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
matrix<float> adams5_system_solve_matrix ( float h, int n, vector<float> y_0, matrix<float> F )
    {
    vector<float> k1, k2, k3, k4, k5;
    matrix<float> Y ( n + 1, y_0.size() );

    Y.set_submatrix ( rk_system_solve_matrix ( h, 4, y_0, F ), 0, 5, 0, y_0.size() );

    for ( int i = 4; i < n; ++i )
        {
        k1 = adams_bashford_koeff ( 0 ) * prod ( F, row ( Y,  i ) );
        k2 = adams_bashford_koeff ( 1 ) * prod ( F, row ( Y, i - 1 ) );
        k3 = adams_bashford_koeff ( 2 ) * prod ( F, row ( Y, i - 2 ) );
        k4 = adams_bashford_koeff ( 3 ) * prod ( F, row ( Y, i - 3 ) );
        k5 = adams_bashford_koeff ( 4 ) * prod ( F, row ( Y, i - 4 ) );
        Y.set_row ( Y.get_row ( i ) + h * ( k1 - k2 + k3 - k4 + k5), i+1 );
        }


    return Y;
    }

//TODO add omp for scalar-vector mult
matrix<float> adams5_system_solve_matrix_omp ( float h, int n, vector<float> y_0, matrix<float> F )
    {
    vector<float> k1, k2, k3, k4, k5;
    matrix<float> Y ( n + 1, y_0.size() );

    Y.set_submatrix ( rk_system_solve_matrix_omp ( h, 4, y_0, F ), 0, 5, 0, y_0.size() );


    for ( int i = 4; i < n; ++i )
        {
        k1 = adams_bashford_koeff ( 0 ) * prod_omp ( F, row ( Y,  i ) );
        k2 = adams_bashford_koeff ( 1 ) * prod_omp ( F, row ( Y, i - 1 ) );
        k3 = adams_bashford_koeff ( 2 ) * prod_omp ( F, row ( Y, i - 2 ) );
        k4 = adams_bashford_koeff ( 3 ) * prod_omp ( F, row ( Y, i - 3 ) );
        k5 = adams_bashford_koeff ( 4 ) * prod_omp ( F, row ( Y, i - 4 ) );
        Y.set_row ( Y.get_row ( i ) + h * ( k1 - k2 + k3 - k4 + k5), i+1 );
        }
    return Y;
    }

//TODO add omp for scalar-vector mult
matrix<float> adams5_system_solve_matrix_mpi ( float h, int n, vector<float> y_0, matrix<float> F )
    {
    vector<float> k1, k2, k3, k4, k5;
    matrix<float> Y ( n + 1, y_0.size() );


    Y.set_submatrix ( rk_system_solve_matrix_mpi ( h, 4, y_0, F ), 0, 5, 0, y_0.size() );

    for ( int i = 4; i < n; ++i )
        {
        k1 = adams_bashford_koeff ( 0 ) * prod_mpi ( F, row ( Y,  i ) );
        k2 = adams_bashford_koeff ( 1 ) * prod_mpi ( F, row ( Y, i - 1 ) );
        k3 = adams_bashford_koeff ( 2 ) * prod_mpi ( F, row ( Y, i - 2 ) );
        k4 = adams_bashford_koeff ( 3 ) * prod_mpi ( F, row ( Y, i - 3 ) );
        k5 = adams_bashford_koeff ( 4 ) * prod_mpi ( F, row ( Y, i - 4 ) );
        Y.set_row ( Y.get_row ( i ) + h * ( k1 - k2 + k3 - k4 + k5), i+1 );
        }


    return Y;
    }

//Метод Адамса-Моултона
matrix<float> AdMltn_system_solve_matrix ( float h, int n, vector<float> y_0, matrix<float> F )
    {
    vector<float> k1, k2, k3, k4;
    matrix<float> Y ( n + 1, y_0.size() );

    Y.set_submatrix ( rk_system_solve_matrix ( h, 4, y_0, F ), 0, 5, 0, y_0.size() );

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

matrix<float> AdMltn_system_solve_matrix_omp ( float h, int n, vector<float> y_0, matrix<float> F )
    {
    vector<float> k1, k2, k3, k4;
    matrix<float> Y ( n + 1, y_0.size() );

    Y.set_submatrix ( rk_system_solve_matrix_omp ( h, 4, y_0, F ), 0, 5, 0, y_0.size() );

    for ( int i = 4; i < n; ++i )
        {
        k1 = row ( Y, i ) + h * (
                 adams_bashford_koeff ( 0 ) * prod_omp ( F, row ( Y,  i ) ) -
                 adams_bashford_koeff ( 1 ) * prod_omp ( F, row ( Y, i - 1 ) ) +
                 adams_bashford_koeff ( 2 ) * prod_omp ( F, row ( Y, i - 2 ) ) -
                 adams_bashford_koeff ( 3 ) * prod_omp ( F, row ( Y, i - 3 ) ) +
                 adams_bashford_koeff ( 4 ) * prod_omp ( F, row ( Y, i - 4 ) )
             );

        row ( Y, i + 1 ) = row ( Y, i ) + h* (
                               adams_moulton_koeff ( 0 ) * prod_omp ( F, k1 ) -
                               adams_moulton_koeff ( 1 ) * prod_omp ( F, row ( Y, i ) ) +
                               adams_moulton_koeff ( 2 ) * prod_omp ( F, row ( Y, i - 1 ) ) -
                               adams_moulton_koeff ( 3 ) * prod_omp ( F, row ( Y, i - 2 ) ) +
                               adams_moulton_koeff ( 4 ) * prod_omp ( F, row ( Y, i - 3 ) )
                           );
        }


    return Y;
    }

matrix<float> AdMltn_system_solve_matrix_mpi ( float h, int n, vector<float> y_0, matrix<float> F )
    {
    vector<float> k1, k2, k3, k4;
    matrix<float> Y ( n + 1, y_0.size() );

    Y.set_submatrix ( rk_system_solve_matrix_mpi ( h, 4, y_0, F ), 0, 5, 0, y_0.size() );

    for ( int i = 4; i < n; ++i )
        {
        k1 = row ( Y, i ) + h * (
                 adams_bashford_koeff ( 0 ) * prod_mpi ( F, row ( Y,  i ) ) -
                 adams_bashford_koeff ( 1 ) * prod_mpi ( F, row ( Y, i - 1 ) ) +
                 adams_bashford_koeff ( 2 ) * prod_mpi ( F, row ( Y, i - 2 ) ) -
                 adams_bashford_koeff ( 3 ) * prod_mpi ( F, row ( Y, i - 3 ) ) +
                 adams_bashford_koeff ( 4 ) * prod_mpi ( F, row ( Y, i - 4 ) )
             );

        row ( Y, i + 1 ) = row ( Y, i ) + h* (
                               adams_moulton_koeff ( 0 ) * prod_mpi ( F, k1 ) -
                               adams_moulton_koeff ( 1 ) * prod_mpi ( F, row ( Y, i ) ) +
                               adams_moulton_koeff ( 2 ) * prod_mpi ( F, row ( Y, i - 1 ) ) -
                               adams_moulton_koeff ( 3 ) * prod_mpi ( F, row ( Y, i - 2 ) ) +
                               adams_moulton_koeff ( 4 ) * prod_mpi ( F, row ( Y, i - 3 ) )
                           );
        }


    return Y;
    }


matrix<float> eiler_system_solve_matrix ( float h, int n, vector<float> y_0, matrix<float> F )
    {
    matrix<float> Y ( n + 1, y_0.size() );
    vector<float> a ( y_0.size() );
    Y.set_row ( y_0, 0 );
    for ( int i = 0; i < n; ++i )
        {
        Y.set_row ( Y.get_row ( i ) + h * prod ( F, Y.get_row ( i ) ), i+1 );
        }
    return Y;
    }

matrix<float> eiler_system_solve_matrix_omp ( float h, int n, vector<float> y_0, matrix<float> F )
    {
    matrix<float> Y ( n + 1, y_0.size() );
    row ( Y, 0 ) = y_0;

    for ( int i = 0; i < n; ++i )
        Y.set_row ( Y.get_row ( i ) + h * prod_omp ( F, Y.get_row ( i ) ), i+1 );
    return Y;
    }

matrix<float> eiler_system_solve_matrix_mpi ( float h, int n, vector<float> y_0, matrix<float> F )
    {
    matrix<float> Y ( n + 1, y_0.size() );
    row ( Y, 0 ) = y_0;

    for ( int i = 0; i < n; ++i )
        Y.set_row ( Y.get_row ( i ) + h * prod_mpi ( F, Y.get_row ( i ) ), i+1 );
    return Y;
    }

matrix<float> rk_system_solve_matrix ( float h, int n, vector<float> y_0, matrix<float> F )
    {
    matrix<float> Y ( n + 1, y_0.size() );
    vector<float> k1, k2, k3, k4;
    Y.set_row ( y_0, 0 );
    for ( int i = 0; i < n; ++i )
        {
        k1 = prod ( F, Y.get_row ( i ) );
        k2 = prod ( F, Y.get_row ( i ) + h / 2. * k1 );
        k3 = prod ( F, Y.get_row ( i ) + h / 2. * k2 );
        k4 = prod ( F, Y.get_row ( i ) + h      * k3 );
        Y.set_row ( ( Y.get_row ( i ) + h / 6.* ( k1 + 2 * k2 + 2 * k3 + k4 ) ), i+1 );
        }
    return Y;
    }

matrix<float> rk_system_solve_matrix_omp ( float h, int n, vector<float> y_0, matrix<float> F )
    {
    vector<float> k1, k2, k3, k4;
    matrix<float> Y ( n + 1, y_0.size() );
    row ( Y, 0 ) = y_0;

    for ( int i = 0; i < n; ++i )
        {
        k1 = prod_omp ( F, row ( Y, i ) );
        k2 = prod_omp ( F, row ( Y, i ) + h / 2. * k1 );
        k3 = prod_omp ( F, row ( Y, i ) + h / 2. * k2 );
        k4 = prod_omp ( F, row ( Y, i ) + h      * k3 );
        Y.set_row ( row ( Y, i ) + h / 6.* ( k1 + 2 * k2 + 2 * k3 + k4 ), i+1);
        }

    return Y;
    }

matrix<float> rk_system_solve_matrix_mpi ( float h, int n, vector<float> y_0, matrix<float> F )
    {
    vector<float> k1, k2, k3, k4;
    matrix<float> Y ( n + 1, y_0.size() );
    row ( Y, 0 ) = y_0;


    for ( int i = 0; i < n; ++i )
        {
        k1 = prod_mpi ( F, row ( Y, i ) );
        k2 = prod_mpi ( F, row ( Y, i ) + h / 2. * k1 );
        k3 = prod_mpi ( F, row ( Y, i ) + h / 2. * k2 );
        k4 = prod_mpi ( F, row ( Y, i ) + h      * k3 );
        Y.set_row ( row ( Y, i ) + h / 6.* ( k1 + 2 * k2 + 2 * k3 + k4 ), i+1 );
        }

    return Y;
    }
}

