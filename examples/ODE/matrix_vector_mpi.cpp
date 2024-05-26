
#include "multMatrixVector.hpp"
#include <iostream>

using namespace std;

int main ( int argc, char **argv )
    {
    int n = atoi ( argv[1] );
    simfor::vector<float> vector ( n ), res ( n );
    simfor::matrix<float> matrix ( n, n );
    for ( unsigned i = 0; i < n; i++ )
        for ( unsigned j = 0; j < n; j++ )
            {
            matrix ( i, j ) = 1;
            vector ( i ) = 2;
            }
    double start, end;
    start = omp_get_wtime();
    res = simfor::prod_mpi ( matrix, vector );
    end = omp_get_wtime();
    std::cout << end - start << "\n";
    return 0;
    }
