
#include "simfor/multMatrVec.hpp"
#include <boost/numeric/ublas/io.hpp>
#include <boost/mpi.hpp>
#include <iostream>

using namespace std;

int main ( int argc, char **argv )
    {
    namespace mpi = boost::mpi;
    int n = atoi ( argv[1] );
    int p = atoi ( argv[2] );
    simfor::vec vector ( n ), res;
    simfor::matr matrix ( n, n );
    for ( unsigned i = 0; i < n; i++ )
        for ( unsigned j = 0; j < n; j++ )
            {
            matrix ( i, j ) = 1;
            vector ( i ) = 2;
            }
    double start, end;
    start = omp_get_wtime();
    omp_set_num_threads ( p );
    res = simfor::multMatrVec_omp ( matrix, vector );
    end = omp_get_wtime();
    std::cout << end - start << "\t" << "\n";
    }
