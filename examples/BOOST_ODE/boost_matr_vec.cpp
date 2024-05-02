#include "simfor/multMatrVec.hpp"
#include <boost/numeric/ublas/io.hpp>
#include <boost/mpi.hpp>
#include <iostream>

using namespace std;

int main ( int argc, char **argv )
    {
    namespace mpi = boost::mpi;
    int n = atoi( argv[1]);
    simfor::vec vector ( n ), res;
    simfor::matr matrix ( n, n );
    for ( unsigned i = 0; i < n; i++ )
        for ( unsigned j = 0; j < n; j++ )
            {
            matrix ( i, j ) = 1;
            vector ( i ) = 2;
            }
    double t;
    t = clock();
    res = prod( matrix, vector );
    t = ( clock() - t ) / CLOCKS_PER_SEC ;
    std::cout << t << "\t" << "\n";

    t = clock();
    res = simfor::multMatrVec ( matrix, vector );
    t = ( clock() - t ) / CLOCKS_PER_SEC ;
    std::cout << t << "\t" << "\n";
    }

