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
//     t = clock();
//     res = simfor::multMatrVec ( matrix, vector );
//     t = ( clock() - t ) / CLOCKS_PER_SEC ;
//     std::cout << t << "\t" << "\n";
//     t = clock();
//     omp_set_num_threads(4);
//     res = simfor::multMatrVec_omp ( matrix, vector );
//     t = ( clock() - t ) / CLOCKS_PER_SEC ;
//     std::cout << t << "\t" << "\n";
    mpi::environment env;
    mpi::communicator word;
    res = simfor::multMatrVec_mpi ( matrix, vector );

    if ( !word.rank() )
        {

        t = ( clock() - t ) / CLOCKS_PER_SEC ;
        std::cout << t << "\t" << "\n";
        }
    }
