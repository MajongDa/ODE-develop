#include "multMatrixVector.hpp"
#include <iostream>

int main ( int argc, char **argv )
    {
    unsigned n = atoi ( argv[1] );
    simfor::vector<float> v ( n ), res;
    simfor::matrix<float> M ( n, n );
    for ( unsigned i = 0; i < n; i++ )
        for ( unsigned j = 0; j < n; j++ )
            {
            M ( i, j ) = 1;
            v ( i ) = 2;
            }
    double t, min_t = 1000, avg_t=0;

    //burnin
    for ( int i = 0; i < 10; ++i )
        res = simfor::prod ( M, v );

    for ( int i = 0; i < 10; ++i )
        {
        t = clock();
        res = simfor::prod ( M, v );
        avg_t += t = ( clock() - t ) / CLOCKS_PER_SEC ;
        min_t = min_t > t ? t : min_t;
        }
    std::cout << "BOOSTless Matrix x Vector " << avg_t/10  << " min_time "<< min_t << "\t" << "\n";
    return 0;
    }


