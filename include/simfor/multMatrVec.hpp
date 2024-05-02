#ifndef SIMFOR_MATMULT_HPP
#define SIMFOR_MATMULT_HPP

#include "internal/types.hpp"
#include <boost/mpi.hpp>
#include <boost/mpi/collectives.hpp>
#include <omp.h>
#include <iostream>
#include <boost/numeric/ublas/io.hpp>



namespace simfor
{

template<class E1, class E2>
vec multMatrVec ( const E1 &A, const E2 &b )
    {
    int n = A.size2() == b.size() ? b.size() : 0;
    vec v ( n );
    float s = 0;
    for ( unsigned i = 0; i < n; i++ )
        {
        v ( i ) = 0;
        for ( unsigned j = 0; j < n; j++ )
            v ( i ) += A ( i, j ) * b ( j );
        }
    return v;
    }

// template<class E1, class E2>
// vec multMatrVec_omp ( const E1 &A, const E2 &b )
//     {
//     int n = A.size1() == b.size() ? b.size() : 0;
//     vec v ( n, 0 );
//     #pragma omp parallel for collapse(2) schedule(auto)
//     for ( int i = 0; i < n; i++ )
//         {
//         //#pragma omp simd reduction(+:sum)
//         for ( int j = 0; j < b.size(); j++ )
//             v ( i ) += A ( i, j ) * b ( j );
//         }
//     return v;
//     }

template<class E1, class E2>
vec multMatrVec_omp ( const E1 &A, const E2 &b )
    {
    int n = A.size2() == b.size() ? b.size() : 0;
    vec v ( n );

    #pragma omp parallel for shared( A, b, v) schedule(auto)
    for ( int i = 0; i < n; i++ )
        {
        for ( int j = 0; j < b.size(); j++ )
            v ( i ) += A ( i, j ) * b ( j );

        }
    return v;
    }


template<class E1, class E2>
void multMatrVec_mpi ( const E1 &A, const E2 &b, vec &c )
    {
    namespace mpi = boost::mpi;
    int n = A.size2() == b.size() ? b.size() : 0;

    mpi::environment env;
    mpi::communicator world;
    int p = world.size();
    int r = world.rank();
    int cnt = n / p;
    int from = r * cnt;
    int to = n;

    if ( r != p-1 ) to = from + cnt;

    for ( unsigned i = from; i < to; i++ )
        {
        for ( unsigned j = 0; j < n; j++ )
            c ( i ) += A ( i, j ) * b ( j );

        }

    if ( r != 0 )
        {
        mpi::request reqs[3];
        reqs[0] = world.isend ( 0, 1, from );
        reqs[1] = world.isend ( 0, 2, to );
        reqs[2] = world.isend ( 0, 3, &c ( from ), ( to - from ) );
        mpi::wait_all ( reqs, reqs+3 );
        }
    else
        {
        for ( int i = 1; i < p; ++i )
            {
            mpi::request reqs[3];
            reqs[0] = world.irecv ( i, 1, from );
            reqs[1] = world.irecv ( i, 2, to );
            reqs[2] = world.irecv ( i, 3, &c ( from ), ( to - from ) );
            mpi::wait_all ( reqs, reqs+3 );
            }
        }
    }

template<class E1, class E2>
vec multMatrVec_mpi ( const E1 &A, const E2 &b )
    {
    namespace mpi = boost::mpi;
    int n = A.size2() == b.size() ? b.size() : 0;


    mpi::communicator world;

    int p = world.size();
    int r = world.rank();
    int cnt = n / p;
    int from = r * cnt;
    int to = n;
    vec c ( n, 0 );

    if ( r != p-1 ) to = from + cnt;

    for ( unsigned i = from; i < to; i++ )
        {
        for ( unsigned j = 0; j < n; j++ )
            c ( i ) += A ( i, j ) * b ( j );
        }

    if ( r != 0 )
        {
        world.send ( 0, 1, from );
        world.send ( 0, 2, to );
        world.send ( 0, 3, &c ( from ), ( to - from ) );
        }
    else
        {
        for ( int i = 1; i < p; ++i )
            {
            world.recv ( i, 1, from );
            world.recv ( i, 2, to );
            world.recv ( i, 3, &c ( from ), ( to - from ) );
            }
        }
    return c;
    }
}

#endif



