#ifndef SIMFOR_MATMULT_HPP
#define SIMFOR_MATMULT_HPP

#include "internal/types.hpp"
#include <boost/mpi.hpp>
#include <omp.h>



namespace simfor{
    template<class E1, class E2>
    vec multMatrVec ( const E1 &A, const E2 &b )
    {
    int n = A.size1() == b.size() ? b.size() : 0;
    vec v ( n );
    for ( unsigned i = 0; i < n; i++ )
        {
        v ( i ) = 0;
        for ( unsigned j = 0; j < n; j++ )
            v ( i ) += A ( i, j ) * b ( j );
        }
    return v;
    }

template<class E1, class E2>
vec multMatrVec_omp ( const E1 &A, const E2 &b )
    {
    int n = A.size1() == b.size() ? b.size() : 0;
    vec v ( n );
    int i, j;

    #pragma omp parallel for shared( A, b, v) private(i, j) schedule(static)
    for ( i = 0; i < n; i++ )
        {
        float sum = 0;
        #pragma omp simd reduction(+:sum)
        for ( j = 0; j < b.size(); j++ )
            {
                #pragma omp atomic
            sum += A ( i, j ) * b ( j );
            }
        #pragma omp atomic
        v ( i ) += sum;
        }
    return v;
    }


template<class E1, class E2>
void multMatrVec_mpi ( const E1 &A, const E2 &b, vec &c )
    {
    namespace mpi = boost::mpi;
    int n = A.size1() == b.size() ? b.size() : 0;
    mpi::communicator world;

    int p = world.size();
    int r = world.rank();
    int cnt = n / p;
    int from = r * cnt;
    int to = n;

    if ( r != p-1 ) to = from + cnt;

    for ( unsigned i = from; i < to; i++ )
        {
        c ( i ) = 0;
        for ( unsigned j = 0; j < n; j++ )
            c ( i ) += A ( i, j ) * b ( j );

        }

    world.barrier ();

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
    }

    template<class E1, class E2>
vec multMatrVec_mpi ( const E1 &A, const E2 &b )
    {
    namespace mpi = boost::mpi;
    int n = A.size1() == b.size() ? b.size() : 0;
    mpi::communicator world;

    int p = world.size();
    int r = world.rank();
    int cnt = n / p;
    int from = r * cnt;
    int to = n;
    vec c ( n );

    if ( r != p-1 ) to = from + cnt;

    for ( unsigned i = from; i < to; i++ )
        {
        c ( i ) = 0;
        for ( unsigned j = 0; j < n; j++ )
            c ( i ) += A ( i, j ) * b ( j );

        }

    world.barrier ();

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
            return c;
        }
    }
}

#endif
