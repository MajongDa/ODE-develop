#ifndef SIMFOR_MATRIX_VECTOR_HPP
#define SIMFOR_MATRIX_VECTOR_HPP

#include "vector.hpp"
#include "matrix.hpp"
#include <ompi/mpi/cxx/mpicxx.h>
#include <omp.h>


namespace simfor
{

template<class C1, class C2>
vector<C1> prod ( const matrix<C1> &R, const vector<C2> &l )
    {
    int n = R.size2() == l.size() ? l.size() : 0;
    vector<C1> v ( n );
    for ( unsigned i = 0; i < n; i++ )
        {
        v ( i ) = 0;
        for ( unsigned j = 0; j < n; j++ )
            v ( i ) += R ( i, j ) * l ( j );
        }
    return v;
    }

// template<class E1, class E2>
// vector<float> multMatrVec_omp ( const E1 &A, const E2 &b )
//     {
//     int n = A.size1() == b.size() ? b.size() : 0;
//     vector<float> v ( n, 0 );
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
vector<E1> prod_omp ( const matrix<E1> &A, const vector<E2> &b )
    {
    int n = A.size2() == b.size() ? b.size() : 0;
    vector<E1> v ( n );

    #pragma omp parallel for shared( A, b, v) schedule(auto)
    for ( int i = 0; i < n; i++ )
        {
        for ( int j = 0; j < n; j++ )
            v ( i ) += A ( i, j ) * b ( j );

        }
    return v;
    }


template<class E1, class E2>
void prod_mpi ( const matrix<E1> &A, const vector<E2> &b, vector<float> &c )
    {

    int n = A.size2() == b.size() ? b.size() : 0;


    MPI::Intracomm world = MPI::COMM_WORLD;
    MPI::Status status;
    MPI::Init();
    int p = world.Get_size();
    int r = world.Get_rank();
    int cnt = n / p;
    int from = r * cnt;
    int to = n;

    if ( r != p-1 ) to = from + cnt;

    for ( unsigned i = from; i < to; i++ )
        for ( unsigned j = 0; j < n; j++ )
            c ( i ) += A ( i, j ) * b ( j );


    if ( r != 0 )
        {
        world.Send ( &from, 1, MPI::INT, 0, 1);
        world.Send ( &to, 1, MPI::INT, 0, 2 );
        world.Send ( &c ( from ), ( to - from ), MPI::FLOAT, 0, 3);
        }
    else
        {
        for ( int i = 1; i < p; ++i )
            {
            world.Recv ( &from, 1, MPI::INT, i, 1, status );
            world.Recv ( &to, 1, MPI::INT, i, 2, status );
            world.Recv ( &c ( from ), ( to - from ), MPI::FLOAT, i, 3, status );
            }
        }
    MPI::Finalize();
    }

template<class E1, class E2>
vector<float> prod_mpi ( const matrix<E1> &A, const vector<E2> &b )
    {

    int n = A.size2() == b.size() ? b.size() : 0;


    MPI::Intracomm world = MPI::COMM_WORLD;
    MPI::Status status;
    //MPI::Init();
    int p = world.Get_size();
    int r = world.Get_rank();
    int cnt = n / p;
    int from = r * cnt;
    int to = n;
    vector<float> c ( n );

    if ( r != p-1 ) to = from + cnt;

    for ( unsigned i = from; i < to; i++ )
        for ( unsigned j = 0; j < n; j++ )
            c ( i ) += A ( i, j ) * b ( j );


    if ( r != 0 )
        {
        world.Send ( &from, 1, MPI::INT, 0, 1);
        world.Send ( &to, 1, MPI::INT, 0, 2 );
        world.Send ( &c ( from ), ( to - from ), MPI::FLOAT, 0, 3);
        }
    else
        {
        for ( int i = 1; i < p; ++i )
            {
            world.Recv ( &from, 1, MPI::INT, i, 1, status );
            world.Recv ( &to, 1, MPI::INT, i, 2, status );
            world.Recv ( &c ( from ), ( to - from ), MPI::FLOAT, i, 3, status );
            }
        }
    //MPI::Finalize();
    return c;
    }
}
#endif





