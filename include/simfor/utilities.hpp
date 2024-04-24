#ifndef SIMFOR_UTILITIES_HPP
#define SIMFOR_UTILITIES_HPP

#include "internal/types.hpp"
namespace simfor
{
float koef_matrix_settings ( unsigned j, unsigned i, unsigned n )
    {
    float value;
    if ( j == i ) value = -4.f;
    else if ( j > 0 && ( ( j-1 ) == i ) ) value = 1.f;
    else if ( j < ( n-1 ) && ( ( j+1 ) == i ) ) value = 1.f;
    else value = 0.0f;
    return value;
    }

matr odu_matrix_create ( unsigned n )
    {
    matr matr ( n, n );
    for ( int i = 0; i < n; i++ )
        for ( int j = 0; j < n; j++ )
            matr ( i, j ) = koef_matrix_settings ( j, i, n );
    return matr;
    }
}
#endif
