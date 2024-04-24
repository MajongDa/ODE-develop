#include <simfor/elementary.hpp>

namespace simfor
{
float v_norm2 ( vec a, vec b )
    {
    float tb = 0, ta = 0;

    for ( int i = 0; i < a.size (); i++ )
        {
        ta += a[i]*a[i];
        tb += b[i]*b[i];
        }
    return sqrt ( tb-ta );
    }

float v_norm2_omp ( vec a, vec b )
    {
    float tb = 0, ta = 0;
    #pragma omp parallel for private(tb, ta)
    for ( int i = 0; i < a.size (); i++ )
        {
        ta += a[i]*a[i];
        tb += b[i]*b[i];
        }
    return sqrt ( tb-ta );
    }

float m_norm4 ( mat A )
    {
    float res;
    for ( int i = 0; i < A.size1(); i++ )
        for ( int j = 0; j < A.size2(); j++ )
            res += A ( j, i ) *A ( j, i );
    return sqrt ( res );
    }

float m_norm4_omp ( mat A )
    {
    float res;
    #pragma omp parallel for reduction(+:res)
    for ( int i = 0; i < A.size1(); i++ )
        for ( int j = 0; j < A.size2(); j++ )
            res += A ( j, i ) *A ( j, i );
    return sqrt ( res );
    }

float psum ( vec a )
    {
    float res;
    for ( int i = 0; i < a.size(); i++ )
        res += a ( i );
    return res;
    }

float psum_omp ( vec a )
    {
    float res;
    double time;
    time = omp_get_wtime();
    #pragma omp parallel for reduction(+:res)
    for ( int i = 0; i < a.size(); i++ )
        res += a ( i );
    time = time - omp_get_wtime();
    return res;
    }
}
