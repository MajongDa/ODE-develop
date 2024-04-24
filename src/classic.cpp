#include <chrono>
#include <simfor/elementary.hpp>
#include <boost/operators.hpp>

/*test time globs
std::clock_t T1, T2;
double RT1, RT2;
T1 = std::clock();
T2 = std::clock();*/

namespace simfor
{
float scalar_mult ( vec a, vec b )
    {
    float res = 0;
    for ( int i = 0; i < a.size(); ++i )
        res += a ( i ) * b ( i );
    return res;
    }

float scalar_mult_omp ( vec a, vec b )
    {
    unsigned i;
    float res = 0;
    #pragma omp parallel for reduction(+ : res)
    for ( i = 0; i < a.size(); ++i )
        res += a[i] * b[i];
    return res;
    }

float scalar_mult_mpi ( vec a, vec b )
    {
    boost::mpi::communicator world;

    int p = world.size();
    int r = world.rank();
    int n = a.size();
    float s;

    int cnt = n / p;
    int from = r * cnt;
    int to = n;


    //boost::mpi::timer t;
    //t.restart();

    if ( r != p-1 ) to = from + cnt;

    for ( int i = from; i < to; ++i )
        s += a ( i ) * b ( i );
    world.barrier ();

    if ( !r ) reduce ( world, s, std::plus<float>(), 0 );

    return s;

    }
}
