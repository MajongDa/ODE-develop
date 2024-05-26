#include <simfor/odu.hpp>
#include <simfor/internal/types.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cmath>
#include <iostream>
#include <fstream>
#include "matplot/matplot.h"

double odu ( double t, double x )
    {
    return ( 6*x - 13*t*t*t - 22*t*t + 17*t - 11 + sin ( t ) );
    //return t*t - 2*x;
    }
double solution ( double t )
    {
    return ( 13.*t*t*t ) /6 + ( 19.*t*t ) /4 - ( 5.*t ) /4 + 13./8 - cos ( t ) /37. - ( 6*sin ( t ) ) /37. + ( 119.*exp ( 6*t ) ) /296;
    //return ( 3./4 * exp ( -2*t ) + 1./2*t*t - 1./2*t + 1./4 );
    }

simfor::vec absvec ( simfor::vec v )
    {
    simfor::vec vabs ( v.size() );
    for ( unsigned i = 0; i < v.size(); ++i )
        vabs ( i ) = std::abs ( v ( i ) );
    return vabs;
    }

template <class T>
std::vector<T> convert( boost::numeric::ublas::vector<T> a )
    {
    std::vector<T> b(a.size());
    std::copy(a.begin(), a.end(), b.begin());
    return b;
    }


int main()
    {

    double a = 0, b = 5, x0 = solution ( a );
    double h = 1e-4;
    int n = std::ceil ( ( b-a)/h );
    simfor::matr adaptive;
    simfor::matr data ( 4, n+1 );

    for ( int i = 0; i < n+1; ++i )
        {
        data ( 0, i ) = ( a + ( b-a ) / n * i );
        data ( 1, i ) = solution ( data ( 0, i ) );
        }

    using namespace matplot;
    hold(on);

    //Euler method
    row ( data, 2 ) = simfor::eiler_function_solve ( a, b, h, n, odu, x0 );
    row ( data, 3 ) = absvec ( row ( data, 2 ) - row ( data, 1 ) );

    std::cout << "Euler2 abs error " << std::abs (row ( data, 2 ) ( n ) - solution ( b ) ) << "\n";


    row ( data, 2 ) = simfor::rk4_function_solve ( a, b, h, n, odu, x0 );
    row ( data, 3 ) = absvec ( row ( data, 2 ) - row ( data, 1 ) );

    std::cout << "RK4 abs error "
    << std::abs(row ( data, 2 ) ( n ) - solution ( b )) << "\n";

    row ( data, 2 ) = simfor::rk6_function_solve ( a, b, h, n, odu, x0 );
    row ( data, 3 ) = absvec ( row ( data, 2 ) - row ( data, 1 ) );
    auto rk6 = plot ( convert<float> ( row ( data, 0 ) ), convert<float> ( row ( data, 2 ) ) );

    std::cout << "RK6 abs error "
    << std::abs(row ( data, 2 ) ( n ) - solution ( b )) << "\n";

    adaptive = simfor::rk45 ( a, b, odu, x0, 1e-7, 1e-8, 1e-3 );
    auto rkf45 = plot ( convert<float>( row(adaptive, 0 ) ), convert<float>( row(adaptive, 1 ) ), "+" );
    std::cout << "RKF4(5) abs error " <<
    std::abs ( solution ( b ) - adaptive ( 1, adaptive.size2()-1 ) ) << "\n";


    // adaptive = simfor::rk54 ( a, b, odu, x0, 1e-7, 1e-2, 1e-1 );
    // auto dopr54 = plot ( convert<float>( row(adaptive, 0 ) ), convert<float>( row(adaptive, 1 ) ), "o" );
    // std::cout << "RKDoPr5(4) abs error " <<
    // std::abs ( solution ( b ) - adaptive ( 1, adaptive.size2()-1 ) ) << "\n";


    row ( data, 2 ) = simfor::AdMiln_function_solve ( a, b, h, n, odu, x0 );
    row ( data, 3 ) = absvec ( row ( data, 2 ) - row ( data, 1 ) );
    std::cout << "AdMiln3 abs error " <<
    row ( data, 2 ) ( n ) - solution ( b ) << "\n";

    row ( data, 2 ) = simfor::adams5_function_solve ( a, b, h, n, odu, x0 );
    row ( data, 3 ) = absvec ( row ( data, 2 ) - row ( data, 1 ) );
    std::cout << "AdBsh5 abs error " <<
    std::abs ( row ( data, 2 ) ( n ) - solution ( b ) ) << "\n";

    row ( data, 2 ) = simfor::AdMltn_function_solve ( a, b, h, n, odu, x0 );
    row ( data, 3 ) = absvec ( row ( data, 2 ) - row ( data, 1 ) );
    std::cout << "AdMlt5 abs error " <<
    std::abs ( row ( data, 2 ) ( n ) - solution ( b ) ) << "\n";

    ::matplot::legend ( { rk6, rkf45/*, dopr54*/ }, { "RK6","RK_{}F45", "RK_{}DoPr54"} );
    hold(off);
    show();
    return 0;
    }
