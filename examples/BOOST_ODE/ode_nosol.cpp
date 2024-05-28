#include <simforode/odu.hpp>
#include <simforode/internal/types.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cmath>
#include <iostream>
#include <fstream>
#include "matplot/matplot.h"

double odu ( double t, double x )
    {
    return -x*x - x*x*x + std::cos(5*t);
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

    double a = 0, b = 5, x0 = 3;
    double h = 1e-1;
    int n = std::ceil ( ( b-a)/h )+1;
    simfor::matr adaptive;
    simfor::vec t ( n ), u;

    for ( int i = 0; i < n; ++i )
        {
        t ( i ) = ( a + h * i );
        }

    using namespace matplot;



    hold(on);

    //Euler method
    u = simfor::eiler_function_solve ( a, b, h, n, odu, x0 );
    auto eiler = plot ( convert<float> ( t ), convert<float> (  u ), "*" );


    u = simfor::rk4_function_solve ( a, b, h, n, odu, x0 );
     std::cout << t << "\n";
     std::cout << u << "\n";




    u = simfor::rk5_function_solve ( a, b, h, n, odu, x0 );

    auto rk5 = plot ( convert<float> ( t ), convert<float> ( u ) , "s" );
hold(off);

    adaptive = simfor::rk45 ( a, b, odu, x0, 1e-4, 1e-6, 1e-1 );
    //auto rkf45 = plot ( convert<float>( row(adaptive, 0 ) ), convert<float>( row(adaptive, 1 ) ), "+" );
     std::cout << row ( adaptive, 1 ) << "\n";

    adaptive = simfor::rk54 ( a, b, odu, x0, 1e-3, 1e-6, 1e-1 );
    auto dopr54 = plot ( convert<float>( row(adaptive, 0 ) ), convert<float>( row(adaptive, 1 ) ), "o" );
    std::cout << row ( adaptive, 1 ) << "\n";


    u = simfor::AdMiln_function_solve ( a, b, h, n, odu, x0 );


    u = simfor::adams5_function_solve ( a, b, h, n, odu, x0 );


    u = simfor::AdMltn_function_solve ( a, b, h, n, odu, x0 );



    ::matplot::legend ( { eiler, rk5, /*rkf45,*/  dopr54}, { "Метод Эйлера", "Метод Рунге-Кутта 5 порядка", "Метод Дормана-Принса 5(4) порядка"} );

    show();
    return 0;
    }

