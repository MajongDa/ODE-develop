#include <simforode/odu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cmath>
#include <iostream>
#include <fstream>
#include <matplot/matplot.h>

simfor::vec lorenz ( double t, simfor::vec a )
    {
    double sigma = 10, beta = 8/3, rho = 28;
    simfor::vec df ( a.size() );

    df ( 0 ) = -sigma*a ( 0 ) + sigma*a ( 1 );
    df ( 1 ) = rho*a ( 0 ) - a ( 1 ) - a ( 0 ) * a ( 2 );
    df ( 2 ) = -beta*a ( 2 ) + a ( 0 ) *a ( 1 );

    return df;
    }

simfor::vec Van_der_Pol(double t, simfor::vec a)
    {
    double mu = 10;
    simfor::vec oscil ( a.size() );
    oscil ( 0 ) = a ( 1 );
    oscil ( 1 ) = mu*( 1 - a(0)*a(0))*a(1) - a(0);

    return oscil;

    }

void write_matr ( const std::string file, const simfor::matr A )
    {
    std::ofstream data ( file, std::ofstream::out );
    for ( int i = 0; i < A.size1(); ++i )
        {
        for( int j = 0; j < A.size2(); ++j)
            data << A( i, j ) << " ";
        data << std::endl;
        }
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
    using namespace matplot;
    simfor::matr f;
    simfor::vec iv ( 3 );


    iv ( 0 ) = 10; iv ( 1 ) = iv ( 2 ) = 1;
    double a = 0., b = 10., h = 1e-2;



    int n = std::ceil ( ( b-a)/h );
    std::vector<double> t(n);
    for(size_t i=0; i < t.size(); ++i) t[i] = i*h;

    // f = simfor::eiler_system_solve_vector ( a, b, h, n, lorenz, iv );
    //
    // write_matr("lorenz_attr_Eul.txt", f);
    //
    // f = simfor::rk4_system_function_solve_vector ( a, b, h, n, lorenz, iv );
    //
    // write_matr("lorenz_attr_RK.txt", f);
    //
    // f = simfor::adams5_function_solve_vector( a, b, h, n, lorenz, iv );
    //
    // write_matr("lorenz_attr_AdBash.txt", f);
    //
    // f = simfor::AdMltn_function_solve_vector( a, b, h, n, lorenz, iv );
    // write_matr("lorenz_attr_AdMltn.txt", f);
    // //plot3 (  convert<double>( column( f, 0)), convert<double>( column( f, 1)), convert<double>( column( f, 2)) );





    //Van_der_Pol oscilator

    std::vector<double> u, v;

    a = 0; b = 100; h = 1e-3; n = std::ceil ( ( b-a)/h );
    t = linspace(a, b, n+1);

    //initial conditions at t=0
    simfor::vec vdp(2);
    vdp(0) = 2; vdp(1) = 0;

    hold(on);
    //
    // f = simfor::rk6_system_function_solve_vector ( a, b, h, n, Van_der_Pol, vdp );
    // u = convert<double> ( column(f, 0 ) );
    // v = convert<double> ( column(f, 1 ) );
    //
    //
    //
    // //auto rk4 = plot ( t, v );
    // plot ( v );
    //
    //
    // f = simfor::rk6_system_function_solve_vector ( a, b, h, n, Van_der_Pol, vdp );
    // u = convert<double> ( column(f, 0 ) );
    // v = convert<double> ( column(f, 1 ) );
    //
    // ///auto rk6 = plot ( u, v );
    //
    // f = simfor::AdMltn_function_solve_vector( a, b, h, n, Van_der_Pol, vdp );
    // u = convert<double> ( column(f, 0 ) );
    // v = convert<double> ( column(f, 1 ) );
    //
    // //auto AdBash = plot ( u, v );
    //
    //
    f = simfor::AdMltn_function_solve_vector( a, b, h, n, Van_der_Pol, vdp );
    std::cout << n << " " << f.size1() << "\n";
    write_matr("VdP_osc_AdMltn.txt", f);
    u = convert<double> ( column(f, 0 ) );
    v = convert<double> ( column(f, 1 ) );

    // plot ( t, v );
    // plot ( t, u );
    plot ( u, v );


    //auto AdMltn = plot ( u, v );

    f = simfor::rk54 ( a, b, Van_der_Pol, vdp, 1e-6, 1e-7, 1e-3);
    write_matr("VdP_osc_RKF45.txt", f);
    std::cout << n << " " << f.size1() << "\n";
    t = convert<double> ( column(f, 0 ) );
    u = convert<double> ( column(f, 1 ) );
    v = convert<double> ( column(f, 2 ) );
    // plot ( t, v );
    // plot ( t, u );
    plot ( u, v );

    //::matplot::legend ( {rk4, rk6, AdBash, AdMltn}, {"rk4", "rk6", "AdBash", "AdMltn"} );
    show();
    return 0;

    }
