#include <simfor/odu.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <cmath>
#include <iostream>
#include <fstream>

simfor::vec lorenz ( float t, simfor::vec a )
    {
    float sigma = 10, beta = 8/3, rho = 28;
    simfor::vec df ( a.size() );

    df ( 0 ) = -sigma*a ( 0 ) + sigma*a ( 1 );
    df ( 1 ) = rho*a ( 0 ) - a ( 1 ) - a ( 0 ) * a ( 2 );
    df ( 2 ) = -beta*a ( 2 ) + a ( 0 ) *a ( 1 );

    return df;
    }

simfor::vec Van_der_Pol(float t, simfor::vec a)
    {
    float mu = 10;
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

int main()
    {
    simfor::matr f;
    simfor::vec v ( 3 );

    v ( 0 ) = 10; v ( 1 ) = v ( 2 ) = 1;
    float a = 0., b = 100., h;

    int n = 10000;
    h = 0.01 ;

    f = simfor::eiler_system_solve_vector ( a, b, h, n, lorenz, v );

    write_matr("lorenz_attr_Eul.txt", f);

    f = simfor::rk_system_function_solve_vector ( a, b, h, n, lorenz, v );

    write_matr("lorenz_attr_RK.txt", f);

    f = simfor::adams5_function_solve_vector( a, b, h, n, lorenz, v );

    write_matr("lorenz_attr_AdBash.txt", f);

//TODO find root of instability
//     f = simfor::AdMltn_function_solve_vector( a, b, h, n, lorenz, v );
//     write_matr("lorenz_attr_AdMltn.txt", f);

    //Van_der_Pol oscilator

    a = 0; b = 100; n = 100000; h = 1e-3;
    //initial conditions at t=0
    simfor::vec vdp(2);
    vdp(0) = 2; vdp(1) = 0;

    f = simfor::eiler_system_solve_vector ( a, b, h, n, Van_der_Pol, vdp );

    write_matr("VdP_osc_Eul.txt", f);

    f = simfor::rk_system_function_solve_vector ( a, b, h, n, Van_der_Pol, vdp );

    write_matr("VdP_osc_RK.txt", f);

    f = simfor::adams5_function_solve_vector( a, b, h, n, Van_der_Pol, vdp );

    write_matr("VdP_osc_AdBash.txt", f);






    return 0;

    }
