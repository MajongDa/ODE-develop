#ifndef SIMFOR_MATRIX_HPP
#define SIMFOR_MATRIX_HPP

#include <iostream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <cmath>
#include "vector.hpp"

#if __cplusplus >= 201103L

#include <initializer_list>

#define _NULL nullptr
#else
#define _NULL NULL
#endif

namespace simfor
{


template <class T>
class matrix
    {
    private:
        T *_data;
        unsigned _size1;
        unsigned _size2;
    public:

        matrix()
            {
            _size1 = _size2 = 0;
            _data = _NULL;
            }

        matrix ( unsigned size1 )
            {
            _size1 = _size2 = size1;
            _data = ( T * ) calloc ( size1*size1, sizeof ( T ) );
            }

        explicit matrix ( unsigned size1, unsigned size2 )
            {
            _size1 = size1;
            _size2 = size2;
            _data = ( T * ) calloc ( size1*size2, sizeof ( T ) );
            }

        matrix ( T *mat, unsigned size )
            {
            _size1 = _size2 = size;
            _data = ( T * ) calloc ( size*size, sizeof ( T ) );
            std::memcpy ( _data, mat, size * sizeof ( T ) );
            }

        matrix ( T *mat, unsigned size1, unsigned size2 )
            {
            _size1 = size1;
            _size2 = size2;
            _data = ( T * ) calloc ( size1*size2, sizeof ( T ) );
            std::memcpy ( _data, mat, size1*size2 * sizeof ( T ) );
            }

        matrix ( unsigned size1, unsigned size2, float a )
            {
            _size1 = size1;
            _size2 = size2;
            _data = ( T * ) malloc ( size1*size2 * sizeof ( T ) );
            for ( T *i = _data; i < _data + size1*size2; ++i )
                *i = a;
            }

        matrix ( const matrix &A )
            {
            _size1 = A._size1;
            _size2 = A._size2;

            _data = ( T * ) malloc ( A._size1 * A._size2 * sizeof ( T ) );
            std::memcpy ( _data, A._data, A._size1 * A._size2 * sizeof ( T ) );

            }

        ~matrix()
            {
            _size1 = 0;
            _size2 = 0;

            if ( _data == _NULL )
                return;
            free ( _data );
            _data = _NULL;
            }

        unsigned size1() const
            {
            return _size1;
            }

        unsigned size2() const
            {
            return _size2;
            }

        void set_row ( const vector<T> &v, unsigned row )
            {
            if ( _data == _NULL )
                throw std::invalid_argument ( "ERR_DATA_IS_NULLPTR" );
            else if ( row > _size1 )
                throw std::invalid_argument ( "ERR_ROW_INDEX" );
            else if ( v.size() > _size2 )
                throw std::invalid_argument ( "ERR_VECTOR_SIZE" );
            for ( int i = 0; i < _size2; ++i )
                {
                _data[row*_size2 + i] = v ( i );
                }

            }

        void set_submatrix ( const matrix &SUB, unsigned start1, unsigned stop1, unsigned start2, unsigned stop2 )
            {
            if ( _data == _NULL )
                throw std::invalid_argument ( "ERR_DATA_IS_NULLPTR" );
            else if ( start1 > _size1 )
                throw std::invalid_argument ( "ERR_SUBRANGE_SIZE1" );
            else if ( ( stop1 - start1 < 0 ) || ( stop1 - start1 > _size1 ) )
                throw std::invalid_argument ( "ERR_SUBRANGE_DIMENSION1" );
            else if ( start2 > _size2 )
                throw std::invalid_argument ( "ERR_SUBRANGE_SIZE2" );
            else if ( ( stop2 - start2 < 0 ) || ( stop2 - start2 > _size2 ) )
                throw std::invalid_argument ( "ERR_SUBRANGE_DIMENSION2" );
            for ( unsigned i = 0; i < stop1; ++i )
                for ( unsigned j = 0; j < stop2; ++j )
                    _data[ ( start1+i ) *_size2 + start2+j] = SUB ( i, j );
            }

        vector<T> get_row ( unsigned row )
            {
            vector<T> v ( _size2 );
            for ( int i = 0; i < _size2; ++i )
                v ( i ) = _data[row*_size2+i];
            return v;
            }

        matrix get_submatrix ( unsigned start1, unsigned stop1, unsigned start2, unsigned stop2 )
            {
            if ( _data == _NULL )
                throw std::invalid_argument ( "ERR_DATA_IS_NULLPTR" );
            else if ( start1 > _size1 )
                throw std::invalid_argument ( "ERR_SUBRANGE_SIZE1" );
            else if ( ( stop1 - start1 < 0 ) || ( stop1 - start1 > _size1 ) )
                throw std::invalid_argument ( "ERR_SUBRANGE_DIMENSION1" );
            else if ( start2 > _size2 )
                throw std::invalid_argument ( "ERR_SUBRANGE_SIZE2" );
            else if ( ( stop2 - start2 < 0 ) || ( stop2 - start2 > _size2 ) )
                throw std::invalid_argument ( "ERR_SUBRANGE_DIMENSION2" );
            matrix<T> SUB ( stop1-start1, stop2-start2 );
            for ( unsigned i = 0; i < stop1; ++i )
                for ( unsigned j = 0; j < stop2; ++j )
                    SUB ( i, j ) = _data[ ( start1+i ) *_size2 + start2+j];
            return SUB;
            }

        matrix &operator= ( const matrix &M )
            {
            _size1 = M.size1();
            _size2 = M.size2();
            _data = ( T* ) realloc ( _data, _size1 * _size2 * sizeof ( T ) );
            std::memcpy ( _data, M._data, _size1 * _size2 * sizeof ( T ) );

            return ( *this );
            }

        matrix operator- ( const matrix &M )
            {
            if ( _data == _NULL )
                {
                throw std::invalid_argument ( "ERR_LHS_DATA_IS_NULLPTR" );
                }
            else if ( M._data == _NULL )
                {
                throw std::invalid_argument ( "ERR_RHS_DATA_IS_NULLPTR" );
                }
            else if ( ( _size1 != M.size1() ) || ( _size2 != M.size2 () ) )
                {
                throw std::invalid_argument ( "ERR_DIFFERENT_DIMENSIONS" );
                }
            matrix tmp ( _size1, _size2 );
            float *endptr = _data + _size1*_size2;
            for ( float *i = _data, *j = M._data, *tmp_prt = tmp._data; i < endptr; i++, j++, tmp_prt++ )
                *tmp_prt = *i - *j;
            return tmp;
            }

        matrix operator+ ( const matrix &M )
            {
            if ( _data == _NULL )
                {
                throw std::invalid_argument ( "ERR_LHS_DATA_IS_NULLPTR" );
                }
            else if ( M._data == _NULL )
                {
                throw std::invalid_argument ( "ERR_RHS_DATA_IS_NULLPTR" );
                }
            else if ( ( _size1 != M.size1() ) || ( _size2 != M.size2 () ) )
                {
                throw std::invalid_argument ( "ERR_DIFFERENT_DIMENSIONS" );
                }
            matrix tmp ( _size1, _size2 );
            float *endptr = _data + _size1*_size2;
            for ( float *i = _data, *j = M._data, *tmp_prt = tmp._data; i < endptr; i++, j++, tmp_prt++ )
                *tmp_prt = *i + *j;
            return tmp;
            }

        matrix operator* ( float a ) const
            {
            if ( _data == _NULL )
                {
                throw std::invalid_argument ( "ERR_LHS_DATA_IS_NULLPTR" );
                }
            matrix tmp ( _size1, _size2 );
            for ( int i = 0; i < _size1*_size2; ++i )
                tmp(i) = _data[i] * a;
            return tmp;
            }


        matrix &operator-= ( const matrix &v )
            {
            matrix<T> S( *this-v );
            *this = S;
            return *this;
            }

        matrix &operator+= ( const matrix &v )
            {
            matrix<T> S( *this+v );
            *this = S;
            return *this;
            }

        matrix &operator*= ( float a )
            {
            matrix<T> S( *this*a );
            *this = S;
            return *this;
            }

        T &operator() ( unsigned i )
            {
            return _data[i];
            }

        T &operator() ( unsigned i ) const
            {
            return _data[i];
            }

        T &operator() ( unsigned i, unsigned j )
            {
            return _data[i*_size2 + j];
            }

        T &operator() ( unsigned i, unsigned j ) const
            {
            return _data[i*_size2 + j];
            }

        friend std::ostream &operator<< ( std::ostream &out, const matrix &M )
            {
            char tmp_str[20];
            out << M.size1() << " x " << M.size2() << '\n';
            for ( int i = 0; i < M.size1(); ++i )
                {
                for ( int j = 0; j < M.size2(); ++j )
                    {
                    sprintf ( tmp_str, "%6.5g ", fabs ( M ( i, j ) ) < 10e-5 ? 0. : M ( i, j ) );
                    out << tmp_str << "\t";
                    }
                out << '\n';
                }
            return out;
            }

        friend matrix operator* ( float a, const matrix &v )
            {
            return v*a;
            }
        friend vector<T> row ( matrix<T> M, unsigned row )
            {
            return M.get_row ( row );
            }


    };
}
#endif





