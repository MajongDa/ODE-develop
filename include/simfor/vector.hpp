#ifndef SIMFOR_VECTOR_HPP
#define SIMFOR_VECTOR_HPP

#include <istream>
#include <fstream>
#include <cstring>
#include <sstream>
#include <iomanip>
#include <cstdio>
#include <cmath>

#if __cplusplus >= 201103L

#include <initializer_list>

#define _NULL nullptr
#else
#define _NULL NULL
#endif

namespace simfor
{


template <class T>
class vector
    {
    private:
        T *_data;
        unsigned _size;
    public:

        vector()
            {
            _size = 0;
            _data = _NULL;
            }

        explicit vector ( unsigned size )
            {
            _size = size;
            _data = ( T * ) calloc ( size, sizeof ( T ) );
            }

        vector ( T *mat, unsigned size )
            {
            _size = size;
            _data = ( T * ) calloc ( size, sizeof ( T ) );
            std::memcpy ( _data, mat, size * sizeof ( T ) );
            }

        vector ( unsigned size, float a )
            {
            _size = size;
            _data = ( T * ) malloc ( size * sizeof ( T ) );
            for ( T *i = _data; i < _data + size; ++i )
                *i = a;
            }

        explicit vector ( std::ifstream &file )
            {
            float a;
            T *x = _NULL;
            unsigned n1;
            if ( file.is_open() )
                {
                file >> n1;
                if ( file.eof() && ( n1 != 0 ) )
                    {
                    throw std::invalid_argument ( "ERR_NOT_VALID_FILE_DATA" );
                    }
                x = new float[n1];
                for ( unsigned j = 0; j < n1; j++ )
                    {
                    file >> a;
                    x[j] = a;
                    if ( file.eof() && ( j + 1 != n1 ) )
                        {
                        throw std::invalid_argument ( "ERR_NOT_VALID_FILE_DATA" );
                        }
                    }

                _size = n1;
                _data = ( T * ) malloc ( n1 * sizeof ( T ) );
                std::memcpy ( _data, x, n1 * sizeof ( T ) );
                }
            else
                {
                throw std::invalid_argument ( "ERR_FILE_NOT_OPEN" );
                }
            }

        vector ( const vector &A )
            {
            _size = A._size;
            _data = ( T * ) malloc ( A.size() * sizeof ( T ) );
            std::memcpy ( _data, A._data, A.size() * sizeof ( T ) );

            }

        ~vector()
            {
            _size = 0;
            free ( _data );
            _data = _NULL;
            }




        unsigned size() const
            {
            return _size;
            }

        void set_subrange ( const vector vec, unsigned start, unsigned stop )
            {
            if ( _data == _NULL )
                throw std::invalid_argument ( "ERR_MISSING_FIRST_OPERAND" );
            else if ( stop - start > _size )
                throw std::invalid_argument ( "ERR_SUBRANGE_SIZE" );
            for ( unsigned i = 0; i < stop; ++i )
                _data[start + i] = vec[i];
            }

        vector get_subrange ( unsigned start, unsigned stop )
            {
            if ( _data == _NULL )
                {
                throw std::invalid_argument ( "ERR_DATA_IS_NULLPTR" );
                }
            else if ( start > _size )
                {
                throw std::invalid_argument ( "ERR_NON_EXISTENT_SUBVECTOR" );
                }
            else if ( stop-start > _size )
                {
                throw std::invalid_argument ( "ERR_NON_EXISTENT_MATRIX_ROW" );
                }
            vector<T> v ( stop-start );
            for ( int i = 0; i < stop; i++ )
                {
                v[i] = _data[start+i];
                }
            return v;
            }

        vector &operator= ( const vector &v )
            {
            _size = v._size;
            _data = ( T* ) realloc ( _data, v._size * sizeof ( T ) );
            std::memcpy ( _data, v._data, v._size * sizeof ( T ) );
            return *this;
            }

        vector operator- ( const vector &v )
            {
            if ( _data == _NULL )
                {
                throw std::invalid_argument ( "ERR_LHS_DATA_IS_NULLPTR" );
                }
            else if ( v._data == _NULL )
                {
                throw std::invalid_argument ( "ERR_RHS_DATA_IS_NULLPTR" );
                }
            else if ( _size != v._size )
                {
                throw std::invalid_argument ( "ERR_DIFFERENT_DIMENSIONS" );
                }
            vector tmp ( _size );
            float *endptr = _data + _size;
            for ( int i = 0; i < _size; i++ )
                tmp ( i ) = _data[i]-v ( i );
            return tmp;
            }

        vector operator+ ( const vector &v )
            {
            if ( _data == _NULL )
                {
                throw std::invalid_argument ( "ERR_LHS_DATA_IS_NULLPTR" );
                }
            else if ( v._data == _NULL )
                {
                throw std::invalid_argument ( "ERR_RHS_DATA_IS_NULLPTR" );
                }
            else if ( _size != v._size )
                {
                throw std::invalid_argument ( "ERR_DIFFERENT_DIMENSIONS" );
                }
            vector tmp ( _size );
            for ( int i = 0; i < _size; i++ )
                tmp ( i ) = _data[i]+v ( i );
            return tmp;
            }

        template < class C>
        vector operator* ( const C &a )
            {
            if ( _data == _NULL )
                {
                throw std::invalid_argument ( "ERR_LHS_DATA_IS_NULLPTR" );
                }
            vector tmp ( _size );
            for (int i = 0; i < _size; ++i)
                tmp(i) = _data[i] *a;
            return tmp;
            }

        template < class C>
        vector operator* ( const C &a ) const
            {
            if ( _data == _NULL )
                {
                throw std::invalid_argument ( "ERR_LHS_DATA_IS_NULLPTR" );
                }
            vector tmp ( _size );
            for (int i = 0; i < _size; ++i)
                tmp(i) = _data[i] *a;
            return tmp;
            }


        vector &operator-= ( const vector &v )
            {
            vector<T> s( *this-v );
            *this = s;
            return *this;
            }

        vector &operator+= ( const vector &v )
            {
            vector<T> s( *this+v );
            *this = s;
            return *this;
            }

        vector &operator*= ( float a )
            {
            vector<T> s( *this*a );
            *this = s;
            return *this;
            }

        // T *operator[] ( unsigned i )
        //     {
        //     return &_data[i];
        //     }
        //
        // T *operator[] ( unsigned i ) const
        //     {
        //     return &_data[i];
        //     }

        T &operator() ( unsigned i )
            {
            return _data[i];
            }

        T operator() ( unsigned i ) const
            {
            return _data[i];
            }

        friend std::ostream &operator<< ( std::ostream &out, const vector &v )
            {
            out << v._size << '\n';
            for ( int i = 0; i < v._size; ++i )
                {
                out << std::setprecision ( 3 ) << v ( i ) << "\t";
                }
            out << '\n';
            return out;
            }
        template < class C>
        friend vector<T> operator* ( const C &a, const vector<T> &v )
            {
            return v * a;
            }
    };
}
#endif


