#ifndef SIMFOR_TYPES_HPP
#define SIMFOR_TYPES_HPP

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <boost/numeric/ublas/operation_blocked.hpp>
#include <boost/numeric/ublas/operation.hpp>

namespace simfor
{
using matr = boost::numeric::ublas::matrix<float>;
using vec = boost::numeric::ublas::vector<float>;
using scalar_vector = boost::numeric::ublas::scalar_vector<float>;
}

#endif //SIMFOR_TYPES_HPP
