#ifndef SIMFOR_ELEMENTARY_HPP
#define SIMFOR_ELEMENTARY_HPP

#include "internal/types.hpp"
#include <omp.h>
#include <boost/mpi.hpp>

namespace simfor {

float v_norm2(vec a, vec b);

float v_norm2_omp(vec a, vec b);

float m_norm4(matr A);

float m_norm4_omp(matr A);

float psum(vec a);

float psum_omp(vec a);

} // namespace simfor

#endif
