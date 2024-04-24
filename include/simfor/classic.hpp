#ifndef SIMFOR_CLASSIC_HPP
#define SIMFOR_CLASSIC_HPP

#include <omp.h>
#include "internal/types.hpp"

namespace simfor{
    float scalar_mult(vec a, vec b);

    float scalar_mult_omp(vec a, vec b);

    float scalar_mult_mpi(vec a, vec b);
}

#endif
