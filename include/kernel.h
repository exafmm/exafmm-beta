#ifndef kernel_h
#define kernel_h
#include <cmath>
#include "types.h"

namespace exafmm {
  class Kernel {
  public:
    static vec3 Xperiodic;                                      //!< Periodic coordinate offset
    static real_t eps2;                                         //!< Epslion squared
    static complex_t wavek;                                     //!< Helmholtz wave number
  };
}
#endif
