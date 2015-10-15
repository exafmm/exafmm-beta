#ifndef kernel_h
#define kernel_h
#include <cmath>
#include "types.h"

namespace exafmm {
  namespace kernel {
    extern real_t eps2;                                         //!< Epslion squared
    extern complex_t wavek;                                     //!< Helmholtz wave number
    extern vec3 Xperiodic;                                      //!< Periodic coordinate offset

    void setup();                                               //!< Setup phase for kernels
    void P2P(C_iter Ci, C_iter Cj, bool mutual);                //!< P2P kernel between cells Ci and Cj
    void P2P(C_iter Ci, C_iter Cj);                             //!< P2P kernel between cells Ci and Cj for direct summation
    void P2P(C_iter C);                                         //!< P2P kernel for cell C
    void P2M(C_iter C);                                         //!< P2M kernel for cell C
    void M2M(C_iter Ci, C_iter C0);                             //!< M2M kernel for one parent cell Ci
    void M2L(C_iter Ci, C_iter Cj, bool mutual);                //!< M2L kernel between cells Ci and Cj
    void L2L(C_iter Ci, C_iter C0);                             //!< L2L kernel for one child cell Ci
    void L2P(C_iter Ci);                                        //!< L2P kernel for cell Ci
  }
}
#endif
