#include "kernel.h"

namespace exafmm {
  class EmptyKernel : public Kernel {
  public:
    static void P2P(C_iter , C_iter , bool ) {}
    static void P2P(C_iter ) {}
    static void setup() {}
    static void P2M(C_iter ) {}
    static void M2M(C_iter , C_iter ) {}
    static void M2L(C_iter , C_iter , bool ) {}
    static void L2L(C_iter , C_iter ) {}
    static void L2P(C_iter ) {}
  };
}
