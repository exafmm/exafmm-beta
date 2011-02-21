#ifndef kernel_h
#define kernel_h
#define KERNEL
#include "types.h"
#undef KERNEL

class Kernel {
public:
  void setDevice();

  void P2P(B_iter B0, B_iter BN);

  void P2P(B_iter Bi0, B_iter BiN, B_iter Bj0, B_iter BjN);

  void P2M(C_iter C);

  void M2M(C_iter CI, C_iter CJ);

  void M2L(C_iter CI, C_iter CJ);

  void L2L(C_iter CI, C_iter CJ);

  void L2P(C_iter C);

  void M2P(C_iter CI, C_iter CJ);
};

#endif
