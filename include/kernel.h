#ifndef kernel_h
#define kernel_h
#define KERNEL
#include "types.h"
#undef KERNEL

class Kernel {
protected:
  B_iter BI0;
  B_iter BIN;
  B_iter BJ0;
  B_iter BJN;
  C_iter CI;
  C_iter CJ;
public:
  void initialize();

  void P2M();

  void M2M();

  void M2L();

  void M2P();

  void P2P();

  void L2L();

  void L2P();

  void finalize();
};

#endif
