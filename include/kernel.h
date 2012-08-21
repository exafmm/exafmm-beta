#ifndef kernel_h
#define kernel_h
#include "sort.h"

class Kernel : public Sort {
protected:
  C_iter Ci0;
  C_iter Cj0;

public:
  void P2P(C_iter Ci, C_iter Cj, bool mutual) const;
  void P2P(C_iter C) const;
  void P2M(C_iter C, real_t &Rmax) const;
  void M2M(C_iter Ci, real_t &Rmax) const;
  void M2L(C_iter Ci, C_iter Cj, bool mutual) const;
  void M2P(C_iter Ci, C_iter Cj, bool mutual) const;
  void L2L(C_iter Ci) const;
  void L2P(C_iter Ci) const;
};

#endif
