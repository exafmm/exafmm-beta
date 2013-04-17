#ifndef kernel_h
#define kernel_h
#include <cmath>
#include "types.h"

class Kernel {
protected:
  vec3 Xperiodic;
  C_iter Ci0;                                                   //!< Begin iterator for target cells
  C_iter Cj0;                                                   //!< Begin iterator for source cells

public:
  void P2P(C_iter Ci, C_iter Cj, bool mutual) const;
  void P2P(C_iter C) const;
  void P2M(C_iter C) const;
  void M2M(C_iter Ci) const;
  void M2L(C_iter Ci, C_iter Cj, bool mutual) const;
  void L2L(C_iter Ci) const;
  void L2P(C_iter Ci) const;
};

#endif
