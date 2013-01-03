#ifndef kernel_h
#define kernel_h
#include <cmath>
#include "types.h"

class Kernel {
private:
  real_t    *factorial;                                         //!< Factorial
  real_t    *prefactor;                                         //!< \f$ \sqrt{ \frac{(n - |m|)!}{(n + |m|)!} } \f$
  real_t    *Anm;                                               //!< \f$ (-1)^n / \sqrt{ \frac{(n + m)!}{(n - m)!} } \f$
  complex_t *Cnm;                                               //!< M2L translation matrix \f$ C_{jn}^{km} \f$

protected:
  vec3 Xperiodic;
  C_iter Ci0;                                                   //!< Begin iterator for target cells
  C_iter Cj0;                                                   //!< Begin iterator for source cells

private:
  void preCalculation();
  void postCalculation();

public:
//! Constructor
  Kernel() {
    preCalculation();
  }
//! Destructor
  ~Kernel() {
    postCalculation();
  }
  void P2P(C_iter Ci, C_iter Cj, bool mutual) const;
  void P2P(C_iter C) const;
#if EVAL_ERROR_KAHAN
  void P2PKahan(C_iter Ci, C_iter Cj) const;
#endif
  void P2M(C_iter C, real_t &Rmax) const;
  void M2M(C_iter Ci, real_t &Rmax) const;
  void M2L(C_iter Ci, C_iter Cj, bool mutual) const;
  void L2L(C_iter Ci) const;
  void L2P(C_iter Ci) const;
};

#endif
