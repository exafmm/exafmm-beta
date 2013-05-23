#ifndef kernel_h
#define kernel_h
#include <cmath>
#include "types.h"

class Kernel {
 protected:
  vec3 Xperiodic;                                               //!< Coordinate offset for periodic B.C.

 private:
//! Calculate Bmax
  real_t getBmax(vec3 const &X, C_iter C) const {
    real_t rad = C->R;                                          // Radius of cell
    real_t dx = rad + std::abs(X[0]-C->X[0]);                   // Add x distance from center of mass
    real_t dy = rad + std::abs(X[1]-C->X[1]);                   // Add y distance from center of mass
    real_t dz = rad + std::abs(X[2]-C->X[2]);                   // Add z distance from center of mass
    return std::sqrt(dx * dx + dy * dy + dz * dz);              // Return scalar distance
  }

 protected:
//! Set center of expansion to center of mass
  void setBmax(C_iter C, C_iter C0) const {
    real_t m = 0;                                               // Initialize mass
    vec3 X = 0;                                                 // Initialize coordinates
    for (B_iter B=C->BODY; B!=C->BODY+C->NCBODY; B++) {         // Loop over bodies
      m += B->SRC;                                              //  Accumulate mass
      X += B->X * B->SRC;                                       //  Accumulate dipole
    }                                                           // End loop over bodies
    for (C_iter c=C0+C->CHILD; c!=C0+C->CHILD+C->NCHILD; c++) { // Loop over child cells
      m += std::abs(c->M[0]);                                   //  Accumulate mass
      X += c->X * std::abs(c->M[0]);                            //  Accumulate dipole
    }                                                           // End loop over child cells
    X /= m;                                                     // Center of mass
    C->R = getBmax(X,C);                                        // Use Bmax as cell radius
  }

 public:
  Kernel() : Xperiodic(0) {}                                    //!< Constructor
  void P2P(C_iter Ci, C_iter Cj, bool mutual) const;            //!< P2P kernel between cells Ci and Cj
  void P2P(C_iter C) const;                                     //!< P2P kernel for cell C
  void P2M(C_iter C) const;                                     //!< P2M kernel for cell C
  void M2M(C_iter Ci, C_iter C0) const;                         //!< M2M kernel for one parent cell Ci
  void M2L(C_iter Ci, C_iter Cj, bool mutual) const;            //!< M2L kernel between cells Ci and Cj
  void L2L(C_iter Ci, C_iter C0) const;                         //!< L2L kernel for one child cell Ci
  void L2P(C_iter Ci) const;                                    //!< L2P kernel for cell Ci
};
#endif
