/*
Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
#define KERNEL
#include "kernel.h"
#undef KERNEL

template<>
void Kernel<Laplace>::P2P(C_iter Ci, C_iter Cj) const {         // Laplace P2P kernel on CPU
  for( B_iter Bi=Ci->LEAF; Bi!=Ci->LEAF+Ci->NDLEAF; ++Bi ) {    // Loop over target bodies
    for( B_iter Bj=Cj->LEAF; Bj!=Cj->LEAF+Cj->NDLEAF; ++Bj ) {  //  Loop over source bodies
      vect dist = Bi->X - Bj->X - Xperiodic;                    //   Distance vector from source to target
      real R2 = norm(dist) + EPS2;                              //   R^2
      real invR = 1 / std::sqrt(R2);                            //   1 / R
      if( R2 == 0 ) invR = 0;                                   //   Exclude self interaction
      real invR3 = Bj->SRC * invR * invR * invR;                //   charge / R^3
      Bi->TRG[0] += Bj->SRC * invR;                             //   potential
      Bi->TRG[1] -= dist[0] * invR3;                            //   x component of force
      Bi->TRG[2] -= dist[1] * invR3;                            //   y component of force
      Bi->TRG[3] -= dist[2] * invR3;                            //   z component of force
    }                                                           //  End loop over source bodies
  }                                                             // End loop over target bodies
}

template<>
void Kernel<Laplace>::EwaldReal(C_iter Ci, C_iter Cj) const {   // Ewald real part on CPU
  for( B_iter Bi=Ci->LEAF; Bi!=Ci->LEAF+Ci->NDLEAF; ++Bi ) {    // Loop over target bodies
    for( B_iter Bj=Cj->LEAF; Bj!=Cj->LEAF+Cj->NDLEAF; ++Bj ) {  //  Loop over source bodies
      vect dist = Bi->X - Bj->X - Xperiodic;                    //   Distance vector from source to target
      for( int d=0; d<3; d++ ) {                                //   Loop over dimensions
        if( dist[d] < -R0 ) {
          dist[d] += 2 * R0;
        }
        if( dist[d] >= R0 ) {
          dist[d] -= 2 * R0;
        }
      }                                                         //   End loop over dimensions
      real R2 = norm(dist);                                     //   R^2
      if( R2 != 0 ) {                                           //   Exclude self interaction
        real R2s = R2 * ALPHA * ALPHA;                          //    (R * alpha)^2
        real Rs = std::sqrt(R2s);                               //    R * alpha
        real invRs = 1 / Rs;                                    //    1 / (R * alpha)
        real invR2s = invRs * invRs;                            //    1 / (R * alpha)^2
        real invR3s = invR2s * invRs;                           //    1 / (R * alpha)^3
        real dtmp = Bj->SRC * (M_2_SQRTPI * exp(-R2s) * invR2s + erfc(Rs) * invR3s);
        dtmp *= ALPHA * ALPHA * ALPHA;                          //    Scale temporary value
        Bi->TRG[0] += Bj->SRC * erfc(Rs) * invRs * ALPHA;       //    Ewald real potential
        Bi->TRG[1] -= dist[0] * dtmp;                           //    x component of Ewald real force
        Bi->TRG[2] -= dist[1] * dtmp;                           //    y component of Ewald real force
        Bi->TRG[3] -= dist[2] * dtmp;                           //    z component of Ewald real force
      }                                                         //   End if for self interaction
    }                                                           //  End loop over source bodies
  }                                                             // End loop over target bodies
}

template<>
void Kernel<VanDerWaals>::P2P(C_iter Ci, C_iter Cj) const {     // Van der Waals P2P kernel on CPU
  for( B_iter Bi=Ci->LEAF; Bi!=Ci->LEAF+Ci->NDLEAF; ++Bi ) {    // Loop over target bodies
    int atypei = Bi->SRC;                                       //  Atom type of target
    for( B_iter Bj=Cj->LEAF; Bj!=Cj->LEAF+Cj->NDLEAF; ++Bj ) {  //  Loop over source bodies
      int atypej = Bj->SRC;                                     //   Atom type of source
      vect dist = Bi->X - Bj->X - Xperiodic;                    //   Distance vector from source to target
      real R2 = norm(dist);                                     //   R squared
      if( R2 != 0 ) {                                           //   Exclude self interaction
        real rs = RSCALE[atypei*ATOMS+atypej];                  //    r scale
        real gs = GSCALE[atypei*ATOMS+atypej];                  //    g scale
        real R2s = R2 * rs;                                     //    R^2 * r scale
        if( R2MIN <= R2s && R2s < R2MAX ) {                     //    Exclude outlier values
          real invR2 = 1.0 / R2s;                               //     1 / R^2
          real invR6 = invR2 * invR2 * invR2;                   //     1 / R^6
          real dtmp = gs * invR6 * invR2 * (2.0 * invR6 - 1.0); //     g scale / R^2 * (2 / R^12 + 1 / R^6)
          Bi->TRG[0] += gs * invR6 * (invR6 - 1.0);             //     Van der Waals potential
          Bi->TRG[1] -= dist[0] * dtmp;                         //     x component of Van der Waals force
          Bi->TRG[2] -= dist[1] * dtmp;                         //     y component of Van der Waals force
          Bi->TRG[3] -= dist[2] * dtmp;                         //     z component of Van der Waals force
        }                                                       //    End if for outlier values
      }                                                         //   End if for self interaction 
    }                                                           //  End loop over source bodies
  }                                                             // End loop over target bodies
}
