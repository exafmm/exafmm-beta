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
void Kernel<Laplace>::P2P(C_iter C) const {
  B_iter B = C->BODY;
  int n = C->NDBODY;
  int i = 0;
  for ( ; i<n; i++) {
    real pot = 0;
    vect acc = 0;
    for (int j=i+1; j<n; j++) {
      vect dX = B[i].X - B[j].X;
      real R2 = norm(dX) + EPS2;
      if (R2 != 0) {
        real invR2 = 1.0 / R2;
        real invR = B[i].SRC * B[j].SRC * sqrtf(invR2);
        dX *= invR2 * invR;
        pot += invR;
        acc += dX;
        B[j].TRG[0] += invR;
        B[j].TRG[1] += dX[0];
        B[j].TRG[2] += dX[1];
        B[j].TRG[3] += dX[2];
      }
    }
    B[i].TRG[0] += pot;
    B[i].TRG[1] -= acc[0];
    B[i].TRG[2] -= acc[1];
    B[i].TRG[3] -= acc[2];
  }
}

template<>
void Kernel<Laplace>::P2P(C_iter Ci, C_iter Cj) const {         // Laplace P2P kernel on CPU
#if 1
  B_iter Bi = Ci->BODY;
  B_iter Bj = Cj->BODY;
  int ni = Ci->NDBODY;
  int nj = Cj->NDBODY;
  int i = 0;
  for ( ; i<ni; i++) {
    real pot = 0;
    vect acc = 0;
    for (int j=0; j<nj; j++) {
      vect dX = Bi[i].X - Bj[j].X - Xperiodic;
      real R2 = norm(dX) + EPS2;
      if (R2 != 0) {
        real invR2 = 1.0f / R2;
        real invR = Bi[i].SRC * Bj[j].SRC * sqrtf(invR2);
        dX *= invR2 * invR;
        pot += invR;
        acc += dX;
      }
    }
    Bi[i].TRG[0] += pot;
    Bi[i].TRG[1] -= acc[0];
    Bi[i].TRG[2] -= acc[1];
    Bi[i].TRG[3] -= acc[2];
  }
#else
  for( B_iter Bi=Ci->BODY; Bi!=Ci->BODY+Ci->NDBODY; ++Bi ) {    // Loop over target bodies
    real P0 = 0;                                                //  Initialize potential
    vect F0 = 0;                                                //  Initialize force
    for( B_iter Bj=Cj->BODY; Bj!=Cj->BODY+Cj->NDBODY; ++Bj ) {  //  Loop over source bodies
      vect dist = Bi->X - Bj->X - Xperiodic;                    //   Distance vector from source to target
      real R2 = norm(dist) + EPS2;                              //   R^2
      real invR2 = 1.0 / R2;                                    //   1 / R^2
      if( R2 == 0 ) invR2 = 0;                                  //   Exclude self interaction
      real invR = Bi->SRC * Bj->SRC * std::sqrt(invR2);         //   potential
      dist *= invR2 * invR;                                     //   force
      P0 += invR;                                               //   accumulate potential
      F0 += dist;                                               //   accumulate force
    }                                                           //  End loop over source bodies
    Bi->TRG[0] += P0;                                           //  potential
    Bi->TRG[1] -= F0[0];                                        //  x component of force
    Bi->TRG[2] -= F0[1];                                        //  y component of force
    Bi->TRG[3] -= F0[2];                                        //  z component of force
  }                                                             // End loop over target bodies
#endif
}

template<>
void Kernel<VanDerWaals>::P2P(C_iter Ci, C_iter Cj) const {     // Van der Waals P2P kernel on CPU
  for( B_iter Bi=Ci->BODY; Bi!=Ci->BODY+Ci->NDBODY; ++Bi ) {    // Loop over target bodies
    int atypei = int(Bi->SRC);                                  //  Atom type of target
    for( B_iter Bj=Cj->BODY; Bj!=Cj->BODY+Cj->NDBODY; ++Bj ) {  //  Loop over source bodies
      int atypej = int(Bj->SRC);                                //   Atom type of source
      vect dist = Bi->X - Bj->X - Xperiodic;                    //   Distance vector from source to target
      real R2 = norm(dist);                                     //   R squared
      if( R2 != 0 ) {                                           //   Exclude self interaction
        real rs = RSCALE[atypei*ATOMS+atypej];                  //    r scale
        real gs = GSCALE[atypei*ATOMS+atypej];                  //    g scale
        real R2s = R2 * rs;                                     //    R^2 * r scale
        if( R2MIN <= R2 && R2 < R2MAX ) {                       //    Exclude outlier values
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
