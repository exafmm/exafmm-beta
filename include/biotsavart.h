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
#ifndef biotsavart_h
#define biotsavart_h

template<>
void Kernel<BiotSavart>::P2P_CPU(C_iter Ci, C_iter Cj, vect Xperiodic) {// Biot-Savart P2P kernel on CPU
  for( B_iter Bi=Ci->LEAF; Bi!=Ci->LEAF+Ci->NDLEAF; ++Bi ) {    // Loop over target bodies
    for( B_iter Bj=Cj->LEAF; Bj!=Cj->LEAF+Cj->NDLEAF; ++Bj ) {  //  Loop over source bodies
      vect dist = Bi->X - Bj->X - Xperiodic;                    //   Distance vector from source to target
      real S2 = 2 * Bj->SRC[3] * Bj->SRC[3];                    //    2 * sigma^2
      real R2  = norm(dist) + EPS2;                             //    R^2 + epsilon^2
      real invR = 1 / std::sqrt(R2);                            //    1 / R
      if( R2 == 0 ) invR = 0;                                   //    Exclude self interaction
      real invR3 = invR * invR * invR;                          //    1 / R^3
      real RS = R2 / S2;                                        //    R^2 / (2 * simga^2)
      real cutoff = 0.25 / M_PI * invR3 * (erf( std::sqrt(RS) ) //    cutoff function
                  - std::sqrt(4 / M_PI * RS) * exp(-RS));
      Bi->TRG[0] += (dist[1] * Bj->SRC[2] - dist[2] * Bj->SRC[1]) * cutoff;// x component of curl G * cutoff
      Bi->TRG[1] += (dist[2] * Bj->SRC[0] - dist[0] * Bj->SRC[2]) * cutoff;// y component of curl G * cutoff
      Bi->TRG[2] += (dist[0] * Bj->SRC[1] - dist[1] * Bj->SRC[0]) * cutoff;// z component of curl G * cutoff
    }                                                           //  End loop over source bodies
  }                                                             // End loop over target bodies
}

#endif
