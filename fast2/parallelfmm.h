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
#ifndef parallelfmm_h
#define parallelfmm_h
#include "partition.h"

const int MPIDIM[3] = {4,2,2};

//! Handles all the communication of local essential trees
template<Equation equation>
class ParallelFMM : public Partition<equation> {
public:
  using Kernel<equation>::X0;                                   //!< Center of root cell
  using Kernel<equation>::R0;                                   //!< Radius of root cel
  using Evaluator<equation>::MAXLEVEL;                          //!< Max depth of tree
  using Partition<equation>::print;                             //!< Print in MPI

private:
  vec<3,int> Morton2Ivec(int n) const {
    int d = 0, level = 0;
    vec<3,int> Ivec = 0;
    while( n != 0 ) {
      Ivec[d] += (n % 2) * (1 << level);
      n >>= 1;
      d = (d+1) % 3;
      if( d == 0 ) level++;
    }
    return Ivec;
  }

  int Ivec2Morton(vec<3,int> Ivec, int level) const {
    int n = 0;
    for( int l=0; l!=level; ++l ) {
      n += Ivec[0] % 2 << (3 * l);
      n += Ivec[1] % 2 << (3 * l + 1);
      n += Ivec[2] % 2 << (3 * l + 2);
      Ivec[0] >>= 1;
      Ivec[1] >>= 1;
      Ivec[2] >>= 1;
    }
    return n;
  }

  void getLocalDomain(vect &xmin, vect &xmax, vec<3,int> &IvecSelf) const {
    assert(MPIDIM[0]*MPIDIM[1]*MPIDIM[2] == MPISIZE);
    int d = 0, n = MPIRANK;
    IvecSelf = Morton2Ivec(n);
    for( d=0; d!=3; ++d ) xmin[d] = X0[d] - R0 +  IvecSelf[d]      * 2 * R0 / MPIDIM[d];
    for( d=0; d!=3; ++d ) xmax[d] = X0[d] - R0 + (IvecSelf[d] + 1) * 2 * R0 / MPIDIM[d];
  }

  void getNeighborList(vec<3,int> &IvecSelf, const vect &xmin, const vect &xmax) const {
    vec<3,int> IvecNeighbor;
    int level = int(log(MPISIZE-1) / M_LN2 / 3) + 1;
    if( MPISIZE == 1 ) level = 0;
    for( int ix=-1; ix<=1; ++ix ) {
      for( int iy=-1; iy<=1; ++iy ) {
        for( int iz=-1; iz<=1; ++iz ) {
          IvecNeighbor[0] = (IvecSelf[0] + ix + MPIDIM[0]) % MPIDIM[0];
          IvecNeighbor[1] = (IvecSelf[1] + iy + MPIDIM[1]) % MPIDIM[1];
          IvecNeighbor[2] = (IvecSelf[2] + iz + MPIDIM[2]) % MPIDIM[2];
          int irank = Ivec2Morton(IvecNeighbor,level);
          if(MPIRANK==0) std::cout << irank << " " << IvecNeighbor << std::endl;
        }
      }
    }
  }

public:
  ParallelFMM() {}
  ~ParallelFMM() {};

  void commBodies(Cells &cells) {
    vec<3,int> IvecSelf = 0;
    vect xmin, xmax;
    getLocalDomain(xmin,xmax,IvecSelf);
    getNeighborList(IvecSelf,xmin,xmax);
  }

};
#endif
