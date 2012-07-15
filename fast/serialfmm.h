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
#ifndef serialfmm_h
#define serialfmm_h
#include "bottomup.h"

class SerialFMM : public BottomUp {
public:
  void topdown(Bodies &bodies, Cells &cells) {
    TOPDOWN = true;
    TopDown::setDomain(bodies);
    TopDown::buildTree();
    TopDown::linkTree(bodies,cells);
    TopDown::upwardPass(cells);
  }

  void bottomup(Bodies &bodies, Cells &cells) {
    TOPDOWN = false;
    BottomUp::setDomain(bodies);
    BottomUp::buildTree(bodies,cells);
    BottomUp::linkTree(cells);
    BottomUp::upwardPass(cells);
  }

  void direct(Bodies &ibodies, Bodies &jbodies) {
    Cells cells;
    cells.resize(2);
    C_iter Ci = cells.begin(), Cj = cells.begin()+1;
    Ci->LEAF = ibodies.begin();
    Ci->NDLEAF = ibodies.size();
    Cj->LEAF = jbodies.begin();
    Cj->NDLEAF = jbodies.size();
    int prange = 0;                                             // Range of periodic images
    for( int i=0; i<IMAGES; i++ ) {                             // Loop over periodic image sublevels
      prange += int(pow(3,i));                                  //  Accumulate range of periodic images
    }                                                           // End loop over perioidc image sublevels
    for( int ix=-prange; ix<=prange; ++ix ) {                   // Loop over x periodic direction
      for( int iy=-prange; iy<=prange; ++iy ) {                 //  Loop over y periodic direction
        for( int iz=-prange; iz<=prange; ++iz ) {               //   Loop over z periodic direction
          Xperiodic[0] = ix * 2 * R0;                           //    Coordinate offset for x periodic direction
          Xperiodic[1] = iy * 2 * R0;                           //    Coordinate offset for y periodic direction
          Xperiodic[2] = iz * 2 * R0;                           //    Coordinate offset for z periodic direction
          P2P(Ci,Cj,false);                                     //    Evaluate P2P kernel
        }                                                       //   End loop over z periodic direction
      }                                                         //  End loop over y periodic direction
    }                                                           // End loop over x periodic direction
    for( B_iter B=ibodies.begin(); B!=ibodies.end(); ++B ) {
      B->TRG /= B->SRC;
    }
  }

  void evaluate(Cells &cells) {
    setRootCell(cells);
    CellQueue cellQueue;
    pushCell(ROOT,cellQueue);
    Xperiodic = 0;
    startTimer("Traverse");
    traverse(cellQueue);
    stopTimer("Traverse",printNow);
    startTimer("Downward pass");
    if( TOPDOWN ) {
      TopDown::downwardPass(cells);
    } else {
      BottomUp::downwardPass(cells);
    }
    stopTimer("Downward pass",printNow);
    if(printNow) printTreeData(cells);
  }

  void evaluate(Cells &icells, Cells &jcells) {
    setRootCell(icells,jcells);
    Pair pair(ROOT,ROOT2);
    PairQueue pairQueue;
    startTimer("Traverse");
    if( IMAGES == 0 ) {
      Xperiodic = 0;
      pairQueue.push_back(pair);
      traverse(pairQueue);
    } else {
      for( int ix=-1; ix<=1; ++ix ) {                           //  Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //   Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz ) {                       //    Loop over z periodic direction
            Xperiodic[0] = ix * 2 * R0;                         //     Coordinate offset for x periodic direction
            Xperiodic[1] = iy * 2 * R0;                         //     Coordinate offset for y periodic direction
            Xperiodic[2] = iz * 2 * R0;                         //     Coordinate offset for z periodic direction
            pairQueue.push_back(pair);                          //     Push pair to queue
            traverse(pairQueue);                                //     Traverse a pair of trees
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      traversePeriodic();                                       //  Traverse tree for periodic images
    }
    stopTimer("Traverse",printNow);
    startTimer("Downward pass");
    if( TOPDOWN ) {
      TopDown::downwardPass(icells);
    } else {
      BottomUp::downwardPass(icells);
    }
    stopTimer("Downward pass",printNow);
    if(printNow) printTreeData(icells);
  }

};

#endif
