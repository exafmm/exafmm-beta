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
#ifndef tree_h
#define tree_h
#include "bottomup.h"

template<Equation equation>
class SerialFMM : public BottomUp<equation> {
public:
  using Kernel<equation>::printNow;                             //!< Switch to print timings
  using Kernel<equation>::startTimer;                           //!< Start timer for given event
  using Kernel<equation>::stopTimer;                            //!< Stop timer for given event
  using Kernel<equation>::preCalculation;                       //!< Precalculate M2L translation matrix
  using Kernel<equation>::postCalculation;                      //!< Free temporary allocations
  using Evaluator<equation>::TOPDOWN;                           //!< Flag for top down tree construction
  using Evaluator<equation>::upwardPass;                        //!< Upward pass to get all multipoles
  using Evaluator<equation>::traverse;                          //!< Traverse tree to get interaction list

  using Kernel<equation>::P2P;

public:
  SerialFMM() {
    preCalculation();
  }
  ~SerialFMM() {
    postCalculation();
  }

  void direct(Bodies &ibodies, Bodies &jbodies) {
    Cells cells;
    cells.resize(2);
    C_iter Ci = cells.begin(), Cj = cells.begin()+1;
    Ci->LEAF = ibodies.begin();
    Ci->NDLEAF = ibodies.size();
    Cj->LEAF = jbodies.begin();
    Cj->NDLEAF = jbodies.size();
    P2P(Ci,Cj);
  }

  void topdown(Bodies &bodies, Cells &cells) {
    TOPDOWN = true;
    TopDown<equation>::setDomain(bodies);
    TopDown<equation>::buildTree();
    TopDown<equation>::linkTree(bodies,cells);
    upwardPass(cells);
  }

  void bottomup(Bodies &bodies, Cells &cells) {
    TOPDOWN = false;
    BottomUp<equation>::setDomain(bodies);
    BottomUp<equation>::buildTree(bodies,cells);
    BottomUp<equation>::linkTree(cells);
    upwardPass(cells);
  }

  void evaluate(Cells &icells, Cells &jcells) {
    startTimer("Traverse");
    traverse(icells,jcells);
    stopTimer("Traverse",printNow);
    startTimer("Downward pass");
    if( TOPDOWN ) {
      TopDown<equation>::downwardPass(icells);
    } else {
      BottomUp<equation>::downwardPass(icells);
    }
    stopTimer("Downward pass",printNow);
  }

};

#endif
