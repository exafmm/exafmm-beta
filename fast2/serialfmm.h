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
  using Kernel<equation>::preCalculation;                       //!< Precalculate M2L translation matrix
  using Kernel<equation>::postCalculation;                      //!< Free temporary allocations
  using Evaluator<equation>::TOPDOWN;                           //!< Flag for top down tree construction
  using Evaluator<equation>::upwardPass;                        //!< Upward pass to get all multipole expansions
  using Evaluator<equation>::traverse;                          //!< Traverse tree to get interaction list
  using Evaluator<equation>::downwardPass;                      //!< Downward pass to evaluate all local expansions

public:
  SerialFMM() {
    preCalculation();
  }
  ~SerialFMM() {
    postCalculation();
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
    traverse(icells,jcells);
    downwardPass(icells);
  }

};

#endif
