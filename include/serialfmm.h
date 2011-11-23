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

//! Serial FMM interface
template<Equation equation>
class SerialFMM : public BottomUp<equation> {
public:
  using Kernel<equation>::sortBodies;                           //!< Sort bodies according to cell index
  using Kernel<equation>::preCalculation;                       //!< Precalculate M2L translation matrix
  using Kernel<equation>::postCalculation;                      //!< Free temporary allocations
  using TreeStructure<equation>::buffer;                        //!< Buffer for MPI communication & sorting
  using TreeStructure<equation>::bodies2twigs;                  //!< Group bodies into twig cells
  using TreeStructure<equation>::twigs2cells;                   //!< Link twigs bottomup to create all cells in tree

public:
//! Constructor
  SerialFMM() : BottomUp<equation>() {
    preCalculation();
  }
//! Destructor
  ~SerialFMM() {
    postCalculation();
  }

//! Topdown tree constructor interface. Input: bodies, Output: cells
  void topdown(Bodies &bodies, Cells &cells) {
    TopDown<equation>::grow(bodies);                            // Grow tree structure topdown

    TopDown<equation>::setIndex();                              // Set index of cells

    buffer.resize(bodies.size());                               // Resize sort buffer
    sortBodies(bodies,buffer,false);                            // Sort bodies in descending order

    Cells twigs;                                                // Twigs are cells at the bottom of tree
    bodies2twigs(bodies,twigs);                                 // Turn bodies to twigs

    Cells sticks;                                               // Sticks are twigs from other processes that are not twigs in the current process
    twigs2cells(twigs,cells,sticks);                            // Turn twigs to cells
  }

//! Bottomup tree constructor interface. Input: bodies, Output: cells
  void bottomup(Bodies &bodies, Cells &cells) {
    BottomUp<equation>::setIndex(bodies);                       // Set index of cells

    buffer.resize(bodies.size());                               // Resize sort buffer
    sortBodies(bodies,buffer,false);                            // Sort bodies in descending order

/*
    prune(bodies);                                              // Prune tree structure bottomup

    BottomUp<equation>::grow(bodies);                           // Grow tree structure at bottom if necessary

    sortBodies(bodies,buffer,false);                            // Sort bodies in descending order
*/

    Cells twigs;                                                // Twigs are cells at the bottom of tree
    bodies2twigs(bodies,twigs);                                 // Turn bodies to twigs

    Cells sticks;                                               // Sticks are twigs from other processes not twigs here
    twigs2cells(twigs,cells,sticks);                            // Turn twigs to cells
  }
};

#endif
