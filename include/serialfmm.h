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
#include "dataset.h"

//! Serial FMM interface
template<Equation kernelName>
class SerialFMM : public BottomUp<kernelName>, public Dataset<kernelName> {
public:
  using Kernel<kernelName>::sortBodies;                         //!< Sort bodies according to cell index
  using Kernel<kernelName>::preCalculation;                     //!< Precalculate M2L translation matrix
  using Kernel<kernelName>::postCalculation;                    //!< Free temporary allocations
  using TreeStructure<kernelName>::buffer;                      //!< Buffer for MPI communication & sorting
  using TreeStructure<kernelName>::bodies2twigs;                //!< Group bodies into twig cells
  using TreeStructure<kernelName>::twigs2cells;                 //!< Link twigs bottomup to create all cells in tree
  using Dataset<kernelName>::initSource;                        //!< Initialize source values
  using Dataset<kernelName>::initTarget;                        //!< Initialize target values

public:
//! Constructor
  SerialFMM() : BottomUp<kernelName>(), Dataset<kernelName>() {
    preCalculation();
  }
//! Destructor
  ~SerialFMM() {
    postCalculation();
  }

//! Random distribution in [-1,1]^3 cube
  void random(Bodies &bodies, int seed=1, int numSplit=1) {
    srand(seed);                                                // Set seed for random number generator
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      if( numSplit != 1 && B-bodies.begin() == int(seed*bodies.size()/numSplit) ) {// Mimic parallel dataset
        seed++;                                                 //   Mimic seed at next rank
        srand(seed);                                            //   Set seed for random number generator
      }                                                         //  Endif for mimicing parallel dataset
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        B->X[d] = rand() / (1. + RAND_MAX) * 2 * M_PI - M_PI;   //   Initialize positions
      }                                                         //  End loop over dimension
    }                                                           // End loop over bodies
    initSource(bodies);                                         // Initialize source values
    initTarget(bodies);                                         // Initialize target values
  }

//! Random distribution on r = 1 sphere
  void sphere(Bodies &bodies, int seed=1, int numSplit=1) {
    srand(seed);                                                // Set seed for random number generator
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      if( numSplit != 1 && B-bodies.begin() == int(seed*bodies.size()/numSplit) ) {// Mimic parallel dataset
        seed++;                                                 //   Mimic seed at next rank
        srand(seed);                                            //   Set seed for random number generator
      }                                                         //  Endif for mimicing parallel dataset
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        B->X[d] = rand() / (1. + RAND_MAX) * 2 - 1;             //   Initialize positions
      }                                                         //  End loop over dimension
      real r = std::sqrt(norm(B->X));                           //  Distance from center
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        B->X[d] /= r * 1.1;                                     //   Normalize positions
      }                                                         //  End loop over dimension
    }                                                           // End loop over bodies
    initSource(bodies);                                         // Initialize source values
    initTarget(bodies);                                         // Initialize target values
  }

//! Uniform distribution on [-1,1]^3 lattice (for debugging)
  void lattice(Bodies &bodies) {
    int level = int(log(bodies.size()*MPISIZE+1.)/M_LN2/3);     // Level of tree
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      int d = 0, l = 0;                                         //  Initialize dimension and level
      int index = MPIRANK * bodies.size() + (B-bodies.begin()); //  Set index of body iterator
      vec<3,int> nx = 0;                                        //  Initialize 3-D cell index
      while( index != 0 ) {                                     //  Deinterleave bits while index is nonzero
        nx[d] += (index % 2) * (1 << l);                        //   Add deinterleaved bit to 3-D cell index
        index >>= 1;                                            //   Right shift the bits
        d = (d+1) % 3;                                          //   Increment dimension
        if( d == 0 ) l++;                                       //   If dimension is 0 again, increment level
      }                                                         //  End while loop for deinterleaving bits
      for( d=0; d!=3; ++d ) {                                   //  Loop over dimensions
        B->X[d] = -1 + (2 * nx[d] + 1.) / (1 << level);         //   Calculate cell center from 3-D cell index
      }                                                         //  End loop over dimensions
    }                                                           // End loop over bodies
    initSource(bodies);                                         // Initialize source values
    initTarget(bodies);                                         // Initialize target values
  }

//! Topdown tree constructor interface. Input: bodies, Output: cells
  void topdown(Bodies &bodies, Cells &cells) {
    TopDown<kernelName>::grow(bodies);                          // Grow tree structure topdown

    TopDown<kernelName>::setIndex();                            // Set index of cells

    buffer.resize(bodies.size());                               // Resize sort buffer
    sortBodies(bodies,buffer,false);                            // Sort bodies in descending order

    Cells twigs;                                                // Twigs are cells at the bottom of tree
    bodies2twigs(bodies,twigs);                                 // Turn bodies to twigs

    Cells sticks;                                               // Sticks are twigs from other processes that are not twigs in the current process
    twigs2cells(twigs,cells,sticks);                            // Turn twigs to cells
  }

//! Bottomup tree constructor interface. Input: bodies, Output: cells
  void bottomup(Bodies &bodies, Cells &cells) {
    BottomUp<kernelName>::setIndex(bodies);                     // Set index of cells

    buffer.resize(bodies.size());                               // Resize sort buffer
    sortBodies(bodies,buffer,false);                            // Sort bodies in descending order

/*
    prune(bodies);                                              // Prune tree structure bottomup

    BottomUp<kernelName>::grow(bodies);                         // Grow tree structure at bottom if necessary

    sortBodies(bodies,buffer,false);                            // Sort bodies in descending order
*/

    Cells twigs;                                                // Twigs are cells at the bottom of tree
    bodies2twigs(bodies,twigs);                                 // Turn bodies to twigs

    Cells sticks;                                               // Sticks are twigs from other processes not twigs here
    twigs2cells(twigs,cells,sticks);                            // Turn twigs to cells
  }
};

#endif
