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
  using Kernel<equation>::R0;                                   //!< Radius of root cell
  using Kernel<equation>::preCalculation;                       //!< Precalculate M2L translation matrix
  using Kernel<equation>::postCalculation;                      //!< Free temporary allocations
  using Evaluator<equation>::TOPDOWN;                           //!< Flag for top down tree construction
  using Evaluator<equation>::Iperiodic;                         //!< Periodic image flag (using each bit for images)
  using Evaluator<equation>::Icenter;                           //!< Periodic image flag at center
  using Evaluator<equation>::upwardPass;                        //!< Upward pass to get all multipole expansions
  using Evaluator<equation>::upwardPeriodic;                    //!< Upward pass for periodic images             
  using Evaluator<equation>::traverse;                          //!< Traverse tree to get interaction list
  using Evaluator<equation>::traversePeriodic;                  //!< Traverse tree for periodic images     
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
    if( IMAGES != 0 ) upwardPeriodic(jcells);
    startTimer("Traverse");
    if( IMAGES == 0 ) {
      Iperiodic = Icenter;
      Xperiodic = 0;
      traverse(icells,jcells);
    } else {
      int I = 0;                                                //  Initialize index of periodic image
      for( int ix=-1; ix<=1; ++ix ) {                           //  Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //   Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz, ++I ) {                  //    Loop over z periodic direction
            Iperiodic = 1 << I;                                 //     Set periodic image flag
            Xperiodic[0] = ix * 2 * R0;                         //     Coordinate offset for x periodic direction
            Xperiodic[1] = iy * 2 * R0;                         //     Coordinate offset for y periodic direction
            Xperiodic[2] = iz * 2 * R0;                         //     Coordinate offset for z periodic direction
            traverse(icells,jcells);                            //     Traverse a pair of trees
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
    }
    stopTimer("Traverse",printNow);
    if( IMAGES != 0 ) traversePeriodic(icells,jcells);
    downwardPass(icells);
  }

};

#endif
