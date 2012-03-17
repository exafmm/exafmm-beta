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

//! Handles all the communication of local essential trees
template<Equation equation>
class ParallelFMM : public Partition<equation> {
public:
  using Kernel<equation>::X0;                                   //!< Center of root cell
  using Kernel<equation>::R0;                                   //!< Radius of root cel
  using Kernel<equation>::Cj0;                                  //!< Begin iterator for source cells
  using Evaluator<equation>::TOPDOWN;                           //!< Flag for top down tree construction
  using Partition<equation>::print;                             //!< Print in MPI
  using Partition<equation>::XMIN;                              //!< Minimum position vector of bodies
  using Partition<equation>::XMAX;                              //!< Maximum position vector of bodies

private:
//! Get disatnce to other domain
  real getDistance(C_iter C) {
    vect dist;                                                  // Distance vector
    for( int d=0; d!=3; ++d ) {                                 // Loop over dimensions
      dist[d] = (C->X[d] + Xperiodic[d] > XMAX[d])*             //  Calculate the distance between cell C and
                (C->X[d] + Xperiodic[d] - XMAX[d])+             //  the nearest point in domain [xmin,xmax]^3
                (C->X[d] + Xperiodic[d] < XMIN[d])*             //  Take the differnece from xmin or xmax
                (C->X[d] + Xperiodic[d] - XMIN[d]);             //  or 0 if between xmin and xmax
    }                                                           // End loop over dimensions
    real R = std::sqrt(norm(dist));                             // Scalar distance
    return R;                                                   // Return scalar distance
  }

//! Determine which cells to send
  void traverseLET(C_iter C) {
    int level = int(log(MPISIZE-1) / M_LN2 / 3) + 1;            // Level of local root cell
    if( MPISIZE == 1 ) level = 0;                               // Account for serial case
    for( C_iter CC=Cj0+C->CHILD; CC!=Cj0+C->CHILD+C->NCHILD; ++CC ) {// Loop over child cells
      // send CC
      if( CC->NCHILD == 0 ) {                                   //  If cell is twig
      // send CC->LEAF
      } else {                                                  //  If cell is not twig
        bool divide = false;                                    //   Initialize logical for dividing
        if( IMAGES == 0 ) {                                     //   If free boundary condition
          Xperiodic = 0;                                        //    Set periodic coordinate offset
          real R = getDistance(CC);                             //    Get distance to other domain
          divide |= 2 * CC->R > THETA * R;                      //    Divide if the cell seems too close
        } else {                                                //   If periodic boundary condition
          for( int ix=-1; ix<=1; ++ix ) {                       //    Loop over x periodic direction
            for( int iy=-1; iy<=1; ++iy ) {                     //     Loop over y periodic direction
              for( int iz=-1; iz<=1; ++iz ) {                   //      Loop over z periodic direction
                Xperiodic[0] = ix * 2 * R0;                     //       Coordinate offset for x periodic direction
                Xperiodic[1] = iy * 2 * R0;                     //       Coordinate offset for y periodic direction
                Xperiodic[2] = iz * 2 * R0;                     //       Coordinate offset for z periodic direction
                real R = getDistance(CC);                       //       Get distance to other domain
                divide |= 2 * CC->R > THETA * R;                //       Divide if cell seems too close
              }                                                 //      End loop over z periodic direction
            }                                                   //     End loop over y periodic direction
          }                                                     //    End loop over x periodic direction
        }                                                       //   Endif for periodic boundary condition
        divide |= CC->R > (R0 / (1 << level));                  //   Divide if cell is larger than local root cell
        if( divide ) {                                          //   If cell has to be divided
          getLET(CC);                                           //    Traverse the tree further
        }                                                       //   Endif for cell division
      }                                                         //  Endif for twig
    }                                                           // End loop over child cells
  }

public:
  ParallelFMM() {}
  ~ParallelFMM() {};

//! Get local essential tree to send to each process
  void getLET(Cells &cells) {
    Cj0 = cells.begin();                                        // Set begin iterator
    C_iter Root;                                                // Root cell
    if( TOPDOWN ) {                                             // If tree was constructed top down
      Root = cells.begin();                                     //  The first cell is the root cell
    } else {                                                    // If tree was constructed bottom up
      Root = cells.end() - 1;                                   //  The last cell is the root cell
    }                                                           // Endif for tree construction
    traverseLET(Root);                                          // Traverse tree to get LET
  }

};
#endif
