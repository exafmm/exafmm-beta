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
#include "evaluator.h"

//! Base class for tree structure
template<Equation equation>
class TreeStructure : public Evaluator<equation> {
public:
  Bodies buffer;                                                //!< Buffer for MPI communication & sorting

  using Kernel<equation>::printNow;                             //!< Switch to print timings
  using Kernel<equation>::startTimer;                           //!< Start timer for given event
  using Kernel<equation>::stopTimer;                            //!< Stop timer for given event
  using Kernel<equation>::sortCells;                            //!< Sort cells according to cell index
  using Kernel<equation>::X0;                                   //!< Center of root cell
  using Kernel<equation>::R0;                                   //!< Radius of root cell
  using Kernel<equation>::NP2P;                                 //!< Number of P2P kernel calls
  using Kernel<equation>::NM2P;                                 //!< Number of M2P kernel calls
  using Kernel<equation>::NM2L;                                 //!< Number of M2L kernel calls
  using Evaluator<equation>::getLevel;                          //!< Get level from cell index
  using Evaluator<equation>::timeKernels;                       //!< Time all kernels for auto-tuning
  using Evaluator<equation>::upwardPeriodic;                    //!< Upward phase for periodic cells
  using Evaluator<equation>::traverse;                          //!< Traverse tree to get interaction list
  using Evaluator<equation>::traversePeriodic;                  //!< Traverse tree for periodic images
  using Evaluator<equation>::neighbor;                          //!< Traverse source tree to get neighbor list
  using Evaluator<equation>::evalP2M;                           //!< Evaluate P2M kernel
  using Evaluator<equation>::evalM2M;                           //!< Evaluate M2M kernel
  using Evaluator<equation>::evalM2L;                           //!< Evaluate M2L kernel
  using Evaluator<equation>::evalM2P;                           //!< Evaluate M2P kernel
  using Evaluator<equation>::evalP2P;                           //!< Evaluate P2P kernel
  using Evaluator<equation>::evalL2L;                           //!< Evaluate L2L kernel
  using Evaluator<equation>::evalL2P;                           //!< Evaluate L2P kernel
  using Evaluator<equation>::evalEwaldReal;                     //!< Evaluate Ewald real part
  using Evaluator<equation>::EwaldWave;                         //!< Evalaute Ewald wave part

private:
//! Get parent cell index from current cell index
  bigint getParent(bigint index) {
    int level = getLevel(index);                                // Get level from cell index
    bigint cOff = ((1 << 3 *  level   ) - 1) / 7;               // Cell index offset of current level
    bigint pOff = ((1 << 3 * (level-1)) - 1) / 7;               // Cell index offset of parent level
    bigint i = ((index-cOff) >> 3) + pOff;                      // Cell index of parent cell
    return i;                                                   // Return cell index of parent cell
  }

//! Merge sticks with cells (levelwise)
  void unique(Cells &cells, Cells &sticks, int begin, int &end) {
    int c_old = begin;                                          // Initialize old cell counter
    for( int c=begin; c!=end; ++c ) {                           // Loop over cells in level
      if( cells[c].ICELL != cells[c_old].ICELL ) {              //  If current cell index is different from previous
        c_old = c;                                              //   Update old cell counter
      } else if( c != c_old ) {                                 //  If cell index is repeated
        if( cells[c].NCHILD != 0 ) {                            //   Stick-cell collision
          cells[c_old].NCHILD = cells[c].NCHILD;                //    Copy number of children
          cells[c_old].NDLEAF = cells[c].NDLEAF;                //    Copy number of leafs
          cells[c_old].PARENT = cells[c].PARENT;                //    Copy parent link
          cells[c_old].CHILD = cells[c].CHILD;                  //    Copy child link
          cells[c_old].LEAF = cells[c].LEAF;                    //    Copy iterator of first leaf
          sticks.push_back(cells[c_old]);                       //    Push stick into vector
        }                                                       //   Endif for collision type
        cells[c_old].M += cells[c].M;                           //   Accumulate multipole
        cells.erase(cells.begin()+c);                           //   Erase colliding cell
        c--;                                                    //   Decrement counter to account for erase
        end--;                                                  //   Decrement end to account for erase
      }                                                         //  Endif for repeated cell index
    }                                                           // End loop over cells in level
  }

//! Form parent-child mutual link
  void linkParent(Cells &cells, int &begin, int &end) {
    Cell parent;                                                // Parent cell
    Cells parents;                                              // Parent cell vector;
    int oldend = end;                                           // Save old end counter
    parent.ICELL = getParent(cells[begin].ICELL);               // Set cell index
    parent.M = 0;                                               // Initialize multipole coefficients
    parent.L = 0;                                               // Initlalize local coefficients
    parent.NDLEAF = parent.NCHILD = 0;                          // Initialize NDLEAF & NCHILD
    parent.LEAF = cells[begin].LEAF;                            // Set pointer to first leaf
    parent.CHILD = begin;                                       // Link to child
    getCenter(parent);                                          // Set cell center and radius
    for( int i=begin; i!=oldend; ++i ) {                        // Loop over cells at this level
      if( getParent(cells[i].ICELL) != parent.ICELL ) {         //  If it belongs to a new parent cell
        cells.push_back(parent);                                //   Push cells into vector
        end++;                                                  //   Increment cell counter
        parent.ICELL = getParent(cells[i].ICELL);               //   Set cell index
        parent.M = 0;                                           //   Initialize multipole coefficients
        parent.L = 0;                                           //   Initialize local coefficients
        parent.NDLEAF = parent.NCHILD = 0;                      //   Initialize NDLEAF & NCHILD
        parent.LEAF = cells[i].LEAF;                            //   Set pointer to first leaf
        parent.CHILD = i;                                       //   Link to child
        getCenter(parent);                                      //   Set cell center and radius
      }                                                         //  Endif for new parent cell
      for( int c=0; c!=cells[i].NCHILD; ++c ) {                 //  Loop over child cells
        (cells.begin()+cells[i].CHILD+c)->PARENT = i;           //   Link child to current
      }                                                         //  End loop over child cells
      cells[i].PARENT = end;                                    //  Link to current to parent
      parent.NDLEAF += cells[i].NDLEAF;                         //  Add nleaf of child to parent
      parents.push_back(parent);                                //  Push parent cell into vector
      parent = parents.back();                                  //  Copy information from vector
      parents.pop_back();                                       //  Pop parent cell from vector
      parent.NCHILD++;                                          //  Increment child counter
    }                                                           // End loop over cells at this level
    cells.push_back(parent);                                    // Push cells into vector
    end++;                                                      // Increment cell counter
    begin = oldend;                                             // Set new begin index to old end index
  }

protected:
//! Get cell center and radius from cell index
  void getCenter(Cell &cell) {
    int level = getLevel(cell.ICELL);                           // Get level from cell index
    bigint index = cell.ICELL - ((1 << 3*level) - 1) / 7;       // Subtract cell index offset of current level
    cell.R = R0 / (1 << level);                                 // Cell radius
    int d = level = 0;                                          // Initialize dimension and level
    vec<3,int> nx = 0;                                          // Initialize 3-D cell index
    while( index != 0 ) {                                       // Deinterleave bits while index is nonzero
      nx[d] += (index % 2) * (1 << level);                      //  Add deinterleaved bit to 3-D cell index
      index >>= 1;                                              //  Right shift the bits
      d = (d+1) % 3;                                            //  Increment dimension
      if( d == 0 ) level++;                                     //  If dimension is 0 again, increment level
    }                                                           // End while loop for deinterleaving bits
    for( d=0; d!=3; ++d ) {                                     // Loop over dimensions
      cell.X[d] = (X0[d]-R0) + (2 *nx[d] + 1) * cell.R;         //  Calculate cell center from 3-D cell index
    }                                                           // End loop over dimensions
  }

//! Group bodies into twig cells
  void bodies2twigs(Bodies &bodies, Cells &twigs) {
    startTimer("Bodies2twigs ");                                // Start timer
    int nleaf = 0;                                              // Initialize number of leafs
    bigint index = bodies[0].ICELL;                             // Initialize cell index
    B_iter firstLeaf = bodies.begin();                          // Initialize body iterator for first leaf
    Cell cell;                                                  // Cell structure
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      if( B->ICELL != index ) {                                 //  If it belongs to a new cell
        cell.NDLEAF = nleaf;                                    //   Set number of leafs
        cell.NCHILD = 0;                                        //   Set number of child cells
        cell.ICELL  = index;                                    //   Set cell index
        cell.LEAF   = firstLeaf;                                //   Set pointer to first leaf
        getCenter(cell);                                        //   Set cell center and radius
        twigs.push_back(cell);                                  //   Push cells into vector
        firstLeaf = B;                                          //   Set new first leaf
        nleaf = 0;                                              //   Reset number of bodies
        index = B->ICELL;                                       //   Set new cell
      }                                                         //  Endif for new cell
      nleaf++;                                                  //  Increment body counter
    }                                                           // End loop over bodies
    cell.NDLEAF = nleaf;                                        // Set number of leafs
    cell.NCHILD = 0;                                            // Set number of child cells
    cell.ICELL  = index;                                        // Set cell index
    cell.LEAF   = firstLeaf;                                    // Set pointer to first leaf
    getCenter(cell);                                            // Set cell center and radius
    twigs.push_back(cell);                                      // Push cells into vector
    stopTimer("Bodies2twigs ",printNow);                        // Stop timer & print
    evalP2M(twigs);                                             // Evaluate all P2M kernels
  }

//! Link twigs bottomup to create all cells in tree
  void twigs2cells(Cells &twigs, Cells &cells, Cells &sticks) {
    int begin = 0, end = 0;                                     // Initialize range of cell vector
    int level = getLevel(twigs.back().ICELL);                   // Initialize level of tree
    startTimer("Sort resize  ");                                // Start timer
    Cells cbuffer;                                              // Sort buffer for cells
    cbuffer.resize(2*twigs.size());                             // Resize sort buffer for cells
    stopTimer("Sort resize  ");                                 // Stop timer
    while( !twigs.empty() ) {                                   // Keep poppig twigs until the vector is empty
      while( getLevel(twigs.back().ICELL) != level ) {          //  While cell belongs to a higher level
        sortCells(cells,cbuffer,false,begin,end);               //   Sort cells at this level
        startTimer("Twigs2cells  ");                            //   Start timer
        unique(cells,sticks,begin,end);                         //   Get rid of duplicate cells
        linkParent(cells,begin,end);                            //   Form parent-child mutual link
        level--;                                                //   Go up one level
        stopTimer("Twigs2cells  ");                             //   Stop timer
      }                                                         //  End while for higher level
      startTimer("Twigs2cells  ");                              //  Start timer
      cells.push_back(twigs.back());                            //  Push cells into vector
      twigs.pop_back();                                         //  Pop twigs from vector
      end++;                                                    //  Increment cell counter
      stopTimer("Twigs2cells  ");                               //  Stop timer
    }                                                           // End while for popping twigs
    for( int l=level; l>0; --l ) {                              // Once all the twigs are done, do the rest
      sortCells(cells,cbuffer,false,begin,end);                 //  Sort cells at this level
      startTimer("Twigs2cells  ");                              //  Start timer
      unique(cells,sticks,begin,end);                           //  Get rid of duplicate cells
      linkParent(cells,begin,end);                              //  Form parent-child mutual link
      stopTimer("Twigs2cells  ");                               //  Stop timer
    }                                                           // End loop over levels
    startTimer("Twigs2cells  ");                                // Start timer
    unique(cells,sticks,begin,end);                             // Just in case there is a collision at root
    stopTimer("Twigs2cells  ",printNow);                        // Stop timer & print
    evalM2M(cells,cells);                                       // Evaluate all M2M kernels
  }

public:
//! Downward phase (M2L,M2P,P2P,L2L,L2P evaluation)
  void downward(Cells &cells, Cells &jcells, bool periodic=true) {
#if HYBRID
    timeKernels();                                              // Time all kernels for auto-tuning
#endif
    for( C_iter C=cells.begin(); C!=cells.end(); ++C ) C->L = 0;// Initialize local coefficients
    if( IMAGES != 0 ) {                                         // If periodic boundary condition
      startTimer("Upward P     ");                              //  Start timer
      upwardPeriodic(jcells);                                   //  Upward phase for periodic images
      stopTimer("Upward P     ",printNow);                      //  Stop timer & print
    }                                                           // Endif for periodic boundary condition
    startTimer("Traverse     ");                                // Start timer
    traverse(cells,jcells);                                     // Traverse tree to get interaction list
    stopTimer("Traverse     ",printNow);                        // Stop timer & print
    if( IMAGES != 0 && periodic ) {                             // If periodic boundary condition
      startTimer("Traverse P   ");                              // Start timer
      traversePeriodic(cells,jcells);                           // Traverse tree for periodic images
      stopTimer("Traverse P   ",printNow);                      // Stop timer & print
    }                                                           // Endif for periodic boundary condition
    evalM2L(cells);                                             // Evaluate queued M2L kernels (only GPU)
    evalM2P(cells);                                             // Evaluate queued M2P kernels (only GPU)
    evalP2P(cells);                                             // Evaluate queued P2P kernels (only GPU)
    evalL2L(cells);                                             // Evaluate all L2L kernels
    evalL2P(cells);                                             // Evaluate all L2P kernels
    if(printNow) std::cout << "P2P: "  << NP2P
                           << " M2P: " << NM2P
                           << " M2L: " << NM2L << std::endl;
  }

//! Calculate Ewald summation
  void Ewald(Bodies &bodies, Cells &cells, Cells &jcells) {
    startTimer("Ewald wave   ");                                // Start timer
    EwaldWave(bodies);                                          // Ewald wave part
    stopTimer("Ewald wave   ",printNow);                        // Stop timer & print
    startTimer("Ewald real   ");                                // Start timer
    neighbor(cells,jcells);                                     // Neighbor calculation for real part
    evalEwaldReal(cells);                                       // Evaluate queued Ewald real kernels (only GPU)
    stopTimer("Ewald real   ",printNow);                        // Stop timer & print
  }
};

#endif
