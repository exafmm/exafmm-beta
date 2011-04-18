#ifndef tree_h
#define tree_h
#include "sort.h"
#include "evaluator.h"

class TreeStructure : public Evaluator, public Sort {           // Base class for tree structure
public:
  Bodies buffer;                                                // Buffer for comm & sort

private:
  bigint getParent(bigint index) {                              // Get parent cell index from current cell index
    int level = getLevel(index);                                // Get level from cell index
    bigint cOff = ((1 << 3 *  level   ) - 1) / 7;               // Cell index offset of current level
    bigint pOff = ((1 << 3 * (level-1)) - 1) / 7;               // Cell index offset of parent level
    bigint i = ((index-cOff) >> 3) + pOff;                      // Cell index of parent cell
    return i;                                                   // Return cell index of parent cell
  }

  void unique(Cells &cells, Cells &sticks, int begin, int &end) {// Merge sticks with cells (levelwise)
    int c_old = begin;                                          // Initialize old cell counter
    for( int c=begin; c!=end; ++c ) {                           // Loop over cells in level
      if( cells[c].ICELL != cells[c_old].ICELL ) {              //  If current cell index is different from previous
        c_old = c;                                              //   Update old cell counter
      } else if( c != c_old ) {                                 //  If cell index is repeated
        if( cells[c].NCHILD != 0 ) {                            //   Stick-cell collision
          cells[c_old].NCHILD = cells[c].NCHILD;                //    Copy number of children
          cells[c_old].NLEAF  = cells[c].NLEAF;                 //    Copy number of leafs
          cells[c_old].PARENT = cells[c].PARENT;                //    Copy parent link
          for( int i=0; i!=cells[c].NCHILD; ++i ) {             //    Loop over children
            cells[c_old].CHILD[i] = cells[c].CHILD[i];          //     Copy child link
          }                                                     //    End loop over children
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

  void linkParent(Cells &cells, int &begin, int &end) {         // Form parent-child mutual link
    Cell parent;                                                // Parent cell
    Cells parents;                                              // Parent cell vector;
    int oldend = end;                                           // Save old end counter
    parent.ICELL = getParent(cells[begin].ICELL);               // Set cell index
    parent.M = parent.L = 0;                                    // Initlalize multipole & local coefficients
    parent.NLEAF = parent.NCHILD = 0;                           // Initialize NLEAF & NCHILD
    parent.LEAF = cells[begin].LEAF;                            // Set pointer to first leaf
    getCenter(parent);                                          // Set cell center and radius
    for( int i=begin; i!=oldend; ++i ) {                        // Loop over cells at this level
      if( getParent(cells[i].ICELL) != parent.ICELL ) {         //  If it belongs to a new parent cell
        cells.push_back(parent);                                //   Push cells into vector
        end++;                                                  //   Increment cell counter
        parent.ICELL = getParent(cells[i].ICELL);               //   Set cell index
        parent.M = parent.L = 0;                                //   Initialize multipole & local coefficients
        parent.NLEAF = parent.NCHILD = 0;                       //   Initialize NLEAF & NCHILD
        parent.LEAF = cells[i].LEAF;                            //   Set pointer to first leaf
        getCenter(parent);                                      //   Set cell center and radius
      }                                                         //  Endif for new parent cell
      for( int c=0; c!=cells[i].NCHILD; ++c ) {                 //  Loop over child cells
        (cells.begin()+cells[i].CHILD[c])->PARENT = i;          //   Link child to current
      }                                                         //  End loop over child cells
      cells[i].PARENT = end;                                    //  Link to current to parent
      parent.NLEAF += cells[i].NLEAF;                           //  Add nleaf of child to parent
      parents.push_back(parent);                                //  Push parent cell into vector
      parent = parents.back();                                  //  Copy information from vector
      parents.pop_back();                                       //  Pop parent cell from vector
      parent.CHILD[parent.NCHILD] = i;                          //  Link to child
      parent.NCHILD++;                                          //  Increment child counter
    }                                                           // End loop over cells at this level
    cells.push_back(parent);                                    // Push cells into vector
    end++;                                                      // Increment cell counter
    begin = oldend;                                             // Set new begin index to old end index
  }

protected:
  int getLevel(bigint index) {                                  // Get level from cell index
    int level = -1;                                             // Initialize level counter
    while( index >= 0 ) {                                       // While cell index is non-negative
      level++;                                                  //  Increment level
      index -= 1 << 3*level;                                    //  Subtract number of cells in that level
    }                                                           // End while loop for cell index
    return level;                                               // Return the level
  }

  void getCenter(Cell &cell) {                                  // Get cell center and radius from cell index
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

  void bodies2twigs(Bodies &bodies, Cells &twigs){              // Group bodies into twig cells
    startTimer("Bodies2twigs ");                                // Start timer
    int nleaf = 0;                                              // Initialize number of leafs
    bigint index = bodies[0].ICELL;                             // Initialize cell index
    B_iter firstLeaf = bodies.begin();                          // Initialize body iterator for first leaf
    Cell cell;                                                  // Cell structure
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      if( B->ICELL != index ) {                                 //  If it belongs to a new cell
        cell.NLEAF  = nleaf;                                    //   Set number of leafs
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
    cell.NLEAF  = nleaf;                                        // Set number of leafs
    cell.NCHILD = 0;                                            // Set number of child cells
    cell.ICELL  = index;                                        // Set cell index
    cell.LEAF   = firstLeaf;                                    // Set pointer to first leaf
    getCenter(cell);                                            // Set cell center and radius
    twigs.push_back(cell);                                      // Push cells into vector
    stopTimer("Bodies2twigs ",printNow);                        // Stop timer & print
    evalP2M(twigs);                                             // Evaluate P2M kernel
  }

  void twigs2cells(Cells &twigs, Cells &cells, Cells &sticks) { // Link twigs bottomup to create all cells in tree
    startTimer("Twigs2cells  ");                                // Start timer
    int begin = 0, end = 0;                                     // Initialize range of cell vector
    int level = getLevel(twigs.back().ICELL);                   // Initialize level of tree
    while( !twigs.empty() ) {                                   // Keep poppig twigs until the vector is empty
      while( getLevel(twigs.back().ICELL) != level ) {          //  While cell belongs to a higher level
        sortCells(cells,false,begin,end);                       //   Sort cells at this level
        unique(cells,sticks,begin,end);                         //   Get rid of duplicate cells
        linkParent(cells,begin,end);                            //   Form parent-child mutual link
        level--;                                                //   Go up one level
      }                                                         //  End while for higher level
      cells.push_back(twigs.back());                            // Push cells into vector
      twigs.pop_back();                                         // Pop twigs from vector
      end++;                                                    // Increment cell counter
    }                                                           // End while for popping twigs
    for( int l=level; l>0; --l ) {                              // Once all the twigs are done, do the rest
      sortCells(cells,false,begin,end);                         //  Sort cells at this level
      unique(cells,sticks,begin,end);                           //  Get rid of duplicate cells
      linkParent(cells,begin,end);                              //  Form parent-child mutual link
    }                                                           // End loop over levels
    unique(cells,sticks,begin,end);                             // Just in case there is a collision at root
    stopTimer("Twigs2cells  ",printNow);                        // Stop timer & print
    evalM2M(cells);                                             // Evaluate M2M kernel
  }

public:
  void downward(Cells &cells, Cells &jcells, int method, bool periodic=true) {// Downward phase
    for( C_iter C=cells.begin(); C!=cells.end(); ++C ) C->L = 0;// Initialize local coefficients
    if( IMAGES != 0 ) {                                         // If periodic boundary condition
      startTimer("Upward P     ");                              //  Start timer
      upwardPeriodic(jcells);                                   //  Upward phase for periodic images
      stopTimer("Upward P     ",printNow);                      //  Stop timer & print
    }                                                           // Endif for periodic boundary condition
    startTimer("Traverse     ");                                // Start timer
    traverse(cells,jcells,method);                              // Traverse tree to get interaction list
    stopTimer("Traverse     ",printNow);                        // Stop timer & print
    if( IMAGES != 0 && periodic ) {                             // If periodic boundary condition
      startTimer("Traverse P   ");                              // Start timer
      traversePeriodic(cells,jcells,method);                    // Traverse tree for periodic images
      stopTimer("Traverse P   ",printNow);                      // Stop timer & print
    }                                                           // Endif for periodic boundary condition
    evalM2L(cells);                                             // Evaluate M2L kernel
    evalM2P(cells);                                             // Evaluate M2P kernel
    evalP2P(cells);                                             // Evaluate P2P kernel
    evalL2L(cells);                                             // Evaluate L2L kernel
    evalL2P(cells);                                             // Evaluate L2P kernel
  }
};

#endif
