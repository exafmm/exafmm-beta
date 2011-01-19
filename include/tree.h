#ifndef tree_h
#define tree_h
#include "types.h"
#include "sort.h"

class TreeStructure : public Sort {
protected:
  Bodies &bodies;                                               // Bodies in the tree
  B_iter B;                                                     // Body iterator for target = source
  B_iter Bi;                                                    // Body iterator for target
  B_iter Bj;                                                    // Body iterator for source
  Cells  cells;                                                 // Cells in the tree
  C_iter C;                                                     // Cell iterator
  C_iter C0;                                                    // Cell iterator begin index
  C_iter CN;                                                    // Cell iterator end index
  C_iter CC0;                                                   // Cell iterator begin index (level-wise)
  C_iter CCN;                                                   // Cell iterator end index (level-wise)
  vect X0;                                                      // Center of root cell
  real R0;                                                      // Radius of root cell
public:
  bigint *Ibody;                                                // Morton index of body
  bigint *Icell;                                                // Morton index of cell

  TreeStructure(Bodies &b) : bodies(b),X0(0),R0(0) {            // Constructor
    int const N = bodies.size();                                // Number of bodies
    Ibody = new bigint [N];                                     // Allocate Morton index of body
    Icell = new bigint [N];                                     // Allocate Morton index of cell
    cells.resize(N);                                            // Resize cell vector
    for( C=cells.begin(); C!=cells.end(); ++C ) {               // Loop over all cells
      C->NLEAF = 0;                                             //  Initialize number of leafs per cell
      C->NCHILD = 0;                                            //  Initialize flag for empty child cells
    }                                                           // End loop over all cells
  }

  ~TreeStructure() {                                            // Destructor
    delete[] Ibody;                                             // Delete Morton index of body
    delete[] Icell;                                             // Delete Morton index of body
  }

  vect getX0() {return X0;}                                     // Get center of root cell
  real getR0() {return R0;}                                     // Get radius of root cell

  void setDomain() {                                            // Set center and size of root cell
    vect xmin,xmax;                                             // Min,Max of domain
    B = bodies.begin();                                         // Reset body iterator
    xmin = xmax = B->pos;                                       // Initialize xmin,xmax
    for( B=bodies.begin(); B!=bodies.end(); ++B ) {             // Loop over all bodies
      for( int d=0; d!=3; ++d ) {                               //  Loop over each dimension
        if     (B->pos[d] < xmin[d]) xmin[d] = B->pos[d];       //   Determine xmin
        else if(B->pos[d] > xmax[d]) xmax[d] = B->pos[d];       //   Determine xmax
      }                                                         //  End loop over each dimension
      X0 += B->pos;                                             //  Sum positions
    }                                                           // End loop over all bodies
    X0 /= bodies.size();                                        // Calculate average position
    for( int d=0; d!=3; ++d ) {                                 // Loop over each dimension
      X0[d] = int(X0[d]+.5);                                    //  Shift center to nearest integer
      R0 = std::max(xmax[d]-X0[d],R0);                          //  Calculate max distance from center
      R0 = std::max(X0[d]-xmin[d],R0);                          //  Calculate max distance from center
    }                                                           // End loop over each dimension
    R0 = pow(2.,int(1.+log(R0)/M_LN2));                         // Add some leeway to root radius
  }

  int getLevel(bigint index) {                                  // Get level from Morton index
    int level(-1);                                              // Initialize level counter
    while( index >= 0 ) {                                       // While Morton index is non-negative
      level++;                                                  //  Increment level
      index -= 1 << 3*level;                                    //  Subtract number of cells in that level
    }                                                           // End while loop for Morton index
    return level;                                               // Return the level
  }

  bigint getParent(bigint index) {                              // Get parent Morton index from current index
    int level = getLevel(index);                                // Get level from Morton index
    int cOffset = ((1 << 3*level) - 1)/7;                       // Morton offset of current level
    int pOffset = ((1 << 3*(level-1)) - 1)/7;                   // Morton offset of parent level
    bigint i = ((index-cOffset) >> 3) + pOffset;                // Morton index of parent cell
    return i;                                                   // Return Morton index of parent cell
  }

  void sortCells(Cells &buffer) {                               // Sort cells according to Morton index
    int begin = CC0-C0;                                         // Begin index for current level
    int end = CCN-C0;                                           // End index for current level
    int c = begin;                                              // Initialize counter for Icell
    for( C=CC0; C!=CCN; ++C,++c ) Icell[c] = C->I;              // Fill Icell with Morton index
    sort(Icell,cells,buffer,false,begin,end);                   // Sort cells according to Icell
  }

  void linkParent(Cells &buffer) {                              // Form parent-child mutual link
    CCN = CN;                                                   // Initialize end iterator for this level
    sortCells(buffer);                                          // Sort cells at this level
    CN->I = getParent(CC0->I);                                  // Assign parent cell index
    CN->LEAF = CC0->LEAF;                                       // Assign pointer to first leaf
    for( C=CC0; C!=CCN; ++C ) {                                 // Loop over all cells at this level
      if( getParent(C->I) != CN->I ) {                          //  If it belongs to a new parent cell
        ++CN;                                                   //   Increment cell iterator
        CN->I = getParent(C->I);                                //   Assign parent cell index
        CN->LEAF = C->LEAF;                                     //   Assign pointer to first leaf
      }                                                         //  Endif for new parent cell
      C->PARENT = CN;                                           //  Link to parent
      CN->NLEAF += C->NLEAF;                                    //  Add nleaf of child to parent
      CN->CHILD[CN->NCHILD] = C;                                //  Link to child
      CN->NCHILD++;                                             //  Increment child counter
    }                                                           // End loop over all cells at this level
    ++CN;                                                       // Increment cell iterator
    CC0 = CCN;                                                  // Set new begin iterator for next level
  }

  void link() {                                                 // Link cells to create tree
    CCN = CC0 = CN = C0 = cells.begin();                        // Initialize cell iterators
    int icell(Ibody[0]),size(0),level(getLevel(icell));         // Initialize cell index, size, level
    B_iter begin(bodies.begin());                               // Initialize body iterator
    Cells buffer(cells.size());                                 // Allocate sort buffer for cells
    int b = 0;                                                  // Initialize body counter
    for( B=bodies.begin(); B!=bodies.end(); ++B,++b ) {         // Loop over all bodies
      if( Ibody[b] != icell ) {                                 //  If it belongs to a new cell
        CN->NLEAF = size;                                       //   Set number of leafs
        CN->NCHILD = 0;                                         //   Set number of child cells
        CN->I = icell;                                          //   Set Morton index
        CN->LEAF = begin;                                       //   Set pointer to first leaf
        ++CN;                                                   //   Increment cell iterator
        if( getLevel(Ibody[b]) != level ) {                     //   If cell belongs to a new level
          linkParent(buffer);                                   //    Form parent-child mutual link
          level = getLevel(Ibody[b]);                           //    Set new level
        }                                                       //   Endif for new level
        begin = B;                                              //   Set new begin iterator
        size = 0;                                               //   Reset number of bodies
        icell = Ibody[b];                                       //   Set new cell
      }                                                         //  Endif for new cell
      size++;                                                   //  Increment body counter
    }                                                           // End loop over bodies
    CN->NLEAF = size;                                           // Set number of leafs
    CN->NCHILD = 0;                                             // Set number of child cells
    CN->I = icell;                                              // Set Morton index
    CN->LEAF = begin;                                           // Set pointer to first leaf
    ++CN;                                                       // Increment cell iterator
    for( int l=level; l>0; --l ) {                              // Once all the twigs are done, do the rest
      linkParent(buffer);                                       //  Form parent-child mutual link
    }                                                           // End upward sweep

    for( C=C0; C!=CN; ++C ) {
//      std::cout << C->I << " " << C->NLEAF << std::endl;
    }
  }

};

#endif
