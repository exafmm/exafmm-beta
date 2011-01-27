#ifndef tree_h
#define tree_h
#include "sort.h"
#include "kernel.h"

class TreeStructure : public Sort {
protected:
  Bodies &bodies;                                               // Bodies in the tree
  Cells  cells;                                                 // Cells in the tree
  C_iter C0;                                                    // Cell iterator begin index
  C_iter CN;                                                    // Cell iterator end index
  C_iter CC0;                                                   // Cell iterator begin index (level-wise)
  C_iter CCN;                                                   // Cell iterator end index (level-wise)
  vect   X0;                                                    // Center of root cell
  real   R0;                                                    // Radius of root cell
  Stack  S;                                                     // Stack of interaction pairs
  Kernel K;                                                     // Kernels
public:
  std::vector<bigint> Ibody;                                    // Cell index of body
  std::vector<bigint> Icell;                                    // Cell index

  TreeStructure(Bodies &b) : bodies(b),X0(0),R0(0) {            // Constructor
    int const N = bodies.size();                                // Number of bodies
    Ibody.resize(N);                                            // Allocate cell index of body
    Icell.resize(N);                                            // Allocate cell index
    cells.resize(N);                                            // Resize cell vector
    for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {        // Loop over all cells
      C->NLEAF = 0;                                             //  Initialize number of leafs per cell
      C->NCHILD = 0;                                            //  Initialize flag for empty child cells
    }                                                           // End loop over all cells
  }

  ~TreeStructure() {}                                           // Destructor

  vect getX0() {return X0;}                                     // Get center of root cell
  real getR0() {return R0;}                                     // Get radius of root cell

  void setDomain() {                                            // Set center and size of root cell
    vect xmin,xmax;                                             // Min,Max of domain
    B_iter B = bodies.begin();                                  // Reset body iterator
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
      R0 = std::max(xmax[d] - X0[d], R0);                       //  Calculate max distance from center
      R0 = std::max(X0[d] - xmin[d], R0);                       //  Calculate max distance from center
    }                                                           // End loop over each dimension
    R0 = pow(2.,int(1. + log(R0) / M_LN2));                     // Add some leeway to root radius
  }

  int getLevel(bigint index) {                                  // Get level from cell index
    int level(-1);                                              // Initialize level counter
    while( index >= 0 ) {                                       // While cell index is non-negative
      level++;                                                  //  Increment level
      index -= 1 << 3*level;                                    //  Subtract number of cells in that level
    }                                                           // End while loop for cell index
    return level;                                               // Return the level
  }

  void getCenter() {                                            // Get cell center and radius from cell index
    int level = getLevel(CN->I);                                // Get level from cell index
    bigint index = CN->I - ((1 << 3*level) - 1) / 7;            // Subtract cell index offset of current level
    CN->R = R0 / (1 << level);                                  // Cell radius
    int d = level = 0;                                          // Initialize dimension and level
    vec<3,int> nx = 0;                                          // Initialize 3-D cell index
    while( index != 0 ) {                                       // Deinterleave bits while index is nonzero
      nx[d] += (index % 2) * (1 << level);                      //  Add deinterleaved bit to 3-D cell index
      index >>= 1;                                              //  Right shift the bits
      d = (d+1) % 3;                                            //  Increment dimension
      if( d == 0 ) level++;                                     //  If dimension is 0 again, increment level
    }                                                           // End while loop for deinterleaving bits
    for( d=0; d!=3; ++d )                                       // Loop over dimensions
      CN->X[d] = (X0[d]-R0) + (2 *nx[d] + 1) * CN->R;           //  Calculate cell center from 3-D cell index
  }

  bigint getParent(bigint index) {                              // Get parent cell index from current cell index
    int level = getLevel(index);                                // Get level from cell index
    bigint cOff = ((1 << 3 *  level   ) - 1) / 7;               // Cell index offset of current level
    bigint pOff = ((1 << 3 * (level-1)) - 1) / 7;               // Cell index offset of parent level
    bigint i = ((index-cOff) >> 3) + pOff;                      // Cell index of parent cell
    return i;                                                   // Return cell index of parent cell
  }

  void sortCells(Cells &buffer) {                               // Sort cells according to cell index
    int begin = CC0-C0;                                         // Begin index for current level
    int end = CCN-C0;                                           // End index for current level
    int c = begin;                                              // Initialize counter for Icell
    for( C_iter C=CC0; C!=CCN; ++C,++c ) Icell[c] = C->I;       // Fill Icell with cell index
    sort(Icell,cells,buffer,false,begin,end);                   // Sort cells according to Icell
  }

  void linkParent(Cells &buffer) {                              // Form parent-child mutual link
    CCN = CN;                                                   // Initialize end iterator for this level
    sortCells(buffer);                                          // Sort cells at this level
    CN->I = getParent(CC0->I);                                  // Set cell index
    CN->LEAF = CC0->LEAF;                                       // Set pointer to first leaf
    getCenter();                                                // Set cell center and radius
    for( C_iter C=CC0; C!=CCN; ++C ) {                          // Loop over all cells at this level
      if( getParent(C->I) != CN->I ) {                          //  If it belongs to a new parent cell
        ++CN;                                                   //   Increment cell iterator
        CN->I = getParent(C->I);                                //   Set cell index
        CN->LEAF = C->LEAF;                                     //   Set pointer to first leaf
        getCenter();                                            //   Set cell center and radius
      }                                                         //  Endif for new parent cell
      for( int i=0; i!=C->NCHILD; ++i )                         //  Loop over child cells
        C->CHILD[i]->PARENT = C;                                //   Link to child to current
      C->PARENT = CN;                                           //  Link to current to parent
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
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B,++b ) {  // Loop over all bodies
      B->pot -= B->scal / std::sqrt(EPS2);                      //  Initialize body values
      if( Ibody[b] != icell ) {                                 //  If it belongs to a new cell
        CN->NLEAF = size;                                       //   Set number of leafs
        CN->NCHILD = 0;                                         //   Set number of child cells
        CN->I = icell;                                          //   Set cell index
        CN->LEAF = begin;                                       //   Set pointer to first leaf
        getCenter();                                            //   Set cell center and radius
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
    CN->I = icell;                                              // Set cell index
    CN->LEAF = begin;                                           // Set pointer to first leaf
    getCenter();                                                // Set cell center and radius
    ++CN;                                                       // Increment cell iterator
    for( int l=level; l>0; --l ) {                              // Once all the twigs are done, do the rest
      linkParent(buffer);                                       //  Form parent-child mutual link
    }                                                           // End upward sweep
  }

  void P2M() {                                                  // Interface for P2M kernel
    for( C_iter C=C0; C!=CN; ++C ) {                            // Loop over all cells
      C->M = 0;                                                 //  Initialize multipole coefficients
      C->L = 0;                                                 //  Initialize local coefficients
      if( C->NLEAF < NCRIT ) {                                  //  If cell is a twig
        K.P2M(C);                                               //   Evaluate P2M kernel
      }                                                         //  Endif for twigs
    }                                                           // End loop over cells
  }

  void M2M() {                                                  // Interface for M2M kernel
    for( C_iter C=C0; C!=CN-1; ++C ) {                          // Loop over all cells bottomup (except root cell)
      K.M2M(C->PARENT,C);                                       //  Evaluate M2M kernel
    }                                                           // End loop over cells
  }

  void M2P(C_iter CI, C_iter CJ) {                              // Interface for M2P kernel
    vect dist = CI->X - CJ->X;                                  // Distance vector between cells
    real R = std::sqrt(norm(dist));                             // Distance between cells
    if( CI->NCHILD != 0 || CI->R + CJ->R > THETA*R ) {          // If target is not twig or box is too large
      pair P(CI,CJ);                                            //  Form pair of interacting cells
      S.push(P);                                                //  Push interacting pair into stack
    } else {                                                    // If target is twig and box is small enough
      K.M2P(CI,CJ);                                             //  Evaluate M2P kernel
    }                                                           // Endif for interaction
  }

  void treecode(C_iter CI, C_iter CJ) {                         // Tree walk for treecode
    if( CI->NCHILD == 0 && CJ->NCHILD == 0) {                   // If both cells are twigs
      K.P2P(CI->LEAF,CI->LEAF+CI->NLEAF,CJ->LEAF,CJ->LEAF+CJ->NLEAF);// Evaluate P2P kernel
    } else if ( CI->NCHILD != 0 ) {                             // If target is not twig
      for( int i=0; i<CI->NCHILD; i++ )                         //  Loop over child cells of target
        M2P(CI->CHILD[i],CJ);                                   //   Try to evaluate M2P kernel
    } else {                                                    // If target is twig
      for( int i=0; i<CJ->NCHILD; i++ )                         //  Loop over child cells of source
        M2P(CI,CJ->CHILD[i]);                                   //   Try to evaluate M2P kernel
    }                                                           // Endif for type of interaction
  }

  void M2L(C_iter CI, C_iter CJ) {                              // Interface for M2L kernel
    vect dist = CI->X - CJ->X;                                  // Distance vector between cells
    real R = std::sqrt(norm(dist));                             // Distance between cells
    if( CI->R + CJ->R > THETA*R ) {                             // If box is too large
      pair P(CI,CJ);                                            //  Form pair of interacting cells
      S.push(P);                                                //  Push interacting pair into stack
    } else {                                                    // If bos is small enough
      K.M2L(CI,CJ);                                             //  Evaluate M2L kernel
    }                                                           // Endif for interaction
  }

  void FMM(C_iter CI, C_iter CJ) {                              // Tree walk for FMM
    if( CI->NCHILD == 0 && CJ->NCHILD == 0 ) {                  // If both cells are twigs
      K.P2P(CI->LEAF,CI->LEAF+CI->NLEAF,CJ->LEAF,CJ->LEAF+CJ->NLEAF);// Evaluate P2P kernel
    } else if ( CJ->NCHILD == 0 || (CI->NCHILD != 0 && CI->R > CJ->R) ) {// If source is twig or target is larger
      for( int i=0; i<CI->NCHILD; i++ )                         //  Loop over child cells of target
        M2L(CI->CHILD[i],CJ);                                   //   Try to evaluate M2L kernel
    } else {                                                    // If target is twig or source is larger
      for( int i=0; i<CJ->NCHILD; i++ )                         //  Loop over child cells of source
        M2L(CI,CJ->CHILD[i]);                                   //   Try to evaluate M2L kernel
    }                                                           // Endif for type of interaction
  }

  void evaluate(int method) {                                   // Interface for treewalk
    pair P0(CN-1,CN-1);                                         // Form pair of root cells
    S.push(P0);                                                 // Push pair to stack
    while( !S.empty() ) {                                       // While interaction stack is not empty
      pair P = S.top();                                         //  Get interaction pair from top of stack
      S.pop();                                                  //  Pop interaction stack
      switch (method) {                                         //  Swtich between methods
        case 0 : treecode(P.CI,P.CJ); break;                    //   0 : treecode
        case 1 : FMM(P.CI,P.CJ);      break;                    //   1 : FMM
      }                                                         //  End switch between methods
    }                                                           // End while loop for interaction stack
    for( C_iter C=CN-2; C!=C0-1; --C ) {                        // Loop over all cells topdown (except root cell)
      K.L2L(C,C->PARENT);                                       //  Evaluate L2L kernel
      if( C->NLEAF < NCRIT )                                    //  If cell is a twig
        K.L2P(C);                                               //   Evaluate L2P kernel
    }
  }

};

#endif
