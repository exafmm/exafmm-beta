#ifndef tree_h
#define tree_h
#include "sort.h"
#include "kernel.h"

class TreeStructure : virtual public Sort {
protected:
  Bodies &bodies;                                               // Bodies in the tree
  Cells  cells;                                                 // Cells in the tree
  Pairs  pairs;                                                 // Stack of interaction pairs
  C_iter C0;                                                    // cells.begin()
  vect   X0;                                                    // Center of root cell
  real   R0;                                                    // Radius of root cell
  Kernel K;                                                     // Kernels
public:
  std::vector<bigint> Ibody;                                    // Cell index of body
  std::vector<bigint> Icell;                                    // Cell index

  TreeStructure(Bodies &b) : bodies(b),X0(0),R0(0) {            // Constructor
    Ibody.resize(bodies.size());                                // Resize cell index of body
  }

  ~TreeStructure() {}                                           // Destructor

  vect &setX0() {return X0;}                                    // Get center of root cell
  real &setR0() {return R0;}                                    // Get radius of root cell
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
    int level = -1;                                             // Initialize level counter
    while( index >= 0 ) {                                       // While cell index is non-negative
      level++;                                                  //  Increment level
      index -= 1 << 3*level;                                    //  Subtract number of cells in that level
    }                                                           // End while loop for cell index
    return level;                                               // Return the level
  }

  void getCenter(Cell &cell) {                                  // Get cell center and radius from cell index
    int level = getLevel(cell.I);                               // Get level from cell index
    bigint index = cell.I - ((1 << 3*level) - 1) / 7;           // Subtract cell index offset of current level
    cell.R = R0 / (1 << level);                                 // Cell radius
    int d = level = 0;                                          // Initialize dimension and level
    vec<3,int> nx = 0;                                          // Initialize 3-D cell index
    while( index != 0 ) {                                       // Deinterleave bits while index is nonzero
      nx[d] += (index % 2) * (1 << level);                      //  Add deinterleaved bit to 3-D cell index
      index >>= 1;                                              //  Right shift the bits
      d = (d+1) % 3;                                            //  Increment dimension
      if( d == 0 ) level++;                                     //  If dimension is 0 again, increment level
    }                                                           // End while loop for deinterleaving bits
    for( d=0; d!=3; ++d )                                       // Loop over dimensions
      cell.X[d] = (X0[d]-R0) + (2 *nx[d] + 1) * cell.R;         //  Calculate cell center from 3-D cell index
  }

  bigint getParent(bigint index) {                              // Get parent cell index from current cell index
    int level = getLevel(index);                                // Get level from cell index
    bigint cOff = ((1 << 3 *  level   ) - 1) / 7;               // Cell index offset of current level
    bigint pOff = ((1 << 3 * (level-1)) - 1) / 7;               // Cell index offset of parent level
    bigint i = ((index-cOff) >> 3) + pOff;                      // Cell index of parent cell
    return i;                                                   // Return cell index of parent cell
  }

  void sortCells(Cells &buffer, int begin, int end) {           // Sort cells according to cell index
    Icell.resize(cells.size());                                 // Resize vector for cell index
    buffer.resize(cells.size());                                // Resize vector for sort buffer
    for( int i=begin; i!=end; ++i ) Icell[i] = cells[i].I;      // Fill Icell with cell index
    sort(Icell,cells,buffer,false,begin,end);                   // Sort cells according to Icell
  }

  void linkParent(int &begin, int &end) {                       // Form parent-child mutual link
    Cell parent;                                                // Define parent cell structure
    int oldend = end;                                           // Save old end counter
    parent.I = getParent(cells[begin].I);                       // Set cell index
    parent.NLEAF = parent.NCHILD = 0;                           // Initialize NLEAF & NCHILD
    parent.LEAF = cells[begin].LEAF;                            // Set pointer to first leaf
    getCenter(parent);                                          // Set cell center and radius
    for( int i=begin; i!=oldend; ++i ) {                        // Loop over all cells at this level
      if( getParent(cells[i].I) != parent.I ) {                 //  If it belongs to a new parent cell
        cells.push_back(parent);                                //   Push cell structure into vector
        end++;                                                  //   Increment cell counter
        parent.I = getParent(cells[i].I);                       //   Set cell index
        parent.NLEAF = parent.NCHILD = 0;                       //   Initialize NLEAF & NCHILD
        parent.LEAF = cells[i].LEAF;                            //   Set pointer to first leaf
        getCenter(parent);                                      //   Set cell center and radius
      }                                                         //  Endif for new parent cell
      for( int c=0; c!=cells[i].NCHILD; ++c )                   //  Loop over child cells
        (cells.begin()+cells[i].CHILD[c])->PARENT = i;          //   Link child to current
      cells[i].PARENT = end;                                    //  Link to current to parent
      parent.NLEAF += cells[i].NLEAF;                           //  Add nleaf of child to parent
      parent.M     += cells[i].M;                               //  Add multipoles of child to parent
      parent.CHILD[parent.NCHILD] = i;                          //  Link to child
      parent.NCHILD++;                                          //  Increment child counter
    }                                                           // End loop over all cells at this level
    cells.push_back(parent);                                    // Push cell structure into vector
    end++;                                                      // Increment cell counter
    begin = oldend;
  }

  void link() {                                                 // Link cells to create tree
    int begin=0, end=0;                                         // Initialize range of cell vector
    int index=Ibody[0], nleaf=0, level=getLevel(index);         // Initialize cell index, nleaf, level
    B_iter firstLeaf = bodies.begin();                          // Initialize body iterator
    Cell cell;                                                  // Define cell structure
    Cells buffer;                                               // Allocate sort buffer for cells
    BI_iter BI = Ibody.begin();                                 // Initialize body index iterator
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B,++BI ) { // Loop over all bodies
      if( *BI != index ) {                                      //  If it belongs to a new cell
        cell.NLEAF = nleaf;                                     //   Set number of leafs
        cell.NCHILD = 0;                                        //   Set number of child cells
        cell.I = index;                                         //   Set cell index
        cell.LEAF = firstLeaf;                                  //   Set pointer to first leaf
        getCenter(cell);                                        //   Set cell center and radius
        cells.push_back(cell);                                  //   Push cell structure into vector
        end++;                                                  //   Increment cell counter
        while( getLevel(*BI) != level ) {                       //   While cell belongs to a higher level
          sortCells(buffer,begin,end);                          //    Sort cells at this level
          linkParent(begin,end);                                //    Form parent-child mutual link
          level--;                                              //    Go up one level
        }                                                       //   Endif for new level
        firstLeaf = B;                                          //   Set new first leaf
        nleaf = 0;                                              //   Reset number of bodies
        index = *BI;                                            //   Set new cell
      }                                                         //  Endif for new cell
      nleaf++;                                                  //  Increment body counter
    }                                                           // End loop over bodies
    cell.NLEAF = nleaf;                                         // Set number of leafs
    cell.NCHILD = 0;                                            // Set number of child cells
    cell.I = index;                                             // Set cell index
    cell.LEAF = firstLeaf;                                      // Set pointer to first leaf
    getCenter(cell);                                            // Set cell center and radius
    cells.push_back(cell);                                      // Push cell structure into vector
    end++;                                                      // Increment cell counter
    for( int l=level; l>0; --l ) {                              // Once all the twigs are done, do the rest
      sortCells(buffer,begin,end);                              //  Sort cells at this level
      linkParent(begin,end);                                    //  Form parent-child mutual link
    }
    C0 = cells.begin();
  }

  void P2M() {                                                  // Interface for P2M kernel
    for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {        // Loop over all cells
      C->M = 0;                                                 //  Initialize multipole coefficients
      C->L = 0;                                                 //  Initialize local coefficients
      if( C->NLEAF < NCRIT ) {                                  //  If cell is a twig
        K.P2M(C);                                               //   Evaluate P2M kernel
      }                                                         //  Endif for twigs
    }                                                           // End loop over cells
  }

  void M2M() {                                                  // Interface for M2M kernel
    for( C_iter C=cells.begin(); C!=cells.end()-1; ++C ) {      // Loop over all cells bottomup (except root cell)
      K.M2M(C0+C->PARENT,C);                                    //  Evaluate M2M kernel
    }                                                           // End loop over cells
  }

  void M2P(C_iter CI, C_iter CJ) {                              // Interface for M2P kernel
    vect dist = CI->X - CJ->X;                                  // Distance vector between cells
    real R = std::sqrt(norm(dist));                             // Distance between cells
    if( CI->NCHILD != 0 || CI->R + CJ->R > THETA*R ) {          // If target is not twig or box is too large
      Pair pair(CI,CJ);                                         //  Form pair of interacting cells
      pairs.push(pair);                                         //  Push interacting pair into stack
    } else {                                                    // If target is twig and box is small enough
      K.M2P(CI,CJ);                                             //  Evaluate M2P kernel
    }                                                           // Endif for interaction
  }

  void treecode(C_iter CI, C_iter CJ) {                         // Tree walk for treecode
    if( CI->NCHILD == 0 && CJ->NCHILD == 0) {                   // If both cells are twigs
      K.P2P(CI->LEAF,CI->LEAF+CI->NLEAF,CJ->LEAF,CJ->LEAF+CJ->NLEAF);// Evaluate P2P kernel
    } else if ( CI->NCHILD != 0 ) {                             // If target is not twig
      for( int i=0; i<CI->NCHILD; i++ )                         //  Loop over child cells of target
        M2P(C0+CI->CHILD[i],CJ);                                //   Try to evaluate M2P kernel
    } else {                                                    // If target is twig
      for( int i=0; i<CJ->NCHILD; i++ )                         //  Loop over child cells of source
        M2P(CI,C0+CJ->CHILD[i]);                                //   Try to evaluate M2P kernel
    }                                                           // Endif for type of interaction
  }

  void M2L(C_iter CI, C_iter CJ) {                              // Interface for M2L kernel
    vect dist = CI->X - CJ->X;                                  // Distance vector between cells
    real R = std::sqrt(norm(dist));                             // Distance between cells
    if( CI->R + CJ->R > THETA*R ) {                             // If box is too large
      Pair pair(CI,CJ);                                         //  Form pair of interacting cells
      pairs.push(pair);                                         //  Push interacting pair into stack
    } else {                                                    // If bos is small enough
      K.M2L(CI,CJ);                                             //  Evaluate M2L kernel
    }                                                           // Endif for interaction
  }

  void FMM(C_iter CI, C_iter CJ) {                              // Tree walk for FMM
    if( CI->NCHILD == 0 && CJ->NCHILD == 0 ) {                  // If both cells are twigs
      K.P2P(CI->LEAF,CI->LEAF+CI->NLEAF,CJ->LEAF,CJ->LEAF+CJ->NLEAF);// Evaluate P2P kernel
    } else if ( CJ->NCHILD == 0 || (CI->NCHILD != 0 && CI->R > CJ->R) ) {// If source is twig or target is larger
      for( int i=0; i<CI->NCHILD; i++ )                         //  Loop over child cells of target
        M2L(C0+CI->CHILD[i],CJ);                                //   Try to evaluate M2L kernel
    } else {                                                    // If target is twig or source is larger
      for( int i=0; i<CJ->NCHILD; i++ )                         //  Loop over child cells of source
        M2L(CI,C0+CJ->CHILD[i]);                                //   Try to evaluate M2L kernel
    }                                                           // Endif for type of interaction
  }

  void evaluate(int method) {                                   // Interface for treewalk
    C_iter root = cells.end()-1;                                // Iterator for root cell
    Pair pair(root,root);                                       // Form pair of root cells
    pairs.push(pair);                                           // Push pair to stack
    while( !pairs.empty() ) {                                   // While interaction stack is not empty
      pair = pairs.top();                                       //  Get interaction pair from top of stack
      pairs.pop();                                              //  Pop interaction stack
      switch (method) {                                         //  Swtich between methods
        case 0 : treecode(pair.CI,pair.CJ); break;              //   0 : treecode
        case 1 : FMM(pair.CI,pair.CJ);      break;              //   1 : FMM
      }                                                         //  End switch between methods
    }                                                           // End while loop for interaction stack
    for( C_iter C=root-1; C!=cells.begin()-1; --C ) {           // Loop over all cells topdown (except root cell)
      K.L2L(C,C0+C->PARENT);                                    //  Evaluate L2L kernel
      if( C->NLEAF < NCRIT )                                    //  If cell is a twig
        K.L2P(C);                                               //   Evaluate L2P kernel
    }
  }

};

#endif
