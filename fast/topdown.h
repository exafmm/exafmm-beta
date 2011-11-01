#ifndef topdown_h
#define topdown_h
#include "evaluator.h"

class TopDown : public Evaluator {
private:
  Nodes    nodes;
  Leafs    leafs;
  B_iter   BN;
  unsigned NLEAF;

protected:
  int      MAXLEVEL;
  Cells    cells;

private:
  void init(Node &node) {
    node.NOCHILD = true;
    node.NLEAF = 0;
    node.LEAF = NULL;
    for( int b=0; b!=8; ++b ) node.CHILD[b] = -1;
  }

  inline void addChild(int octant, N_iter N) {
    Node child;
    init(child);
    child.LEVEL = N->LEVEL+1;
    child.X = N->X;
    real r = R0 / (1 << child.LEVEL);
    for( int d=0; d!=3; ++d ) {                                 // Loop over dimensions
      child.X[d] += r * (((octant & 1 << d) >> d) * 2 - 1);     //  Calculate new center position
    }                                                           // End loop over dimensions
    N->NOCHILD = false;
    N->CHILD[octant] = nodes.size();
    nodes.push_back(child);
    NCELL++;
  }

  void splitNode(N_iter N) {
    while( N->NLEAF > NCRIT ) {
      int c;
      Leaf *Ln;
      for( Leaf *L=N->LEAF; L; L=Ln ) {
        Ln = L->NEXT;
        int octant = (L->X[0] > N->X[0]) + ((L->X[1] > N->X[1]) << 1) + ((L->X[2] > N->X[2]) << 2);
        if( N->CHILD[octant] == -1 ) {
          addChild(octant,N);
        }
        c = N->CHILD[octant];
        Node *child = &nodes[c];
        L->NEXT = child->LEAF;
        child->LEAF = &*L;
        child->NLEAF++;
      }
      N = nodes.begin()+c;
    }
  }

  void nodes2cells(int i, C_iter C) {
    C->R      = R0 / (1 << nodes[i].LEVEL);
    C->X      = nodes[i].X;
    C->NDLEAF = nodes[i].NLEAF;
    C->LEAF   = BN;
    if( nodes[i].NOCHILD ) {
      C->CHILD = 0;
      C->NCHILD = 0;
      C->NCLEAF = nodes[i].NLEAF;
      for( Leaf *L=nodes[i].LEAF; L; L=L->NEXT ) {
        BN->IBODY = L->I;
        BN++;
      }
    } else {
      C->NCLEAF = 0;
      int nsub=0;
      for( int octant=0; octant!=8; ++octant ) {
        if( nodes[i].CHILD[octant] != -1 ) {
          ++nsub;
        }
      }
      C_iter Ci = CN;
      C->CHILD = Ci - C0;
      C->NCHILD = nsub;
      CN += nsub;
      for( int octant=0; octant!=8; ++octant ) {
        if( nodes[i].CHILD[octant] != -1 ) {
          Ci->PARENT = C - C0;
          nodes2cells(nodes[i].CHILD[octant], Ci++);
        }
      }
    }
  }

protected:
  void sortBodies2(Bodies &bodies) {
    Bodies buffer = bodies;
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      B->X   = buffer[B->IBODY].X;
      B->SRC = buffer[B->IBODY].SRC;
      B->TRG = buffer[B->IBODY].TRG;
    }
  }

  void setDomain(Bodies &bodies) {
    vect xmin, xmax;
    NLEAF = bodies.size();
    leafs.reserve(NLEAF);
    X0 = 0;
    xmax = xmin = bodies.begin()->X;
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      Leaf leaf;                                                //  Leafs are bodies attached to the tree
      leaf.I = B-bodies.begin();                                //  Set leaf index
      leaf.X = B->X;                                            //  Set leaf coordinate
      leafs.push_back(leaf);                                    //  Push leaf into vector
      for( int d=0; d!=3; ++d ) {                               //  Loop over each dimension
        if     (B->X[d] < xmin[d]) xmin[d] = B->X[d];           //   Determine xmin
        else if(B->X[d] > xmax[d]) xmax[d] = B->X[d];           //   Determine xmax
      }                                                         //  End loop over each dimension
      X0 += B->X;                                               //  Sum positions
    }                                                           // End loop over bodies
    X0 /= bodies.size();                                        // Calculate average position
    for( int d=0; d!=3; ++d ) {                                 // Loop over each dimension
      X0[d] = int(X0[d]+.5);                                    //  Shift center to nearest integer
      R0 = std::max(xmax[d] - X0[d], R0);                       //  Calculate max distance from center
      R0 = std::max(X0[d] - xmin[d], R0);                       //  Calculate max distance from center
    }                                                           // End loop over each dimension
    R0 *= 1.000001;                                             // Add some leeway to root radius
  }

  void build() {
    startTimer("Grow tree    ");
    NCELL = 1;
    nodes.reserve(NLEAF);
    Node node;
    init(node);
    node.LEVEL = 0;
    node.X     = X0;
    nodes.push_back(node);
    MAXLEVEL = 0;
    for( L_iter L=leafs.begin(); L!=leafs.end(); ++L ) {
      int i = 0;
      N_iter N = nodes.begin()+i;
      while( !N->NOCHILD ) {
        int octant = (L->X[0] > N->X[0]) + ((L->X[1] > N->X[1]) << 1) + ((L->X[2] > N->X[2]) << 2);
        N->NLEAF++;
        if( N->CHILD[octant] == -1 ) addChild(octant,N);
        N = nodes.begin()+N->CHILD[octant];
      }
      L->NEXT = N->LEAF;
      N->LEAF = &*L;
      N->NLEAF++;
      if( N->NLEAF > NCRIT ) splitNode(N);
      if( MAXLEVEL < N->LEVEL ) MAXLEVEL = N->LEVEL;
    }
    MAXLEVEL++;
    stopTimer("Grow tree    ",printNow);
  }

  void link(Bodies &bodies) {
    startTimer("Link tree    ");
    cells.resize(NCELL);
    C0 = cells.begin();
    BN = bodies.begin();
    CN = C0+1;
    nodes2cells(0,C0);
    CN = C0;
    nodes.clear();
    leafs.clear();
    sortBodies2(bodies);
    stopTimer("Link tree    ",printNow);
  }

public:
  TopDown() {}
  ~TopDown() {}
};

#endif
