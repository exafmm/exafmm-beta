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
#ifndef topdown_h
#define topdown_h
#include "evaluator.h"

class TopDown : public Evaluator {
private:
  Nodes    nodes;
  Leafs    leafs;
  B_iter   BN;
  C_iter   CN;
  unsigned NLEAF;
  unsigned NCELL;

protected:
  int      MAXLEVEL;

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
      int c = 0;
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
      C->CHILD = Ci - Ci0;
      C->NCHILD = nsub;
      CN += nsub;
      for( int octant=0; octant!=8; ++octant ) {
        if( nodes[i].CHILD[octant] != -1 ) {
          Ci->PARENT = C - Ci0;
          nodes2cells(nodes[i].CHILD[octant], Ci++);
        }
      }
    }
  }

  void permuteBodies(Bodies &bodies) {
    Bodies buffer = bodies;
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      B->ICELL = buffer[B->IBODY].ICELL;
      B->X     = buffer[B->IBODY].X;
      B->SRC   = buffer[B->IBODY].SRC;
      B->TRG   = buffer[B->IBODY].TRG;
    }
  }

protected:
  void setDomain(Bodies &bodies) {
    startTimer("Set domain   ");
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
    stopTimer("Set domain   ",printNow);
  }

  void buildTree() {
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

  void linkTree(Bodies &bodies, Cells &cells) {
    startTimer("Link tree    ");
    cells.resize(NCELL);
    Ci0 = cells.begin();
    BN = bodies.begin();
    CN = Ci0 + 1;
    nodes2cells(0,Ci0);
    nodes.clear();
    leafs.clear();
    permuteBodies(bodies);
    stopTimer("Link tree    ",printNow);
  }

  void upwardPass(Cells &cells) {
    startTimer("Upward pass  ");
    setRootCell(cells);
    for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {
      C->M = 0;
      C->L = 0;
    }
    for( C_iter C=cells.end()-1; C!=cells.begin()-1; --C ) {
      real Rmax = 0;
      setCenter(C);
      P2M(C,Rmax);
      M2M(C,Rmax);
    }
#if Cartesian
    for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {
      for( int i=1; i<MTERM; ++i ) C->M[i] /= C->M[0];
    }
#endif
    setRcrit(cells);
    stopTimer("Upward pass  ",printNow);
  }

  void downwardPass(Cells &cells) const {
    for( C_iter C=cells.begin()+1; C!=cells.end(); ++C ) {
      L2L(C);
      L2P(C);
    }
  }

  void printTreeData(Cells &cells) {
    setRootCell(cells);
    std::cout << "------------------------" << std::endl;
    std::cout << "Root center   : " << ROOT->X              << std::endl;
    std::cout << "Root radius   : " << R0                   << std::endl;
    std::cout << "Bodies        : " << ROOT->NDLEAF         << std::endl;
    std::cout << "Cells         : " << cells.size()         << std::endl;
    std::cout << "Tree depth    : " << MAXLEVEL             << std::endl;
    std::cout << "Total charge  : " << std::abs(ROOT->M[0]) << std::endl;
    std::cout << "P2P calls     : " << NP2P                 << std::endl;
    std::cout << "M2P calls     : " << NM2P                 << std::endl;
    std::cout << "M2L calls     : " << NM2L                 << std::endl;
    std::cout << "------------------------" << std::endl;
  }

};

#endif
