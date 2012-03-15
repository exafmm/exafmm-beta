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

template<Equation equation>
class TopDown : public Evaluator<equation> {
private:
  Nodes    nodes;
  Leafs    leafs;
  B_iter   BN;
  C_iter   CN;
  unsigned NLEAF;
  unsigned NCELL;

public:
  using Kernel<equation>::printNow;                             //!< Switch to print timings
  using Kernel<equation>::startTimer;                           //!< Start timer for given event
  using Kernel<equation>::stopTimer;                            //!< Stop timer for given event
  using Kernel<equation>::X0;                                   //!< Center of root cell
  using Kernel<equation>::R0;                                   //!< Radius of root cell
  using Kernel<equation>::Ci0;                                  //!< icells.begin()
  using Evaluator<equation>::MAXLEVEL;                          //!< Max depth of tree

private:
  inline void addChild(int octant, N_iter N) {
    Node child;
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
        if( N->CHILD[octant] == 0 ) {
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
    C->RCRIT  = C->R / THETA;
    C->X      = nodes[i].X;
    C->NDLEAF = nodes[i].NLEAF;
    C->LEAF   = BN;
    float diameter = 2 * C->R;
    int ix = int((nodes[i].X[0] + R0 - X0[0]) / diameter);
    int iy = int((nodes[i].X[1] + R0 - X0[1]) / diameter);
    int iz = int((nodes[i].X[2] + R0 - X0[2]) / diameter);
    C->ICELL  = getMorton(ix,iy,iz,nodes[i].LEVEL);
    if( nodes[i].NOCHILD ) {
      C->NCLEAF = nodes[i].NLEAF;
      for( Leaf *L=nodes[i].LEAF; L; L=L->NEXT ) {
        BN->IBODY = L->I;
        BN++;
      }
    } else {
      int nsub=0;
      for( int octant=0; octant!=8; ++octant ) {
        if( nodes[i].CHILD[octant] != 0 ) {
          ++nsub;
        }
      }
      C_iter Ci = CN;
      C->CHILD = Ci - Ci0;
      C->NCHILD = nsub;
      CN += nsub;
      for( int octant=0; octant!=8; ++octant ) {
        if( nodes[i].CHILD[octant] != 0 ) {
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
  inline int getMorton(int ix, int iy, int iz, const int &level) const {
    int id = 0;
    for( int l=0; l!=level; ++l ) {
      id += ix % 2 << (3 * l);
      id += iy % 2 << (3 * l + 1);
      id += iz % 2 << (3 * l + 2);
      ix >>= 1;
      iy >>= 1;
      iz >>= 1;
    }
    return id;
  }

  void setDomain(Bodies &bodies) {
    startTimer("Set domain");
    vect xmin, xmax;
    NLEAF = bodies.size();
    leafs.reserve(NLEAF);
    X0 = 0;
    xmax = xmin = bodies.begin()->X;
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      Leaf leaf;
      leaf.I = B-bodies.begin();
      leaf.X = B->X;
      leafs.push_back(leaf);
      for( int d=0; d!=3; ++d ) {
        if     (B->X[d] < xmin[d]) xmin[d] = B->X[d];
        else if(B->X[d] > xmax[d]) xmax[d] = B->X[d];
      }
      X0 += B->X;
    }
    X0 /= bodies.size();
    for( int d=0; d!=3; ++d ) {
      X0[d] = int(X0[d]+.5);
      R0 = std::max(xmax[d] - X0[d], R0);
      R0 = std::max(X0[d] - xmin[d], R0);
    }
    R0 *= 1.000001;
    if( IMAGES != 0 ) {
      X0 = 0;
      R0 = M_PI;
    }
    stopTimer("Set domain",printNow);
  }

  void buildTree() {
    startTimer("Grow tree");
    NCELL = 1;
    nodes.reserve(NLEAF);
    Node node;
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
        if( N->CHILD[octant] == 0 ) addChild(octant,N);
        N = nodes.begin()+N->CHILD[octant];
      }
      L->NEXT = N->LEAF;
      N->LEAF = &*L;
      N->NLEAF++;
      if( N->NLEAF > NCRIT ) splitNode(N);
      if( MAXLEVEL < N->LEVEL ) MAXLEVEL = N->LEVEL;
    }
    MAXLEVEL++;
    stopTimer("Grow tree",printNow);
  }

  void linkTree(Bodies &bodies, Cells &cells) {
    startTimer("Link tree");
    cells.resize(NCELL);
    Ci0 = cells.begin();
    BN = bodies.begin();
    CN = Ci0 + 1;
    nodes2cells(0,Ci0);
    nodes.clear();
    leafs.clear();
    permuteBodies(bodies);
    stopTimer("Link tree",printNow);
  }

public:
//! Constructor
  TopDown() : nodes(), leafs(), BN(), CN(), NLEAF(0), NCELL(0) {}

};

#endif
