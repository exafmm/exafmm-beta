#ifndef topdown_h
#define topdown_h
#include "evaluator.h"

class TreeBuilder : public Evaluator {
private:
  Nodes    nodes;
  Leafs    leafs;
  B_iter   BN;
  C_iter   CN;
  unsigned NLEAF;
  unsigned NCELL;
  int      MAXLEVEL;

private:
  void init(Node &node) {
    node.NOCHILD = true;
    node.NLEAF = 0;
    node.LEAF = NULL;
    for( int b=0; b!=8; ++b ) node.CHILD[b] = -1;
  }

  inline long long getIndex(vect X, int level) {
    float d = 2 * R0 / (1 << level);
    long long ix = int((X[0] + R0 - X0[0]) / d);
    long long iy = int((X[1] + R0 - X0[1]) / d);
    long long iz = int((X[2] + R0 - X0[2]) / d);
    long long id = ((1 << 3 * level) - 1) / 7;
    assert( level < 22 );
    for( int l=0; l<level; l++ ) {
      id += (ix & 1) << (3 * l);
      id += (iy & 1) << (3 * l + 1);
      id += (iz & 1) << (3 * l + 2);
      ix >>= 1;
      iy >>= 1;
      iz >>= 1;
    }
    return id;
  }

  inline void addChild(int octant, int n) {
    N_iter N = nodes.begin()+n;
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

  void splitNode(int n) {
    if( nodes[n].NLEAF > NCRIT ) {
      int c = 0;
      Leaf *Ln;
      for( Leaf *L=nodes[n].LEAF; L; L=Ln ) {
	N_iter N = nodes.begin()+n;
        Ln = L->NEXT;
        int octant = (L->X[0] > N->X[0]) + ((L->X[1] > N->X[1]) << 1) + ((L->X[2] > N->X[2]) << 2);
        if( N->CHILD[octant] == -1 ) {
          addChild(octant,n);
        }
        c = nodes[n].CHILD[octant];
        Node *child = &nodes[c];
        L->NEXT = child->LEAF;
        child->LEAF = &*L;
        child->NLEAF++;
      }
      n = c;
    }
  }

  void nodes2cells(int i, C_iter C) {
    C->R      = R0 / (1 << nodes[i].LEVEL);
    C->X      = nodes[i].X;
    C->NDLEAF = nodes[i].NLEAF;
    C->LEAF   = BN;
    C->ICELL  = getIndex(C->X,nodes[i].LEVEL);
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
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      B->ICELL = buffer[B->IBODY].ICELL;
      B->X     = buffer[B->IBODY].X;
      B->SRC   = buffer[B->IBODY].SRC;
      B->TRG   = buffer[B->IBODY].TRG;
    }
  }

protected:
  void setLeafs(Bodies &bodies) {
    startTimer("Set leafs");
    NLEAF = bodies.size();
    leafs.reserve(NLEAF);
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      Leaf leaf;
      leaf.I = B-bodies.begin();
      leaf.X = B->X;
      leafs.push_back(leaf);
    }
    stopTimer("Set leafs",printNow);
  }

  void growTree() {
    startTimer("Grow tree");
    NCELL = 1;
    nodes.reserve(NLEAF);
    Node node;
    init(node);
    node.LEVEL = 0;
    node.X     = X0;
    nodes.push_back(node);
    MAXLEVEL = 0;
    for( L_iter L=leafs.begin(); L!=leafs.end(); ++L ) {
      int n = 0;
      while( !nodes[n].NOCHILD ) {
	N_iter N = nodes.begin() + n;
        int octant = (L->X[0] > N->X[0]) + ((L->X[1] > N->X[1]) << 1) + ((L->X[2] > N->X[2]) << 2);
        N->NLEAF++;
        if( nodes[n].CHILD[octant] == -1 ) addChild(octant,n);
        n = nodes[n].CHILD[octant];
      }
      L->NEXT = nodes[n].LEAF;
      nodes[n].LEAF = &*L;
      nodes[n].NLEAF++;
      if( nodes[n].NLEAF > NCRIT ) splitNode(n);
      if( MAXLEVEL < nodes[n].LEVEL ) MAXLEVEL = nodes[n].LEVEL;
    }
    MAXLEVEL++;
    stopTimer("Grow tree",printNow);
  }

  void linkTree(Bodies &bodies, Cells &cells) {
    startTimer("Link tree");
    cells.clear();
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

  void printTreeData(Cells &cells) {
    setRootCell(cells);
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "Root center          : " << X0                   << std::endl;
    std::cout << "Root radius          : " << R0                   << std::endl;
    std::cout << "Bodies               : " << Ci0->NDLEAF          << std::endl;
    std::cout << "Cells                : " << cells.size()         << std::endl;
    std::cout << "Tree depth           : " << MAXLEVEL             << std::endl;
    std::cout << "Total charge         : " << std::abs(Ci0->M[0])  << std::endl;
    std::cout << "P2P calls            : " << NP2P                 << std::endl;
    std::cout << "M2P calls            : " << NM2P                 << std::endl;
    std::cout << "M2L calls            : " << NM2L                 << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;
  }
};

#endif
