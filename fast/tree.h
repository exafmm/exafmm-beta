#ifndef tree_h
#define tree_h
#include <evaluator.h>

class TreeConstructor : public Evaluator {
private:
  Nodes nodes;
  Cells cells;
  Leafs leafs;
  Bodies BODIES;

private:
  void init(Node &N) {
    N.ICHILD = 0;
    N.NLEAF = 0;
    N.LEAF = NULL;
    for( int b=0; b!=8; ++b ) N.CHILD[b] = -1;
  }

  int getOctant(const Leaf *L, N_iter N) const {
    int octant = 0;
    for( int d=0; d!=3; ++d ) {
      octant += (L->X[d] > N->X[d]) << d;
    }
    return octant;
  }

  inline void addChild(int octant, N_iter node) {
    assert(nodes.size() < NLEAF);
    Node child;
    init(child);
    child.LEVEL = node->LEVEL+1;
    real r = R0 / (1 << child.LEVEL);
    child.X = node->X;
    for( int d=0; d!=3; ++d ) {                                 // Loop over dimensions
      child.X[d] += r * (((octant & 1 << d) >> d) * 2 - 1);     //  Calculate new center position
    }                                                           // End loop over dimensions
    node->ICHILD |= 1 << octant;
    node->CHILD[octant] = nodes.size();
    nodes.push_back(child);
    NCELL++;
  }

  void splitNode(N_iter node) {
    while( node->NLEAF > NCRIT ) {
      int c;
      Leaf *Ln;
      for( Leaf *Li=node->LEAF; Li; Li=Ln ) {
        Ln = Li->NEXT;
        int octant = getOctant(Li,node);
        if( !(node->ICHILD & (1 << octant)) ) {
          addChild(octant,node);
        }
        c = node->CHILD[octant];
        Node *child = &nodes[c];
        Li->NEXT = child->LEAF;
        child->LEAF = Li;
        child->NLEAF++;
      }
      node = nodes.begin()+c;
    }
  }

  void nodes2cells(int i, C_iter C) {
    C->R      = R0 / ( 1 << nodes[i].LEVEL );
    C->X      = nodes[i].X;
    C->NDLEAF = nodes[i].NLEAF;
    C->LEAF   = BN;
    if( nodes[i].ICHILD == 0 ) {
      C->CHILD = 0;
      C->NCHILD = 0;
      C->NCLEAF = nodes[i].NLEAF;
      for( Leaf *Li=nodes[i].LEAF; Li; Li=Li->NEXT ) {
        BN->IBODY = Li->I;
        BN->X = Li->X;
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
      if(nsub) {
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
      } else {
        C->CHILD = 0;
        C->NCHILD = 0;
      }
    }
  }

public:
  TreeConstructor() {}
  ~TreeConstructor() {}

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
    R0 += 1e-5;                                                 // Add some leeway to root radius
  }

  void bodies2leafs(Bodies &bodies) {
    for( B_iter B=BODIES.begin(); B!=BODIES.end(); ++B ) {      // Loop over bodies
      B->SRC[0] = bodies[B->IBODY].SRC[0];
      B->TRG = 0;
    }
  }

  void leafs2bodies(Bodies &bodies) {
    for( B_iter B=BODIES.begin(); B!=BODIES.end(); ++B ) {      // Loop over bodies
      bodies[B->IBODY].TRG[0] = B->TRG[0];
      bodies[B->IBODY].TRG[1] = B->TRG[1];
      bodies[B->IBODY].TRG[2] = B->TRG[2];
      bodies[B->IBODY].TRG[3] = B->TRG[3];
    }
  }

  void build() {
    NCELL = 1;
    nodes.reserve(NLEAF);
    Node node;
    init(node);
    node.LEVEL = 0;
    node.X     = X0;
    nodes.push_back(node);
    LEVEL = 0;
    for( L_iter Li=leafs.begin(); Li!=leafs.end(); ++Li ) {
      int i = 0;
      N_iter N = nodes.begin()+i;
      while( N->ICHILD != 0 ) {
        int octant = getOctant(&*Li,N);
        N->NLEAF++;
        if( N->CHILD[octant] == -1 ) {
          addChild(octant,N);
        }
        i = N->CHILD[octant];
        N = nodes.begin()+i;
      }
      Li->NEXT = N->LEAF;
      N->LEAF = &*Li;
      N->NLEAF++;
      if( N->NLEAF > NCRIT ) splitNode(N);
      if( LEVEL < N->LEVEL ) LEVEL = N->LEVEL;
    }
    LEVEL++;
  }

  void link() {
    BODIES.resize(NLEAF);
    cells.resize(NCELL);
    C0 = cells.begin();
    B0 = BODIES.begin();
    CN = C0+1;
    BN = B0;
    nodes2cells(0,C0);
  }

  void exact(Bodies &bodies) {
    bodies2leafs(bodies);
#if IeqJ
    P2P(C0);
#else
    P2P(C0,C0,false);
#endif
    for( B_iter B=BODIES.begin(); B!=BODIES.end(); ++B ) {      // Loop over bodies
      B->TRG /= B->SRC[0];
    }
    leafs2bodies(bodies);
  }

  void approximate(Bodies &bodies) {
    double tic,toc;
    tic = get_time();
    bodies2leafs(bodies);
    upward();
    toc = get_time();
    std::cout << "upward : " << toc-tic << std::endl;
    tic = get_time();
#if IeqJ
    traverse();
#else
    traverse(false);
#endif
    toc = get_time();
    std::cout << "intrct : " << toc-tic << std::endl;
    tic = get_time();
    for( C_iter C=C0+C0->CHILD; C!=C0+C0->CHILD+C0->NCHILD; ++C ) {
      downward(C);
    }
    toc = get_time();
    std::cout << "downwd : " << toc-tic << std::endl;
#ifndef MANY
    write();
#endif
    leafs2bodies(bodies);
  }
};

#endif
