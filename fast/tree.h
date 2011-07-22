#ifndef tree_h
#define tree_h
#include <evaluator.h>

class TreeBuilder : public Evaluator {
private:
  Nodes nodes;
  Cells cells;
  B_iter B0, BN;
  Leaf *L0, *LN;
  Cell *CN;
  vect XAVE, XMIN, XMAX;

private:
  void init(Node &N) {
    N.ICHILD = 0;
    N.NLEAF = 0;
    N.LEAF = NULL;
    for( int b=0; b!=8; ++b ) N.CHILD[b] = -1;
  }

  inline real root_radius(const vect& x) const {
    real R,D=zero;
    for(int d=0; d!=3; ++d) {
      R=std::max(std::abs(XMAX[d]-x[d]),std::abs(XMIN[d]-x[d]));
      if(R>D) D=R;
    }
    R=pow(2.0,int(1.0+log(D)/M_LN2));
    return R;
  }

  int getOctant(const Leaf *L, N_iter N) const {
    int octant = 0;
    for( int d=0; d!=3; ++d ) {
      octant += (L->X[d] > N->X[d]) << d;
    }
    return octant;
  }

  void set_domain(Bodies &bodies) {
    L0 = new Leaf [bodies.size()];
    Leaf *Li = L0;
    XAVE = zero;
    XMAX = XMIN = bodies.begin()->X;
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      Li->I = B-bodies.begin();
      Li->X = B->X;
      Li->X.min_max(XMIN,XMAX);
      XAVE += Li->X;
      Li++;
    }
    LN = Li;
    XAVE /= real(LN-L0);
  }

  inline void addChild(int octant, N_iter node) {
    assert(nodes.size() < NLEAF);
    Node child;
    init(child);
    child.LEVEL = node->LEVEL+1;
    real r = RAD / (1 << child.LEVEL);
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

  void nodes2cells(int i, Cell *C) {
    C->R      = RAD / ( 1 << nodes[i].LEVEL );
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
        Cell *Ci = CN;
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
  TreeBuilder(Bodies &bodies) : Evaluator(bodies) {}
  ~TreeBuilder() {
    delete[] L0;
  }

  void build() {
    set_domain(BODIES);
    vect X0(zero);
    for(int d=0; d!=3; ++d) X0[d]=int(XAVE[d]+0.5);
    NCELL = 1;
    NLEAF = BODIES.size();
    RAD   = root_radius(X0);
    nodes.reserve(NLEAF);
    Node node;
    init(node);
    node.LEVEL = 0;
    node.X     = X0;
    nodes.push_back(node);
    LEVEL = 0;
    for( Leaf *Li=L0; Li!=LN; ++Li ) {
      int i = 0;
      N_iter N = nodes.begin()+i;
      while( N->ICHILD != 0 ) {
        int octant = getOctant(Li,N);
        N->NLEAF++;
        if( N->CHILD[octant] == -1 ) {
          addChild(octant,N);
        }
        i = N->CHILD[octant];
        N = nodes.begin()+i;
      }
      Li->NEXT = N->LEAF;
      N->LEAF = Li;
      N->NLEAF++;
      if( N->NLEAF > NCRIT ) splitNode(N);
      if( LEVEL < N->LEVEL ) LEVEL = N->LEVEL;
    }
    LEVEL++;
  }

  void link() {
    LEAFS.resize(NLEAF);
    C0 = new Cell [NCELL];
    B0 = LEAFS.begin();
    CN = C0+1;
    BN = B0;
    nodes2cells(0,C0);
  }
};

#endif
