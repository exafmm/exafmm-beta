#ifndef tree_h
#define tree_h
#include <evaluator.h>

class TreeBuilder {
private:
  struct Dot {
    int I;
    vect X;
    Dot *NEXT;
  };

  struct Node {
    int  LEVEL;
    int  ICHILD;
    int  NLEAF;
    int  CHILD[8];
    vect X;
    Dot  *DOTS;
  };

  Node *N0;
  int   NN;
  Dot  *D0, *DN;
  B_iter B0, BN;
  Cell *C0, *CN;
  vect XAVE, XMIN, XMAX;

public:
  int  LEVEL;
  int  NLEAF;
  int  NCELL;
  real RAD;

private:
  void init(int i) {
    N0[i].ICHILD = 0;
    N0[i].NLEAF = 0;
    N0[i].DOTS = NULL;
    for( int b=0; b!=8; ++b ) N0[i].CHILD[b] = -1;
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

  int getOctant(const Dot *D, const int i) const {
    int octant = 0;
    for( int d=0; d!=3; ++d ) {
      octant += (D->X[d] > N0[i].X[d]) << d;
    }
    return octant;
  }

  void set_domain(Bodies &bodies) {
    D0 = new Dot [bodies.size()];
    Dot *Di = D0;
    XAVE = zero;
    XMAX = XMIN = bodies.begin()->X;
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      Di->I = B-bodies.begin();
      Di->X = B->X;
      Di->X.min_max(XMIN,XMAX);
      XAVE += Di->X;
      Di++;
    }
    DN = Di;
    XAVE /= real(DN-D0);
  }

  inline int addChild(int octant, int i) {
    ++NCELL;
    assert(NCELL < NLEAF);
    int c = NN;
    init(NN++);
    N0[c].LEVEL = N0[i].LEVEL+1;
    real r = RAD / (1 << N0[c].LEVEL);
    N0[c].X = N0[i].X;
    for( int d=0; d!=3; ++d ) {                                 // Loop over dimensions
      N0[c].X[d] += r * (((octant & 1 << d) >> d) * 2 - 1);        //  Calculate new center position
    }                                                           // End loop over dimensions
    return c;
  }

  void splitNode(int &i) {
    while( N0[i].NLEAF > NCRIT ) {
      int c;
      Dot *Dn;
      for( Dot *Di=N0[i].DOTS; Di; Di=Dn ) {
        Dn = Di->NEXT;
        int octant = getOctant(Di,i);
        c = N0[i].CHILD[octant];
        if( c == -1 ) c = addChild(octant,i);
        Di->NEXT = N0[c].DOTS;
        N0[c].DOTS = Di;
        N0[c].NLEAF++;
        N0[i].CHILD[octant] = c;
        N0[i].ICHILD |= 1 << octant;
      }
      i = c;
    }
  }

  void link_cells_N(const int i, Cell *C) {
    C->R      = RAD / ( 1 << N0[i].LEVEL );
    C->X      = N0[i].X;
    C->NDLEAF = N0[i].NLEAF;
    C->LEAF   = BN;
    if( N0[i].ICHILD == 0 ) {
      C->CHILD = 0;
      C->NCHILD = 0;
      C->NCLEAF = N0[i].NLEAF;
      for(Dot *Di=N0[i].DOTS; Di; Di=Di->NEXT) {
        BN->IBODY = Di->I;
        BN->X = Di->X;
        BN++;
      }
    } else {
      C->NCLEAF = 0;
      int nsub=0;
      for( int d=0; d!=8; ++d ) if(N0[i].CHILD[d] != -1) {
        ++nsub;
      }
      if(nsub) {
        Cell *Ci = CN;
        C->CHILD = Ci - C0;
        C->NCHILD = nsub;
        CN += nsub;
        for( int d=0; d!=8; ++d ) {
          if(N0[i].CHILD[d] != -1) {
            Ci->PARENT = C - C0;
            link_cells_N(N0[i].CHILD[d], Ci++);
          }
        }
      } else {
        C->CHILD = 0;
        C->NCHILD = 0;
      }
    }
  }

public:
  TreeBuilder() : N0(0), C0(NULL) {}
  ~TreeBuilder() {
    delete[] N0;
    delete[] D0;
  }

  void build(Bodies& bodies) {
    set_domain(bodies);
    vect X0(zero);
    for(int d=0; d!=3; ++d) X0[d]=int(XAVE[d]+0.5);
    NCELL     = 1;
    NLEAF     = bodies.size();
    N0        = new Node [NLEAF];
    NN        = 0;
    init(NN++);
    RAD       = root_radius(X0);
    N0->LEVEL = 0;
    N0->X     = X0;
    LEVEL = 0;
    for(Dot *Di=D0; Di!=DN; ++Di) {
      int i = 0;
      while( N0[i].ICHILD != 0 ) {
        int octant = getOctant(Di,i);
        N0[i].NLEAF++;
        if( N0[i].CHILD[octant] == -1 ) {
          N0[i].ICHILD |= 1 << octant;
          N0[i].CHILD[octant] = addChild(octant,i);
        }
        i = N0[i].CHILD[octant];
      }
      Di->NEXT = N0[i].DOTS;
      N0[i].DOTS = Di;
      N0[i].NLEAF++;
      if( N0[i].NLEAF > NCRIT ) splitNode(i);
      if( LEVEL < N0[i].LEVEL ) LEVEL = N0[i].LEVEL;
    }
    LEVEL++;
  }

  void link(Cell *c0, Bodies &leafs) {
    C0 = c0;
    B0 = leafs.begin();
    CN = C0+1;
    BN = B0;
    link_cells_N(0,C0);
  }
};

#endif
