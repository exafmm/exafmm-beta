#ifndef tree_h
#define tree_h
#include <evaluator.h>

struct Dot {
  int I;
  vect X;
  Dot *NEXT;
};

struct Node : public Dot {
  int LEVEL;
  int ICHILD;
  int NLEAF;
  Node *CHILD[8];
  Node *DOTS;

  void mark_as_box   (int const&i)       { ICHILD |= (1<<i); }
  bool marked_as_box (int const&i) const { return ICHILD & (1<<i); }
  bool marked_as_dot (int const&i) const { return !marked_as_box(i); }

  int octant(const Node*D) const {
    int oct = 0;
    for( int d=0; d!=3; ++d ) {
      oct += (D->X[d] > X[d]) << d;
    }
    return oct;
  }

  bool is_twig() const {
    return ICHILD == 0;
  }

  Node *init() {
    ICHILD   = 0;
    NLEAF    = 0;
    DOTS     = 0;
    for( int i=0; i!=8; ++i ) CHILD[i] = NULL;
    return this;
  }
};

class TreeBuilder {
public:
  int  LEVEL;
  int  NLEAF;
  int  NCELL;
  real RAD;
  Node *START;
  Node *N0, *NN;
  Dot  *D0, *DN;
  B_iter B0, BN;
  Cell *C0, *CN;
  vect XAVE, XMIN, XMAX;

private:
  inline real root_radius(const vect& x) const {
    real R,D=zero;
    for(int d=0; d!=3; ++d) {
      R=std::max(std::abs(XMAX[d]-x[d]),std::abs(XMIN[d]-x[d]));
      if(R>D) D=R;
    }
    R=pow(2.0,int(1.0+log(D)/M_LN2));
    return R;
  }

  void set_domain(Bodies &bodies) {
    D0     = new Dot [bodies.size()];
    Dot*Di = D0;
    XAVE = zero;
    XMAX = XMIN = bodies.begin()->X;
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      Di->I = B-bodies.begin();
      Di->X = B->X;
      Di->X.min_max(XMIN,XMAX);
      XAVE += Di->X;
      Di++;
    }
    DN    = Di;
    XAVE /= real(DN-D0);
  }

  inline Node* make_subbox(const Node*B, int i) {
    ++NCELL;
    assert(NCELL < NLEAF);
    Node* N = (NN++)->init();
    N->LEVEL    = B->LEVEL+1;
    real r = RAD / (1 << N->LEVEL);
    N->X = B->X;
    if(i&1) N->X[0] += r;  else  N->X[0] -= r;
    if(i&2) N->X[1] += r;  else  N->X[1] -= r;
    if(i&4) N->X[2] += r;  else  N->X[2] -= r;
    return N;
  }

  void split_box(Node* N) {
    int NUM[8];
    int b,ne;
    Node *sub=0;
    Node *Di,*Dn;
    do {
      for(b=0; b!=8; ++b) NUM[b] = 0;
      for(Di=N->DOTS; Di; Di=Dn) {
        Dn = static_cast<Node*>(Di->NEXT);
        b = N->octant(Di);
        Di->NEXT = N->CHILD[b];
        N->CHILD[b] = Di;
        NUM[b]++;
      }
      N->DOTS = 0;
      for(ne=b=0; b!=8; ++b) if(NUM[b]) {
        ne++;
        sub = make_subbox(N,b);
        sub->DOTS = N->CHILD[b];
        sub->NLEAF = NUM[b];
        N->CHILD[b] = sub;
        N->mark_as_box(b);
      }
      N = sub;
    } while(ne==1);
  }

  void adddot_N(Node *const&Di) {
    for(Node*N=N0;;) {
      if(N->is_twig()) {
        Di->NEXT = N->DOTS;
        N->DOTS = Di;
        N->NLEAF++;
        if(N->NLEAF>NCRIT) split_box(N);
        if(LEVEL < N->LEVEL) LEVEL = N->LEVEL;
        return;
      } else {
        int b = N->octant(Di);
        N->NLEAF++;
        if( N->CHILD[b] == NULL ) {
          N->mark_as_box(b);
          Node *sub = make_subbox(N,b);
          N->CHILD[b] = sub;
          N = sub;
        } else
          N = N->CHILD[b];
      }
    }
  }

  void link_cells_N(const Node* N, Cell *C) {
    C->LEVEL  = N->LEVEL;
    C->X      = N->X;
    C->NDLEAF = N->NLEAF;
    C->LEAF   = BN;
    if(N->is_twig()) {
      C->FCCELL = C0;
      C->NCCELL = 0;
      C->NCLEAF = N->NLEAF;
      for(Node *Di=N->DOTS; Di; Di=static_cast<Node*>(Di->NEXT)) {
        BN->IBODY = Di->I;
        BN->X = Di->X;
        BN++;
      }
    } else {
      C->NCLEAF = 0;
      int nsub=0;
      for( int i=0; i!=8; ++i ) if(N->CHILD[i]) {
        if(N->marked_as_box(i)) ++nsub;
        else {
          BN->IBODY = N->CHILD[i]->I;
          BN->X = N->CHILD[i]->X;
          BN++;
          C->NCLEAF++;
        }
      }
      if(nsub) {
        Cell *Ci = CN;
        C->FCCELL = Ci;
        C->NCCELL = nsub;
        CN += nsub;
        for( int i=0; i!=8; ++i ) {
          if(N->CHILD[i] && N->marked_as_box(i)) {
            Ci->PARENT = C;
            link_cells_N(N->CHILD[i], Ci++);
          }
        }
      } else {
        C->FCCELL = C0;
        C->NCCELL = 0;
      }
    }
  }

public:
  TreeBuilder(Bodies& bodies) : RAD(0), N0(0), C0(NULL) {
    set_domain(bodies);
    vect X0(zero);
    for(int d=0; d!=3; ++d) X0[d]=int(XAVE[d]+0.5);
    NCELL     = 1;
    NLEAF     = bodies.size();
    START     = new Node [NLEAF];
    NN        = START;
    N0        = (NN++)->init();
    RAD       = root_radius(X0);
    N0->LEVEL = 0;
    N0->X     = X0;
  }
  ~TreeBuilder() {
    delete[] START;
    delete[] D0;
  }

  void build() {
    double tic,toc;
    tic = get_time();
    LEVEL = 0;
    for(Dot *Di=D0; Di!=DN; ++Di) {
      adddot_N(static_cast<Node*>(Di));
    }
    LEVEL++;
    toc = get_time();
    std::cout << "build  : " << toc-tic << std::endl;
  }

  void link(Cell *c0, Bodies &leafs) {
    double tic,toc;
    tic = get_time();
    C0 = c0;
    B0 = leafs.begin();
    CN = C0+1;
    BN = B0;
    link_cells_N(N0,C0);
    toc = get_time();
    std::cout << "link   : " << toc-tic << std::endl;
  }
};

#endif // tree_h
