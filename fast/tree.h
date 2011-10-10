#ifndef tree_h
#define tree_h
#include "evaluator.h"

class TopDown : public Evaluator {
private:
  Nodes nodes;
  Leafs leafs;
  B_iter BN;

private:
  void init(Node &node) {
    node.ICHILD = 0;
    node.NLEAF = 0;
    node.LEAF = NULL;
    for( int b=0; b!=8; ++b ) node.CHILD[b] = -1;
  }

  inline void addChild(int octant, N_iter N) {
    assert(nodes.size() < NLEAF);
    Node child;
    init(child);
    child.LEVEL = N->LEVEL+1;
    real r = R0 / (1 << child.LEVEL);
    child.X = N->X;
    for( int d=0; d!=3; ++d ) {                                 // Loop over dimensions
      child.X[d] += r * (((octant & 1 << d) >> d) * 2 - 1);     //  Calculate new center position
    }                                                           // End loop over dimensions
    N->ICHILD |= 1 << octant;
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
        if( !(N->ICHILD & (1 << octant)) ) {
          addChild(octant,N);
        }
        c = N->CHILD[octant];
        Node *child = &nodes[c];
        L->NEXT = child->LEAF;
        child->LEAF = L;
        child->NLEAF++;
      }
      N = nodes.begin()+c;
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
      for( Leaf *L=nodes[i].LEAF; L; L=L->NEXT ) {
        BN->IBODY = L->I;
        BN->X = L->X;
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

protected:
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

  void build() {
    startTimer("Grow tree    ");
    NCELL = 1;
    nodes.reserve(NLEAF);
    Node node;
    init(node);
    node.LEVEL = 0;
    node.X     = X0;
    nodes.push_back(node);
    LEVEL = 0;
    for( L_iter L=leafs.begin(); L!=leafs.end(); ++L ) {
      int i = 0;
      N_iter N = nodes.begin()+i;
      while( N->ICHILD != 0 ) {
        int octant = (L->X[0] > N->X[0]) + ((L->X[1] > N->X[1]) << 1) + ((L->X[2] > N->X[2]) << 2);
        N->NLEAF++;
        if( N->CHILD[octant] == -1 ) {
          addChild(octant,N);
        }
        i = N->CHILD[octant];
        N = nodes.begin()+i;
      }
      L->NEXT = N->LEAF;
      N->LEAF = &*L;
      N->NLEAF++;
      if( N->NLEAF > NCRIT ) splitNode(N);
      if( LEVEL < N->LEVEL ) LEVEL = N->LEVEL;
    }
    LEVEL++;
    stopTimer("Grow tree    ",printNow);
  }

  void link() {
    startTimer("Link tree    ");
    BODIES.resize(NLEAF);
    cells.resize(NCELL);
    C0 = cells.begin();
    BN = BODIES.begin();
    CN = C0+1;
    nodes2cells(0,C0);
    CN = C0;
    nodes.clear();
    leafs.clear();
    stopTimer("Link tree    ",printNow);
  }

  void upward() {
    for( C_iter C=C0; C!=C0+NCELL; ++C ) {
      C->M = 0;
      C->L = 0;
    }
    for( C_iter C=C0+NCELL-1; C!=C0-1; --C ) {
      setCenter(C);
      P2M(C);
      M2M(C);
    }
#if CART
#elif SPHE
#else
    for( C_iter C=C0; C!=C0+NCELL; ++C ) {
      C->M[1] *= 0.5 / C->M[0];
      C->M[2] *= 1.0 / C->M[0];
      C->M[3] *= 1.0 / C->M[0];
      C->M[4] *= 0.5 / C->M[0];
      C->M[5] *= 1.0 / C->M[0];
      C->M[6] *= 0.5 / C->M[0];
    }
#endif
    set_rcrit();
  }

public:
  TopDown() {}
  ~TopDown() {}
};

class BottomUp : public TopDown {
private:
  int getMaxLevel(Bodies &bodies) {
    const long N = bodies.size();
    int level;
    level = N >= NCRIT ? 1 + int(log(N / NCRIT)/M_LN2/3) : 0;
    return level;
  }

  inline void getIndex() {
    float d = 2 * R0 / (1 << LEVEL);
    for( B_iter B=BODIES.begin(); B!=BODIES.end(); ++B ) {
      int ix = (B->X[0] + R0 - X0[0]) / d;
      int iy = (B->X[1] + R0 - X0[1]) / d;
      int iz = (B->X[2] + R0 - X0[2]) / d;
      int id = 0;
      for( int l=0; l!=LEVEL; ++l ) {
        id += ix % 2 << (3 * l);
        id += iy % 2 << (3 * l + 1);
        id += iz % 2 << (3 * l + 2);
        ix >>= 1;
        iy >>= 1;
        iz >>= 1;
      }
      B->ICELL = id;
    }
  }

  void bodies2twigs() {
    int I = -1;
    C_iter C;
    cells.reserve(1 << (3 * LEVEL));
    float d = 2 * R0 / (1 << LEVEL);
    for( B_iter B=BODIES.begin(); B!=BODIES.end(); ++B ) {
      int IC = B->ICELL;
      int ix = (B->X[0] + R0 - X0[0]) / d;
      int iy = (B->X[1] + R0 - X0[1]) / d;
      int iz = (B->X[2] + R0 - X0[2]) / d;
      if( IC != I ) {
        Cell cell;
        cell.NCHILD = 0;
        cell.NCLEAF = 0;
        cell.NDLEAF = 0;
        cell.CHILD  = 0;
        cell.LEAF   = B;
        cell.X[0]   = d * (ix + .5) + X0[0] - R0;
        cell.X[1]   = d * (iy + .5) + X0[1] - R0;
        cell.X[2]   = d * (iz + .5) + X0[2] - R0;
        cell.R      = d * .5;
        cells.push_back(cell);
        C = cells.end()-1;
        I = IC;
      }
      C->NCLEAF++;
      C->NDLEAF++;
    }
  }

  void twigs2cells() {
    int begin = 0, end = cells.size();
    float d = 2 * R0 / (1 << LEVEL);
    for( int l=0; l!=LEVEL; ++l ) {
      int div = (8 << (3 * l));
      int I = -1;
      int p = end - 1;
      d *= 2;
      for( int c=begin; c!=end; ++c ) {
        B_iter B = cells[c].LEAF;
        int IC = B->ICELL / div;
        int ix = (B->X[0] + R0 - X0[0]) / d;
        int iy = (B->X[1] + R0 - X0[1]) / d;
        int iz = (B->X[2] + R0 - X0[2]) / d;
        if( IC != I ) {
          Cell cell;
          cell.NCHILD = 0;
          cell.NCLEAF = 0;
          cell.NDLEAF = 0;
          cell.CHILD  = c;
          cell.LEAF   = cells[c].LEAF;
          cell.X[0]   = d * (ix + .5) + X0[0] - R0;
          cell.X[1]   = d * (iy + .5) + X0[1] - R0;
          cell.X[2]   = d * (iz + .5) + X0[2] - R0;
          cell.R      = d * .5;
          cells.push_back(cell);
          p++;
          I = IC;
        }
        cells[p].NCHILD++;
        cells[p].NDLEAF += cells[c].NDLEAF;
        cells[c].PARENT = p;
      }
      begin = end;
      end = cells.size();
    }
  }

protected:
  void setDomain(Bodies &bodies) {
    BODIES = bodies;
    LEVEL = getMaxLevel(bodies);
    vect xmin, xmax;
    X0 = 0;
    xmax = xmin = bodies.begin()->X;
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
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

  void build() {
    startTimer("Grow tree    ");
    getIndex();
    Bodies buffer = BODIES;
    sortBodies(BODIES,buffer);
    bodies2twigs();
    stopTimer("Grow tree    ",printNow);
  }

  void link() {
    startTimer("Link tree    ");
    twigs2cells();
    NCELL = cells.size();
    C0 = cells.begin();
    CN = cells.end()-1;
    stopTimer("Link tree    ",printNow);
  }

  void upward() {
    for( C_iter C=C0; C!=C0+NCELL; ++C ) {
      C->M = 0;
      C->L = 0;
    }
    for( C_iter C=C0; C!=C0+NCELL; ++C ) {
      setCenter(C);
      P2M(C);
      M2M(C);
    }
#if CART
#elif SPHE
#else
    for( C_iter C=C0; C!=C0+NCELL; ++C ) {
      C->M[1] *= 0.5 / C->M[0];
      C->M[2] *= 1.0 / C->M[0];
      C->M[3] *= 1.0 / C->M[0];
      C->M[4] *= 0.5 / C->M[0];
      C->M[5] *= 1.0 / C->M[0];
      C->M[6] *= 0.5 / C->M[0];
    }
#endif
    set_rcrit();
  }

public:
  BottomUp() {}
  ~BottomUp() {}
};

class TreeConstructor : public BottomUp {
private:
  void bodies2leafs(Bodies &bodies) {
    for( B_iter B=BODIES.begin(); B!=BODIES.end(); ++B ) {      // Loop over bodies
      B->SRC = bodies[B->IBODY].SRC;
      B->TRG = 0;
    }
  }

  void leafs2bodies(Bodies &bodies) {
    for( B_iter B=BODIES.begin(); B!=BODIES.end(); ++B ) {      // Loop over bodies
      bodies[B->IBODY].TRG = B->TRG;
    }
  }

public:
  void topdown(Bodies &bodies) {
    TopDown::setDomain(bodies);
    TopDown::build();
    TopDown::link();
    startTimer("Upward       ");
    bodies2leafs(bodies);
    TopDown::upward();
    stopTimer("Upward       ",printNow);
  }

  void bottomup(Bodies &bodies) {
    BottomUp::setDomain(bodies);
    BottomUp::build();
    BottomUp::link();
    startTimer("Upward       ");
    bodies2leafs(bodies);
    BottomUp::upward();
    stopTimer("Upward       ",printNow);
  }

  void exact(Bodies &bodies) {
    bodies2leafs(bodies);
#if IEQJ
    P2P(CN);
#else
    P2P(CN,CN,false);
#endif
    for( B_iter B=BODIES.begin(); B!=BODIES.end(); ++B ) {      // Loop over bodies
      B->TRG /= B->SRC[0];
    }
    leafs2bodies(bodies);
  }

  void approximate(Bodies &bodies) {
    startTimer("Traverse     ");
#if IEQJ
    traverse();
#else
    traverse(false);
#endif
    stopTimer("Traverse     ",printNow);
    startTimer("Downward     ");
    for( C_iter C=C0+CN->CHILD; C!=C0+CN->CHILD+CN->NCHILD; ++C ) {
      downward(C);
    }
#ifndef MANY
    write();
#endif
    leafs2bodies(bodies);
    cells.clear();
    stopTimer("Downward     ",printNow);
  }
};

#endif
