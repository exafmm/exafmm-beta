#ifndef kernel_h
#define kernel_h
#include <body.h>

class Kernel {
protected:
  unsigned NLEAFS;
  unsigned NCELLS;
  real RAD;

public:
  Leaf     *L0;
  Cell     *C0;
  bodies   &BODIES;

private:
  inline real getBmax(vect const&X, Cell *C) {
    real rad = RAD / ( 1 << C->LEVEL );
    real dx = rad+std::abs(X[0]-C->X[0]);
    real dy = rad+std::abs(X[1]-C->X[1]);
    real dz = rad+std::abs(X[2]-C->X[2]);
    return std::sqrt( dx*dx + dy*dy + dz*dz );
  }

  inline real getCenter(Cell *C) {
    real m(zero);
    vect X(zero);
    for( Cell *c=C->FCCELL; c!=C->FCCELL+C->NCCELL; ++c ) {
      m += c->M[0];
      X += c->X * c->M[0];
    }
    for( Leaf *l=C->FCLEAF; l!=C->FCLEAF+C->NCLEAF; ++l ) {
      m += l->Q;
      X += l->X * l->Q;
    }
    X /= m;
    real bmax = getBmax(X,C);
    C->X = X;
    return bmax;
  }

  inline void set_D(real&X, real D[5]) const {
    D[0] *= std::sqrt(X);
    D[1]  =     X * D[0];
    D[2]  = 3 * X * D[1];
    D[3]  = 5 * X * D[2];
    D[4]  = 7 * X * D[3];
  }

  inline void set_dPhi(Lset &C, vect const&v, real const D[4]) const {
    C[ 0] = D[0];
    real t =-D[1];
    C[ 1] = v[0]*t;
    C[ 2] = v[1]*t;
    C[ 3] = v[2]*t;
    t     = D[2]*v[0];
    C[ 4] = t*v[0] - D[1];
    C[ 5] = t*v[1];
    C[ 6] = t*v[2];
    t     = D[2]*v[1];
    C[ 7] = t*v[1] - D[1];
    C[ 8] = t*v[2];
    t     = D[2]*v[2];
    C[ 9] = t*v[2] - D[1];
    real D2_3 = 3*D[2];
    t     = D[3]*v[0]*v[0];
    C[10] = (D2_3-t)*v[0];
    C[11] = (D[2]-t)*v[1];
    C[12] = (D[2]-t)*v[2];
    t     =  D[3]*v[1]*v[1];
    C[13] = (D[2]-t)*v[0];
    C[16] = (D2_3-t)*v[1];
    C[17] = (D[2]-t)*v[2];
    t     =  D[3]*v[2]*v[2];
    C[15] = (D[2]-t)*v[0];
    C[18] = (D[2]-t)*v[1];
    C[19] = (D2_3-t)*v[2];
    C[14] =-D[3]*v[0]*v[1]*v[2];
  }

  inline void add_C_C2C(Lset &cd, Lset const&f, Mset const&m) const {
    cd += f;
    cd[0] += f[4]*m[1] + f[9] *m[6] + f[7] *m[4] +
        2*( f[5] *m[2] + f[6] *m[3] + f[8] *m[5] );
    cd[1]+= f[10]*m[1] + f[15]*m[6] + f[13]*m[4] +
        2*( f[11]*m[2] + f[12]*m[3] + f[14]*m[5] );
    cd[2]+= f[11]*m[1] + f[18]*m[6] + f[16]*m[4] +
        2*( f[13]*m[2] + f[14]*m[3] + f[17]*m[5] );
    cd[3]+= f[12]*m[1] + f[19]*m[6] + f[17]*m[4] +
        2*( f[14]*m[2] + f[15]*m[3] + f[18]*m[5] );
  }

  inline void flip_sign_odd (Lset &a) const {
    for( int i=1; i!=4; ++i ) a[i] = -a[i];
    for( int i=10; i!=20; ++i ) a[i] = -a[i];
  }

public:
  Kernel(bodies &b, real rad, unsigned nleafs, unsigned ncells) :
    NLEAFS     (nleafs),
    NCELLS     (ncells),
    RAD        (rad),
    BODIES     ( b )
  {
    L0 = new Leaf [NLEAFS];
    C0 = new Cell [NCELLS];
  }
  ~Kernel() {
    delete[] L0;
    delete[] C0;
  }

  void P2M(Cell *C, real &dmax, real &bmax) {
    bmax = getCenter(C);
    for( Leaf *l=C->FCLEAF; l!=C->FCLEAF+C->NCLEAF; ++l ) {
      vect dX = l->X - C->X;
      if(norm(dX)>dmax) dmax=norm(dX);
      real tmp = l->Q * dX[0];
      C->M[0] += l->Q;
      C->M[1] += dX[0] * tmp;
      C->M[2] += dX[1] * tmp;
      C->M[3] += dX[2] * tmp;
      tmp = l->Q * dX[1];
      C->M[4] += dX[1] * tmp;
      C->M[5] += dX[2] * tmp;
      C->M[6] += l->Q * dX[2] * dX[2];
    }
    if(C->NCLEAF != 0) dmax = std::sqrt(dmax);
    C->RCRIT = std::min(dmax,bmax);
  }

  void M2M(Cell *C, real &dmax, real &bmax) {
    for( Cell *c=C->FCCELL; c!=C->FCCELL+C->NCCELL; ++c ) {
      vect dX = c->X - C->X;
      real Xq = norm(dX);
      real x  = dmax - c->RCRIT;
      if(zero>x || Xq>x*x)
        dmax = std::sqrt(Xq) + c->RCRIT;
      for( int i=1; i!=6; ++i ) C->M[i] += c->M[i];
      real tmp = c->M[0] * dX[0];
      C->M[0] += c->M[0];
      C->M[1] += dX[0] * tmp;
      C->M[2] += dX[1] * tmp;
      C->M[3] += dX[2] * tmp;
      tmp = c->M[0] * dX[1];
      C->M[4] += dX[1] * tmp;
      C->M[5] += dX[2] * tmp;
      C->M[6] += c->M[0] * dX[2] * dX[2];
    }
    C->RCRIT = std::min(dmax,bmax);
  }

  void P2P(Cell *Ci, Cell *Cj, bool mutual=true) const {
    for( Leaf *Li=Ci->FCLEAF; Li!=Ci->FCLEAF+Ci->NDLEAF; ++Li ) {
      real P0(zero);
      vect F0(zero);
      for( Leaf *Lj=Cj->FCLEAF; Lj!=Cj->FCLEAF+Cj->NDLEAF; ++Lj ) {
        vect dR = Li->X - Lj->X;
        real D1 = norm(dR) + EQ;
        real D0 = Li->Q * Lj->Q;
        real XX = 1.0/D1;
        D0 *= std::sqrt(XX);
        D1  = XX * D0;
        dR *= D1;
        P0 -= D0;
        F0 -= dR;
        Lj->TRG[0] -= D0 * mutual;
        Lj->TRG[1] += dR[0] * mutual;
        Lj->TRG[2] += dR[1] * mutual;
        Lj->TRG[3] += dR[2] * mutual;
      }
      Li->TRG[0] += P0;
      Li->TRG[1] += F0[0];
      Li->TRG[2] += F0[1];
      Li->TRG[3] += F0[2];
    }
  }

  void P2P(Cell *C) const {
    unsigned NJ = C->NDLEAF;
    for( Leaf *Li=C->FCLEAF; Li!=C->FCLEAF+C->NDLEAF; ++Li, --NJ ) {
      real P0(zero);
      vect F0(zero);
      for(Leaf *Lj=Li+1; Lj!=Li+NJ; ++Lj) {
        vect dR = Li->X - Lj->X;
        real D1 = norm(dR) + EQ;
        real D0 = Li->Q * Lj->Q;
        real XX = 1.0/D1;
        D0 *= std::sqrt(XX);
        D1  = XX * D0;
        dR *= D1;
        P0 -= D0;
        F0 -= dR;
        Lj->TRG[0] -= D0;
        Lj->TRG[1] += dR[0];
        Lj->TRG[2] += dR[1];
        Lj->TRG[3] += dR[2];
      }
      Li->TRG[0] += P0;
      Li->TRG[1] += F0[0];
      Li->TRG[2] += F0[1];
      Li->TRG[3] += F0[2];
    }
  }

  void M2L(Cell *Ci, Cell *Cj, bool mutual=true) const {
    real D[5];
    vect dX = Ci->X - Cj->X;
    real R2 = norm(dX);
    real XX = 1.0 / (R2 + EQ);
    D[0] = Ci->M[0] * Cj->M[0];
    set_D(XX,D);
    Lset F;
    set_dPhi(F,dX,D);
    add_C_C2C(Ci->L,F,Cj->M);
    if( mutual ) {
      flip_sign_odd(F);
      add_C_C2C(Cj->L,F,Ci->M);
    }
  }

  void L2L(Cell *C) const {
    vect dX = C->X - C->PARENT->X;
    Lset L = C->PARENT->L;
    Lset o;
    o[0] = L[1]*dX[0] + L[3]*dX[2] + L[2]*dX[1];
    o[1] = L[4]*dX[0] + L[6]*dX[2] + L[5]*dX[1];
    o[2] = L[5]*dX[0] + L[8]*dX[2] + L[7]*dX[1];
    o[3] = L[6]*dX[0] + L[9]*dX[2] + L[8]*dX[1];
    o[4] = L[10]*dX[0] + L[12]*dX[2] + L[11]*dX[1];
    o[5] = L[11]*dX[0] + L[14]*dX[2] + L[13]*dX[1];
    o[6] = L[12]*dX[0] + L[15]*dX[2] + L[14]*dX[1];
    o[7] = L[13]*dX[0] + L[17]*dX[2] + L[16]*dX[1];
    o[8] = L[14]*dX[0] + L[18]*dX[2] + L[17]*dX[1];
    o[9] = L[15]*dX[0] + L[19]*dX[2] + L[18]*dX[1];
    for( int i=0; i<10; i++ ) L[i] += o[i];
    o[0] = (o[1]*dX[0] + o[3]*dX[2] + o[2]*dX[1]) / 2;
    o[1] = (o[4]*dX[0] + o[6]*dX[2] + o[5]*dX[1]) / 2;
    o[2] = (o[5]*dX[0] + o[8]*dX[2] + o[7]*dX[1]) / 2;
    o[3] = (o[6]*dX[0] + o[9]*dX[2] + o[8]*dX[1]) / 2;
    L[0]+=o[0] + (dX[0]*o[1]+dX[1]*o[2]+dX[2]*o[3]) / 3;
    L[1]+=o[1];  L[2]+=o[2];  L[3]+=o[3];
    L += C->L / C->M[0];
    C->L = L;
  }

  void L2P(Cell *C) const {
    for( Leaf *l=C->FCLEAF; l!=C->FCLEAF+C->NCLEAF; ++l ) {
      Lset o;
      vect dX = l->X - C->X;
      l->TRG /= l->Q;
      l->TRG[0] -= C->L[0];
      l->TRG[1] += C->L[1];
      l->TRG[2] += C->L[2];
      l->TRG[3] += C->L[3];
      o[0] = C->L[1] *dX[0] + C->L[3] *dX[2] + C->L[2] *dX[1];
      o[1] = C->L[4] *dX[0] + C->L[6] *dX[2] + C->L[5] *dX[1];
      o[2] = C->L[5] *dX[0] + C->L[8] *dX[2] + C->L[7] *dX[1];
      o[3] = C->L[6] *dX[0] + C->L[9] *dX[2] + C->L[8] *dX[1];
      o[4] = C->L[10]*dX[0] + C->L[12]*dX[2] + C->L[11]*dX[1];
      o[5] = C->L[11]*dX[0] + C->L[14]*dX[2] + C->L[13]*dX[1];
      o[6] = C->L[12]*dX[0] + C->L[15]*dX[2] + C->L[14]*dX[1];
      o[7] = C->L[13]*dX[0] + C->L[17]*dX[2] + C->L[16]*dX[1];
      o[8] = C->L[14]*dX[0] + C->L[18]*dX[2] + C->L[17]*dX[1];
      o[9] = C->L[15]*dX[0] + C->L[19]*dX[2] + C->L[18]*dX[1];
      l->TRG[0] -= o[0];
      l->TRG[1] += o[1];
      l->TRG[2] += o[2];
      l->TRG[3] += o[3];
      o[0] = (o[1]*dX[0] + o[3]*dX[2] + o[2]*dX[1]) / 2;
      o[1] = (o[4]*dX[0] + o[6]*dX[2] + o[5]*dX[1]) / 2;
      o[2] = (o[5]*dX[0] + o[8]*dX[2] + o[7]*dX[1]) / 2;
      o[3] = (o[6]*dX[0] + o[9]*dX[2] + o[8]*dX[1]) / 2;
      l->TRG[0] -= o[0] + (dX[0]*o[1]+dX[1]*o[2]+dX[2]*o[3]) / 3;
      l->TRG[1] += o[1];
      l->TRG[2] += o[2];
      l->TRG[3] += o[3];
    }
  }
};

#endif
