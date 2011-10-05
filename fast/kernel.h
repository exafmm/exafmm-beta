#ifndef kernel_h
#define kernel_h
#include "sort.h"

class Kernel : public Sort {
private:
  real DMAX;

protected:
  vect   X0;
  real   R0;
  B_iter B0, BN;
  C_iter C0, CN;

private:
  real getBmax(vect const&X, C_iter C) {
    real rad = C->R;
    real dx = rad+std::abs(X[0]-C->X[0]);
    real dy = rad+std::abs(X[1]-C->X[1]);
    real dz = rad+std::abs(X[2]-C->X[2]);
    return std::sqrt( dx*dx + dy*dy + dz*dz );
  }

  inline void set_dPhi(Lset &D, vect const&dX, real const&m) const {
    real R2 = norm(dX);
    real invR2 = 1.0 / (R2 + EPS2);
    real invR  = m * std::sqrt(invR2);
    real invR3 =     invR2 * invR;
    real invR5 = 3 * invR2 * invR3;
    real invR7 = 5 * invR2 * invR5;
    D[ 0] = invR;
    real t = -invR3;
    D[ 1] = t * dX[0];
    D[ 2] = t * dX[1];
    D[ 3] = t * dX[2];
    t     = invR5 * dX[0];
    D[ 4] = t * dX[0] - invR3;
    D[ 5] = t * dX[1];
    D[ 6] = t * dX[2];
    t     = invR5 * dX[1];
    D[ 7] = t * dX[1] - invR3;
    D[ 8] = t * dX[2];
    t     = invR5 * dX[2];
    D[ 9] = t * dX[2] - invR3;
    real s = 3 * invR5;
    t     = -invR7 * dX[0] * dX[0];
    D[10] = (s     + t) * dX[0];
    D[11] = (invR5 + t) * dX[1];
    D[12] = (invR5 + t) * dX[2];
    t     = -invR7 * dX[1] * dX[1];
    D[13] = (invR5 + t) * dX[0];
    D[16] = (s     + t) * dX[1];
    D[17] = (invR5 + t) * dX[2];
    t     = -invR7 * dX[2] * dX[2];
    D[15] = (invR5 + t) * dX[0];
    D[18] = (invR5 + t) * dX[1];
    D[19] = (s     + t) * dX[2];
    D[14] = -invR7 * dX[0] * dX[1] * dX[2];
  }

  inline void add_C_C2C(Lset &L, Lset const&D, Mset const&m) const {
    L += D;
    L[0] += D[4] *m[1] + D[5] *m[2] + D[6] *m[3] + D[7] *m[4] + D[8] *m[5] + D[9] *m[6];
    L[1] += D[10]*m[1] + D[11]*m[2] + D[12]*m[3] + D[13]*m[4] + D[14]*m[5] + D[15]*m[6];
    L[2] += D[11]*m[1] + D[13]*m[2] + D[14]*m[3] + D[16]*m[4] + D[17]*m[5] + D[18]*m[6];
    L[3] += D[12]*m[1] + D[14]*m[2] + D[15]*m[3] + D[17]*m[4] + D[18]*m[5] + D[19]*m[6];
  }

  inline void flip_sign_odd (Lset &a) const {
    for( int i=1; i!=4; ++i ) a[i] = -a[i];
    for( int i=10; i!=20; ++i ) a[i] = -a[i];
  }

public:
  Kernel() : X0(0), R0(0) {}
  ~Kernel() {}

  void setCenter(C_iter C) {
    DMAX = 0;
    real m = 0;
    vect X = 0;
    for( B_iter B=C->LEAF; B!=C->LEAF+C->NCLEAF; ++B ) {
      m += B->SRC[0];
      X += B->X * B->SRC[0];
    }
    for( C_iter c=C0+C->CHILD; c!=C0+C->CHILD+C->NCHILD; ++c ) {
      m += c->M[0];
      X += c->X * c->M[0];
    }
    X /= m;
    C->R = getBmax(X,C);
    C->X = X;
  }

  void P2P(C_iter Ci, C_iter Cj, bool mutual=true) const {
    for( B_iter BI=Ci->LEAF; BI!=Ci->LEAF+Ci->NDLEAF; ++BI ) {
      real P0 = 0;
      vect F0 = 0;
      for( B_iter BJ=Cj->LEAF; BJ!=Cj->LEAF+Cj->NDLEAF; ++BJ ) {
        vect dR = BI->X - BJ->X;
        real D1 = norm(dR) + EPS2;
        real D0 = BI->SRC[0] * BJ->SRC[0];
        real XX = 1.0/D1;
        D0 *= std::sqrt(XX);
        D1  = XX * D0;
        dR *= D1;
        P0 -= D0;
        F0 -= dR;
        BJ->TRG[0] -= D0 * mutual;
        BJ->TRG[1] += dR[0] * mutual;
        BJ->TRG[2] += dR[1] * mutual;
        BJ->TRG[3] += dR[2] * mutual;
      }
      BI->TRG[0] += P0;
      BI->TRG[1] += F0[0];
      BI->TRG[2] += F0[1];
      BI->TRG[3] += F0[2];
    }
  }

  void P2P(C_iter C) const {
    unsigned NJ = C->NDLEAF;
    for( B_iter BI=C->LEAF; BI!=C->LEAF+C->NDLEAF; ++BI, --NJ ) {
      real P0 = 0;
      vect F0 = 0;
      for( B_iter BJ=BI+1; BJ!=BI+NJ; ++BJ ) {
        vect dR = BI->X - BJ->X;
        real D1 = norm(dR) + EPS2;
        real D0 = BI->SRC[0] * BJ->SRC[0];
        real XX = 1.0/D1;
        D0 *= std::sqrt(XX);
        D1  = XX * D0;
        dR *= D1;
        P0 -= D0;
        F0 -= dR;
        BJ->TRG[0] -= D0;
        BJ->TRG[1] += dR[0];
        BJ->TRG[2] += dR[1];
        BJ->TRG[3] += dR[2];
      }
      BI->TRG[0] += P0;
      BI->TRG[1] += F0[0];
      BI->TRG[2] += F0[1];
      BI->TRG[3] += F0[2];
    }
  }

  void P2M(C_iter C) {
    for( B_iter B=C->LEAF; B!=C->LEAF+C->NCLEAF; ++B ) {
      vect dX = B->X - C->X;
      real R = std::sqrt(norm(dX));
      if( R > DMAX ) DMAX = R;
      real tmp = B->SRC[0] * dX[0];
      C->M[0] += B->SRC[0];
      C->M[1] += dX[0] * tmp;
      C->M[2] += dX[1] * tmp;
      C->M[3] += dX[2] * tmp;
      tmp = B->SRC[0] * dX[1];
      C->M[4] += dX[1] * tmp;
      C->M[5] += dX[2] * tmp;
      C->M[6] += B->SRC[0] * dX[2] * dX[2];
    }
    C->RCRIT = std::min(C->R,DMAX);
  }

  void M2M(C_iter C) {
    for( C_iter c=C0+C->CHILD; c!=C0+C->CHILD+C->NCHILD; ++c ) {
      vect dX = c->X - C->X;
      real R = std::sqrt(norm(dX)) + c->RCRIT;
      if( R > DMAX ) DMAX = R;
      for( int i=0; i!=6; ++i ) C->M[i] += c->M[i];
      real tmp = c->M[0] * dX[0];
      C->M[1] += dX[0] * tmp;
      C->M[2] += dX[1] * tmp;
      C->M[3] += dX[2] * tmp;
      tmp = c->M[0] * dX[1];
      C->M[4] += dX[1] * tmp;
      C->M[5] += dX[2] * tmp;
      C->M[6] += c->M[0] * dX[2] * dX[2];
    }
    C->RCRIT = std::min(C->R,DMAX);
  }

  void M2L(C_iter Ci, C_iter Cj, bool mutual=true) const {
    vect dX = Ci->X - Cj->X;
    Lset D;
    set_dPhi(D,dX,Ci->M[0]*Cj->M[0]);
    add_C_C2C(Ci->L,D,Cj->M);
    if( mutual ) {
      flip_sign_odd(D);
      add_C_C2C(Cj->L,D,Ci->M);
    }
  }

  void L2L(C_iter C) const {
    vect dX = C->X - (C0+C->PARENT)->X;
    Lset L = (C0+C->PARENT)->L;
    Lset o;
    o[0] = L[1] *dX[0] + L[2] *dX[1] + L[3] *dX[2];
    o[1] = L[4] *dX[0] + L[5] *dX[1] + L[6] *dX[2];
    o[2] = L[5] *dX[0] + L[7] *dX[1] + L[8] *dX[2];
    o[3] = L[6] *dX[0] + L[8] *dX[1] + L[9] *dX[2];
    o[4] = L[10]*dX[0] + L[11]*dX[1] + L[12]*dX[2];
    o[5] = L[11]*dX[0] + L[13]*dX[1] + L[14]*dX[2];
    o[6] = L[12]*dX[0] + L[14]*dX[1] + L[15]*dX[2];
    o[7] = L[13]*dX[0] + L[16]*dX[1] + L[17]*dX[2];
    o[8] = L[14]*dX[0] + L[17]*dX[1] + L[18]*dX[2];
    o[9] = L[15]*dX[0] + L[18]*dX[1] + L[19]*dX[2];
    for( int i=0; i<10; i++ ) L[i] += o[i];
    o[0] = (o[1]*dX[0] + o[2]*dX[1] + o[3]*dX[2]) / 2;
    o[1] = (o[4]*dX[0] + o[5]*dX[1] + o[6]*dX[2]) / 2;
    o[2] = (o[5]*dX[0] + o[7]*dX[1] + o[8]*dX[2]) / 2;
    o[3] = (o[6]*dX[0] + o[8]*dX[1] + o[9]*dX[2]) / 2;
    L[0]+=o[0] + (dX[0]*o[1]+dX[1]*o[2]+dX[2]*o[3]) / 3;
    L[1]+=o[1];  L[2]+=o[2];  L[3]+=o[3];
    L += C->L / C->M[0];
    C->L = L;
  }

  void L2P(C_iter C) const {
    for( B_iter B=C->LEAF; B!=C->LEAF+C->NCLEAF; ++B ) {
      Lset o;
      vect dX = B->X - C->X;
      B->TRG /= B->SRC[0];
      B->TRG[0] -= C->L[0];
      B->TRG[1] += C->L[1];
      B->TRG[2] += C->L[2];
      B->TRG[3] += C->L[3];
      o[0] = C->L[1] *dX[0] + C->L[2] *dX[1] + C->L[3] *dX[2];
      o[1] = C->L[4] *dX[0] + C->L[5] *dX[1] + C->L[6] *dX[2];
      o[2] = C->L[5] *dX[0] + C->L[7] *dX[1] + C->L[8] *dX[2];
      o[3] = C->L[6] *dX[0] + C->L[8] *dX[1] + C->L[9] *dX[2];
      o[4] = C->L[10]*dX[0] + C->L[11]*dX[1] + C->L[12]*dX[2];
      o[5] = C->L[11]*dX[0] + C->L[13]*dX[1] + C->L[14]*dX[2];
      o[6] = C->L[12]*dX[0] + C->L[14]*dX[1] + C->L[15]*dX[2];
      o[7] = C->L[13]*dX[0] + C->L[16]*dX[1] + C->L[17]*dX[2];
      o[8] = C->L[14]*dX[0] + C->L[17]*dX[1] + C->L[18]*dX[2];
      o[9] = C->L[15]*dX[0] + C->L[18]*dX[1] + C->L[19]*dX[2];
      B->TRG[0] -= o[0];
      B->TRG[1] += o[1];
      B->TRG[2] += o[2];
      B->TRG[3] += o[3];
      o[0] = (o[1]*dX[0] + o[3]*dX[2] + o[2]*dX[1]) / 2;
      o[1] = (o[4]*dX[0] + o[6]*dX[2] + o[5]*dX[1]) / 2;
      o[2] = (o[5]*dX[0] + o[8]*dX[2] + o[7]*dX[1]) / 2;
      o[3] = (o[6]*dX[0] + o[9]*dX[2] + o[8]*dX[1]) / 2;
      B->TRG[0] -= o[0] + (dX[0]*o[1]+dX[1]*o[2]+dX[2]*o[3]) / 3;
      B->TRG[1] += o[1];
      B->TRG[2] += o[2];
      B->TRG[3] += o[3];
    }
  }
};

#endif
