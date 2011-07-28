#ifndef kernel_h
#define kernel_h
#include <types.h>

class Kernel {
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

  void P2M(C_iter C) {
    for( B_iter B=C->LEAF; B!=C->LEAF+C->NCLEAF; ++B ) {
      vect dist = C->X - B->X;
      real R = std::sqrt(norm(dist));
      if( R > DMAX ) DMAX = R;
      C->M[0] += B->SRC[0];
      C->M[1] += B->SRC[0] * dist[0];
      C->M[2] += B->SRC[0] * dist[1];
      C->M[3] += B->SRC[0] * dist[2];
      C->M[4] += B->SRC[0] * dist[0] * dist[0] / 2;
      C->M[5] += B->SRC[0] * dist[1] * dist[1] / 2;
      C->M[6] += B->SRC[0] * dist[2] * dist[2] / 2;
      C->M[7] += B->SRC[0] * dist[0] * dist[1];
      C->M[8] += B->SRC[0] * dist[1] * dist[2];
      C->M[9] += B->SRC[0] * dist[2] * dist[0];
    }
    C->RCRIT = std::min(C->R,DMAX);
  }

  void M2M(C_iter CI) {
    for( C_iter CJ=C0+CI->CHILD; CJ!=C0+CI->CHILD+CI->NCHILD; ++CJ ) {
      vect dist = CI->X - CJ->X;
      real R = std::sqrt(norm(dist)) + CJ->RCRIT;
      if( R > DMAX ) DMAX = R;
      CI->M[0] += CJ->M[0];
      CI->M[1] += CJ->M[1] + dist[0] * CJ->M[0];
      CI->M[2] += CJ->M[2] + dist[1] * CJ->M[0];
      CI->M[3] += CJ->M[3] + dist[2] * CJ->M[0];
      CI->M[4] += CJ->M[4] + dist[0] * CJ->M[1] + dist[0] * dist[0]  * CJ->M[0] / 2;
      CI->M[5] += CJ->M[5] + dist[1] * CJ->M[2] + dist[1] * dist[1]  * CJ->M[0] / 2;
      CI->M[6] += CJ->M[6] + dist[2] * CJ->M[3] + dist[2] * dist[2]  * CJ->M[0] / 2;
      CI->M[7] += CJ->M[7] + dist[0] * CJ->M[2] + dist[1] * CJ->M[1] + dist[0] * dist[1] * CJ->M[0];
      CI->M[8] += CJ->M[8] + dist[1] * CJ->M[3] + dist[2] * CJ->M[2] + dist[1] * dist[2] * CJ->M[0];
      CI->M[9] += CJ->M[9] + dist[2] * CJ->M[1] + dist[0] * CJ->M[3] + dist[2] * dist[0] * CJ->M[0];
    }
    CI->RCRIT = std::min(CI->R,DMAX);
  }

  void P2P(C_iter CI, C_iter CJ, bool mutual=true) const {
    for( B_iter BI=CI->LEAF; BI!=CI->LEAF+CI->NDLEAF; ++BI ) {
      real P0 = 0;
      vect F0 = 0;
      for( B_iter BJ=CJ->LEAF; BJ!=CJ->LEAF+CJ->NDLEAF; ++BJ ) {
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

  void M2L(C_iter CI, C_iter CJ, bool mutual=true) const {
    vect dist = CI->X - CJ->X;
    real invR = 1 / std::sqrt(norm(dist));
    real invR3 = invR * invR * invR;
    real invR5 = invR3 * invR * invR;
    CI->L[0] += CJ->M[0] * invR;
    CI->L[0] += CJ->M[1] * (-dist[0] * invR3);
    CI->L[0] += CJ->M[2] * (-dist[1] * invR3);
    CI->L[0] += CJ->M[3] * (-dist[2] * invR3);
    CI->L[0] += CJ->M[4] * (3 * dist[0] * dist[0] * invR5 - invR3);
    CI->L[0] += CJ->M[5] * (3 * dist[1] * dist[1] * invR5 - invR3);
    CI->L[0] += CJ->M[6] * (3 * dist[2] * dist[2] * invR5 - invR3);
    CI->L[0] += CJ->M[7] * (3 * dist[0] * dist[1] * invR5);
    CI->L[0] += CJ->M[8] * (3 * dist[1] * dist[2] * invR5);
    CI->L[0] += CJ->M[9] * (3 * dist[2] * dist[0] * invR5);
    CI->L[1] += CJ->M[0] * (-dist[0] * invR3);
    CI->L[1] += CJ->M[1] * (3 * dist[0] * dist[0] * invR5 - invR3);
    CI->L[1] += CJ->M[2] * (3 * dist[0] * dist[1] * invR5);
    CI->L[1] += CJ->M[3] * (3 * dist[0] * dist[2] * invR5);
    CI->L[2] += CJ->M[0] * (-dist[1] * invR3);
    CI->L[2] += CJ->M[1] * (3 * dist[1] * dist[0] * invR5);
    CI->L[2] += CJ->M[2] * (3 * dist[1] * dist[1] * invR5 - invR3);
    CI->L[2] += CJ->M[3] * (3 * dist[1] * dist[2] * invR5);
    CI->L[3] += CJ->M[0] * (-dist[2] * invR3);
    CI->L[3] += CJ->M[1] * (3 * dist[2] * dist[0] * invR5);
    CI->L[3] += CJ->M[2] * (3 * dist[2] * dist[1] * invR5);
    CI->L[3] += CJ->M[3] * (3 * dist[2] * dist[2] * invR5 - invR3);
    CI->L[4] += CJ->M[0] * (3 * dist[0] * dist[0] * invR5 - invR3) / 2;
    CI->L[5] += CJ->M[0] * (3 * dist[1] * dist[1] * invR5 - invR3) / 2;
    CI->L[6] += CJ->M[0] * (3 * dist[2] * dist[2] * invR5 - invR3) / 2;
    CI->L[7] += CJ->M[0] * (3 * dist[0] * dist[1] * invR5);
    CI->L[8] += CJ->M[0] * (3 * dist[1] * dist[2] * invR5);
    CI->L[9] += CJ->M[0] * (3 * dist[2] * dist[0] * invR5);
  }

  void L2L(C_iter CI) const {
    C_iter CJ = C0 + CI->PARENT;
    vect dist = CI->X - CJ->X;
    for( int i=0; i<10; ++i )
      CI->L[i] += CJ->L[i];
    CI->L[0] += CJ->L[1] * dist[0];
    CI->L[0] += CJ->L[2] * dist[1];
    CI->L[0] += CJ->L[3] * dist[2];
    CI->L[0] += CJ->L[4] * dist[0] * dist[0] / 2;
    CI->L[0] += CJ->L[5] * dist[1] * dist[1] / 2;
    CI->L[0] += CJ->L[6] * dist[2] * dist[2] / 2;
    CI->L[0] += CJ->L[7] * dist[0] * dist[1];
    CI->L[0] += CJ->L[8] * dist[1] * dist[2];
    CI->L[0] += CJ->L[9] * dist[2] * dist[0];
    CI->L[1] += CJ->L[4] * dist[0];
    CI->L[1] += CJ->L[7] * dist[1];
    CI->L[1] += CJ->L[9] * dist[2];
    CI->L[2] += CJ->L[7] * dist[0];
    CI->L[2] += CJ->L[5] * dist[1];
    CI->L[2] += CJ->L[8] * dist[2];
    CI->L[3] += CJ->L[9] * dist[0];
    CI->L[3] += CJ->L[8] * dist[1];
    CI->L[3] += CJ->L[6] * dist[2];
  }

  void L2P(C_iter CI) const {
    for( B_iter B=CI->LEAF; B!=CI->LEAF+CI->NCLEAF; ++B ) {
      vect dist = B->X - CI->X;
      B->TRG /= B->SRC[0];
      B->TRG[0] -= CI->L[0];
      B->TRG[0] -= CI->L[1] * dist[0];
      B->TRG[0] -= CI->L[2] * dist[1];
      B->TRG[0] -= CI->L[3] * dist[2];
      B->TRG[0] -= CI->L[4] * dist[0] * dist[0] / 2;
      B->TRG[0] -= CI->L[5] * dist[1] * dist[1] / 2;
      B->TRG[0] -= CI->L[6] * dist[2] * dist[2] / 2;
      B->TRG[0] -= CI->L[7] * dist[0] * dist[1];
      B->TRG[0] -= CI->L[8] * dist[1] * dist[2];
      B->TRG[0] -= CI->L[9] * dist[2] * dist[0];
      B->TRG[1] += CI->L[1];
      B->TRG[1] += CI->L[4] * dist[0];
      B->TRG[1] += CI->L[7] * dist[1];
      B->TRG[1] += CI->L[9] * dist[2];
      B->TRG[2] += CI->L[2];
      B->TRG[2] += CI->L[7] * dist[0];
      B->TRG[2] += CI->L[5] * dist[1];
      B->TRG[2] += CI->L[8] * dist[2];
      B->TRG[3] += CI->L[3];
      B->TRG[3] += CI->L[9] * dist[0];
      B->TRG[3] += CI->L[8] * dist[1];
      B->TRG[3] += CI->L[6] * dist[2];
    }
  }
};

#endif
