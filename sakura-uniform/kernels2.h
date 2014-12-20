#include <cmath>
#include "types2.h"
#include "core.h"

#ifdef SAKURA
#define DIM 3
#define NN 26
#define FN 37
#define CN 152
int executeStencilNN(int *nn_codes, int *node_code, int *stencil, int periodic, int lv);
int executeStencilFN(int *fn_codes, int *node_code, int *fn_stencil, int *cn_stencil, int periodic, int lv);
#endif

class Kernel {
public:
  int maxLevel;
  int numBodies;
  int numCells;
  int numLeafs;
  int numNeighbors;
  int numImages;

  real_t X0[3];
  real_t R0;
  real_t (*Ibodies)[4];
  real_t (*Ibodies2)[4];
  real_t (*Jbodies)[4];
  real_t (*Multipole)[MTERM];
  real_t (*Local)[LTERM];
  int (*Leafs)[2];
  int *Index;
  int *Index2;

private:
  inline void getIndex(int *ix, int index) const {
    for_3 ix[d] = 0;
    int d = 0, level = 0;
    while (index != 0) {
      ix[d] += (index & 1) * (1 << level);
      index >>= 1;
      d = (d+1) % 3;
      if (d == 0) level++;
    }
  }

  void getCenter(real_t *dX, int index, int level) const {
    real_t R = R0 / (1 << level);
    int ix[3] = {0, 0, 0};
    getIndex(ix, index);
    for_3 dX[d] = X0[d] - R0 + (2 * ix[d] + 1) * R;
  }

  void P2PSum(int ibegin, int iend, int jbegin, int jend, real_t *Xperiodic) const {
    for (int i=ibegin; i<iend; i++) {
      real_t Po = 0, Fx = 0, Fy = 0, Fz = 0;
      for (int j=jbegin; j<jend; j++) {
	real_t dX[3];
	for_3 dX[d] = Jbodies[i][d] - Jbodies[j][d] - Xperiodic[d];
	real_t R2 = dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2];
	real_t invR2 = R2 == 0 ? 0 : 1.0 / R2;
	real_t invR = Jbodies[j][3] * sqrt(invR2);
	real_t invR3 = invR2 * invR;
	Po += invR;
	Fx += dX[0] * invR3;
	Fy += dX[1] * invR3;
	Fz += dX[2] * invR3;
      }
      Ibodies[i][0] += Po;
      Ibodies[i][1] -= Fx;
      Ibodies[i][2] -= Fy;
      Ibodies[i][3] -= Fz;
    }
  }

  void M2LPeriodic() const {
    real_t M[MTERM];
    for_m M[m] = Multipole[0][m];
    real_t L[LTERM];
    for_l L[l] = 0;
    for (int lev=1; lev<numImages; lev++) {
      real_t diameter = 2 * R0 * pow(3,lev-1);
      int jx[3];
      for (jx[2]=-4; jx[2]<=4; jx[2]++) {
	for (jx[1]=-4; jx[1]<=4; jx[1]++) {
	  for (jx[0]=-4; jx[0]<=4; jx[0]++) {
	    if(jx[0] < -1 || 1 < jx[0] ||
	       jx[1] < -1 || 1 < jx[1] ||
	       jx[2] < -1 || 1 < jx[2]) {
	      real_t dX[3];
	      for_3 dX[d] = jx[d] * diameter;
	      real_t invR2 = 1. / (dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2]);
	      real_t invR  = sqrt(invR2);
	      real_t C[LTERM];
	      getCoef(C,dX,invR2,invR);
	      M2LSum(L,C,M);
	    }
	  }
	}
      }
      real_t M3[MTERM];
      for_m M3[m] = 0;
      int ix[3];
      for( ix[2]=-1; ix[2]<=1; ix[2]++ ) {
	for( ix[1]=-1; ix[1]<=1; ix[1]++ ) {
	  for( ix[0]=-1; ix[0]<=1; ix[0]++ ) {
	    real_t dX[3];
	    for_3 dX[d] = ix[d] * diameter;
	    real_t C[LTERM];
	    C[0] = 1;
	    powerM(C,dX);
	    for_m M3[m] += C[m] * M[0];
	    M2MSum(M3,C,M);
	  }
	}
      }
      for_m M[m] = M3[m];
    }
    for_l Local[0][l] += L[l];
  }
  
protected:
  inline int getKey(int *ix, int level, bool levelOffset=true) const {
    int id = 0;
    if (levelOffset) id = ((1 << 3 * level) - 1) / 7;
    for (int lev=0; lev<level; lev++) {
      for_3 id += ((ix[d] >> lev) & 1) << (3 * lev + d);
    }
    return id;
  }

public:
#ifndef SAKURA
  void P2P() const {
    int nunit = 1 << maxLevel;
    int nmin = 0;
    int nmax = nunit - 1;
    if (numImages != 0) {
      nmin -= nunit;
      nmax += nunit;
    }
#pragma omp parallel for
    for (int i=0; i<numLeafs; i++) {
      int ix[3] = {0, 0, 0};
      getIndex(ix,i);
      int jxmin[3];
      for_3 jxmin[d] = MAX(nmin, ix[d] - numNeighbors);
      int jxmax[3];
      for_3 jxmax[d] = MIN(nmax, ix[d] + numNeighbors);
      int jx[3];
      for (jx[2]=jxmin[2]; jx[2]<=jxmax[2]; jx[2]++) {
        for (jx[1]=jxmin[1]; jx[1]<=jxmax[1]; jx[1]++) {
          for (jx[0]=jxmin[0]; jx[0]<=jxmax[0]; jx[0]++) {
	    int jxp[3];
	    for_3 jxp[d] = (jx[d] + nunit) % nunit;
            int j = getKey(jxp,maxLevel,false);
	    real_t Xperiodic[3] = {0, 0, 0};
	    for_3 jxp[d] = (jx[d] + nunit) / nunit;
	    for_3 Xperiodic[d] = (jxp[d] - 1) * 2 * R0;
            P2PSum(Leafs[i][0],Leafs[i][1],Leafs[j][0],Leafs[j][1],Xperiodic);
          }
        }
      }
    }
  }
#else
  void P2P(int *nn_stencil) const {
    int nunit = 1 << maxLevel;
#pragma omp parallel for
    for( int i=0; i<numLeafs; i++ ) {
      int ix[3] = {0, 0, 0};
      getIndex(ix,i);
      int jx[3];
      int node_code[DIM];
      node_code[0] = ix[2]; node_code[1] = ix[1]; node_code[2] = ix[0];
      int nn_codes[DIM*NN];
      int nn_count = executeStencilNN(nn_codes, node_code, nn_stencil, 1, maxLevel);
      for(int k=0; k<nn_count; k++){
        jx[2] = nn_codes[k*DIM + 0];
        jx[1] = nn_codes[k*DIM + 1];
        jx[0] = nn_codes[k*DIM + 2];
	int jxp[3];
	for_3 jxp[d] = (jx[d] + nunit) % nunit;
	int j = getKey(jxp,maxLevel,false);
	real_t Xperiodic[3] = {0, 0, 0};
	for_3 jxp[d] = (jx[d] + nunit) / nunit;
	for_3 Xperiodic[d] = (jxp[d] - 1) * 2 * R0;
	P2PSum(Leafs[i][0],Leafs[i][1],Leafs[j][0],Leafs[j][1],Xperiodic);
      }
      jx[0] = ix[0]; jx[1] = ix[1]; jx[2] = ix[2];
      int jxp[3];
      for_3 jxp[d] = (jx[d] + nunit) % nunit;
      int j = getKey(jxp,maxLevel,false);
      real_t Xperiodic[3] = {0, 0, 0};
      for_3 jxp[d] = (jx[d] + nunit) / nunit;
      for_3 Xperiodic[d] = (jxp[d] - 1) * 2 * R0;
      P2PSum(Leafs[i][0],Leafs[i][1],Leafs[j][0],Leafs[j][1],Xperiodic);
    }
  }
#endif

  void P2M() const {
    int levelOffset = ((1 << 3 * maxLevel) - 1) / 7;
#pragma omp parallel for
    for (int i=0; i<numLeafs; i++) {
      real_t center[3];
      getCenter(center,i,maxLevel);
      for (int j=Leafs[i][0]; j<Leafs[i][1]; j++) {
        real_t dX[3];
        for_3 dX[d] = center[d] - Jbodies[j][d];
        real_t M[MTERM];
        M[0] = Jbodies[j][3];
        powerM(M,dX);
        for_m Multipole[i+levelOffset][m] += M[m];
      }
    }
  }

  void M2M() const {
    for (int lev=maxLevel; lev>0; lev--) {
      int childOffset = ((1 << 3 * lev) - 1) / 7;
      int parentOffset = ((1 << 3 * (lev - 1)) - 1) / 7;
      real_t radius = R0 / (1 << lev);
#pragma omp parallel for schedule(static, 8)
      for (int i=0; i<(1 << 3 * lev); i++) {
        int c = i + childOffset;
        int p = (i >> 3) + parentOffset;
        int ix[3];
        ix[0] = 1 - (i & 1) * 2;
        ix[1] = 1 - ((i / 2) & 1) * 2;
        ix[2] = 1 - ((i / 4) & 1) * 2;
        real_t dX[3];
        for_3 dX[d] = ix[d] * radius;
        real_t C[LTERM];
        C[0] = 1;
        powerM(C,dX);
        for_m Multipole[p][m] += C[m] * Multipole[c][0];
        M2MSum(Multipole[p],C,Multipole[c]);
      }
    }
  }

#ifndef SAKURA
  void M2L() const {
    for (int lev=1; lev<=maxLevel; lev++) {
      int levelOffset = ((1 << 3 * lev) - 1) / 7;
      int nunit = 1 << lev;
      real_t diameter = 2 * R0 / (1 << lev);
      int nmin = 0;
      int nmax = (nunit >> 1) - 1;
      if( numImages != 0 ) {
	nmin -= (nunit >> 1);
	nmax += (nunit >> 1);
      }
#pragma omp parallel for
      for (int i=0; i<(1 << 3 * lev); i++) {
        real_t L[LTERM];
        for_l L[l] = 0;
        int ix[3] = {0, 0, 0};
        getIndex(ix,i);
        int jxmin[3];
        for_3 jxmin[d] =  MAX(nmin, (ix[d] >> 1) - numNeighbors) << 1;
        int jxmax[3];
        for_3 jxmax[d] = (MIN(nmax, (ix[d] >> 1) + numNeighbors) << 1) + 1;
        int jx[3];
        for (jx[2]=jxmin[2]; jx[2]<=jxmax[2]; jx[2]++) {
          for (jx[1]=jxmin[1]; jx[1]<=jxmax[1]; jx[1]++) {
            for (jx[0]=jxmin[0]; jx[0]<=jxmax[0]; jx[0]++) {
              if(jx[0] < ix[0]-numNeighbors || ix[0]+numNeighbors < jx[0] ||
                 jx[1] < ix[1]-numNeighbors || ix[1]+numNeighbors < jx[1] ||
                 jx[2] < ix[2]-numNeighbors || ix[2]+numNeighbors < jx[2]) {
		int jxp[3];
		for_3 jxp[d] = (jx[d] + nunit) % nunit;
                int j = getKey(jxp,lev);
                real_t dX[3];
                for_3 dX[d] = (ix[d] - jx[d]) * diameter;
                real_t invR2 = 1. / (dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2]);
                real_t invR  = sqrt(invR2);
                real_t C[LTERM];
                getCoef(C,dX,invR2,invR);
                M2LSum(L,C,Multipole[j]);
              }
            }
          }
        }
        for_l Local[i+levelOffset][l] += L[l];
      }
    }
    if (numImages > 1) {
      M2LPeriodic();
    }
  }
#else
  void M2L(int *common_stencil, int *far_stencil) const {
    for( int lev=1; lev<=maxLevel; lev++ ) {
      int levelOffset = ((1 << 3 * lev) - 1) / 7;
      int nunit = 1 << lev;
      real_t diameter = 2 * R0 / (1 << lev);
#pragma omp parallel for
      for( int i=0; i<(1 << 3 * lev); i++ ) {
        real_t L[LTERM];
        for_l L[l] = 0;
        int ix[3] = {0, 0, 0};
        getIndex(ix,i);
        int jx[3];
        int node_code[DIM];
        node_code[0] = ix[2]; node_code[1] = ix[1]; node_code[2] = ix[0];
        int fn_codes[DIM*(FN+CN)];
        int fn_count = executeStencilFN(fn_codes, node_code, far_stencil, common_stencil, 1, lev);
        for(int k=0; k<fn_count; k++){
          jx[2] = fn_codes[k*DIM + 0];
          jx[1] = fn_codes[k*DIM + 1];
          jx[0] = fn_codes[k*DIM + 2];
	  int jxp[3];
	  for_3 jxp[d] = (jx[d] + nunit) % nunit;
	  int j = getKey(jxp,lev);
	  real_t dX[3];
	  for_3 dX[d] = (ix[d] - jx[d]) * diameter;
	  real_t invR2 = 1. / (dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2]);
	  real_t invR  = sqrt(invR2);
	  real_t C[LTERM];
	  getCoef(C,dX,invR2,invR);
	  M2LSum(L,C,Multipole[j]);
        }
        for_l Local[i+levelOffset][l] += L[l];
      }
    }
    if (numImages > 1) {
      M2LPeriodic();
    }
  }
#endif

  void L2L() const {
    for (int lev=1; lev<=maxLevel; lev++) {
      int childOffset = ((1 << 3 * lev) - 1) / 7;
      int parentOffset = ((1 << 3 * (lev - 1)) - 1) / 7;
      real_t radius = R0 / (1 << lev);
#pragma omp parallel for
      for (int i=0; i<(1 << 3 * lev); i++) {
        int c = i + childOffset;
        int p = (i >> 3) + parentOffset;
        int ix[3];
        ix[0] = (i & 1) * 2 - 1;
        ix[1] = ((i / 2) & 1) * 2 - 1;
        ix[2] = ((i / 4) & 1) * 2 - 1;
        real_t dX[3];
        for_3 dX[d] = ix[d] * radius;
        real_t C[LTERM];
        C[0] = 1;
        powerL(C,dX);
        for_l Local[c][l] += Local[p][l];
        for (int l=1; l<LTERM; l++) Local[c][0] += C[l] * Local[p][l];
        L2LSum(Local[c],C,Local[p]);
      }
    }
  }

  void L2P() const {
    int levelOffset = ((1 << 3 * maxLevel) - 1) / 7;
#pragma omp parallel for
    for (int i=0; i<numLeafs; i++) {
      real_t center[3];
      getCenter(center,i,maxLevel);
      real_t L[LTERM];
      for_l L[l] = Local[i+levelOffset][l];
      for (int j=Leafs[i][0]; j<Leafs[i][1]; j++) {
        real_t dX[3];
        for_3 dX[d] = Jbodies[j][d] - center[d];
        real_t C[LTERM];
        C[0] = 1;
        powerL(C,dX);
        for_4 Ibodies[j][d] += L[d];
        for (int l=1; l<LTERM; l++) Ibodies[j][0] += C[l] * L[l];
        L2PSum(Ibodies[j],C,L);
      }
    }
  }
};
