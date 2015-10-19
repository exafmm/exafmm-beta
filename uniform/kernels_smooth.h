#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <omp.h>

#define EXAFMM_PP 6
const int DP2P = 2; // Use 1 for parallel
const int DM2L = 2; // Use 1 for parallel
const int MTERM = EXAFMM_PP*(EXAFMM_PP+1)*(EXAFMM_PP+2)/6;
const int LTERM = (EXAFMM_PP+1)*(EXAFMM_PP+2)*(EXAFMM_PP+3)/6;
const real_t ALPHA_M = 100;
const real_t ALPHA_L = 100;

#include "core.h"

#define for_3d for (int d=0; d<3; d++)
#define for_4d for (int d=0; d<4; d++)
#define for_m for (int m=0; m<MTERM; m++)
#define for_l for (int l=0; l<LTERM; l++)
#define EXAFMM_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define EXAFMM_MIN(a,b) (((a) < (b)) ? (a) : (b))

namespace exafmm {
  class Kernel {
  public:
    int maxLevel;
    int maxGlobLevel;
    int numBodies;
    int numImages;
    int numCells;
    int numLeafs;
    int numGlobCells;
    int numPartition[10][3];
    int globLevelOffset[10];
    int numSendBodies;
    int numSendCells;
    int numSendLeafs;
    int MPISIZE;
    int MPIRANK;

    real_t X0[3];
    real_t R0;
    real_t RGlob[3];
    int *Index;
    int *Rank;
    int *sendIndex;
    int *recvIndex;
    int (*Leafs)[2];
    int (*sendLeafs)[2];
    int (*recvLeafs)[2];
    real_t (*Ibodies)[4];
    real_t (*Jbodies)[4];
    real_t (*Multipole)[MTERM];
    real_t (*Local)[LTERM];
    real_t (*globMultipole)[MTERM];
    real_t (*globLocal)[LTERM];
    float (*sendJbodies)[4];
    float (*recvJbodies)[4];
    float (*sendMultipole)[MTERM];
    float (*recvMultipole)[MTERM];

  private:
    inline void getIndex(int *ix, int index) const {
      for_3d ix[d] = 0;
      int d = 0, level = 0;
      while( index != 0 ) {
	ix[d] += (index % 2) * (1 << level);
	index >>= 1;
	d = (d+1) % 3;
	if( d == 0 ) level++;
      }
    }

    void getCenter(real_t *dX, int index, int level) const {
      real_t R = R0 / (1 << level);
      int ix[3] = {0, 0, 0};
      getIndex(ix, index);
      for_3d dX[d] = X0[d] - R0 + (2 * ix[d] + 1) * R;
    }

  protected:
    inline int getGlobKey(int *ix, int level) const {
      return ix[0] + (ix[1] + ix[2] * numPartition[level][1]) * numPartition[level][0];
    }

    void P2P(int ibegin, int iend, int jbegin, int jend,
	     real_t *Ximin, real_t *Ximax, real_t *Xjmin, real_t *Xjmax, real_t *periodic) const {
      for( int i=ibegin; i<iend; i++ ) {
	real_t wi = 1;
	for_3d wi *= 1 - erfc(ALPHA_L*(Jbodies[i][d] - Ximin[d])) / 2;
	for_3d wi *= erfc(ALPHA_L*(Jbodies[i][d] - Ximax[d])) / 2;
	real_t Po = 0, Fx = 0, Fy = 0, Fz = 0;
	for( int j=jbegin; j<jend; j++ ) {
	  real_t wj = 1;
	  for_3d wj *= 1 - erfc(ALPHA_M*(Jbodies[j][d] + periodic[d] - Xjmin[d])) / 2;
	  for_3d wj *= erfc(ALPHA_M*(Jbodies[j][d] + periodic[d] - Xjmax[d])) / 2;
	  real_t dX[3];
	  for_3d dX[d] = Jbodies[i][d] - Jbodies[j][d] - periodic[d];
	  real_t R2 = dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2];
	  real_t invR2 = 1.0 / R2;                                  
	  if( R2 == 0 ) invR2 = 0;                                
	  real_t invR = Jbodies[j][3] * sqrt(invR2) * wi * wj;
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

    void P2P() const {
      int ixc[3];
      getGlobIndex(ixc,MPIRANK,maxGlobLevel);
      int nunit = 1 << maxLevel;
      int nunitGlob[3];
      for_3d nunitGlob[d] = nunit * numPartition[maxGlobLevel][d];
      int nxmin[3], nxmax[3];
      for_3d nxmin[d] = -ixc[d] * nunit;
      for_3d nxmax[d] = nunitGlob[d] + nxmin[d] - 1;
      if( numImages != 0 ) {
	for_3d nxmin[d] -= nunitGlob[d];
	for_3d nxmax[d] += nunitGlob[d];
      }
      real_t diameter = 2 * R0 / (1 << maxLevel);
#pragma omp parallel for
      for( int i=0; i<numLeafs; i++ ) {
	int ix[3] = {0, 0, 0};
	getIndex(ix,i);
	real_t Ximin[3], Ximax[3];
	for_3d Ximin[d] = diameter * ix[d];
	for_3d Ximax[d] = diameter * (ix[d] + 1);
	int ixmin[3], ixmax[3];
	for_3d ixmin[d] = EXAFMM_MAX(nxmin[d],ix[d] - 1);
	for_3d ixmax[d] = EXAFMM_MIN(nxmax[d],ix[d] + 1);
	int jxmin[3], jxmax[3];
	for_3d jxmin[d] = EXAFMM_MAX(nxmin[d],ix[d] - DP2P);
	for_3d jxmax[d] = EXAFMM_MIN(nxmax[d],ix[d] + DP2P);
	real_t Xjmin[3], Xjmax[3];
	for_3d Xjmin[d] = diameter * jxmin[d];
	for_3d Xjmax[d] = diameter * (jxmax[d] + 1);
	for_3d jxmin[d] = EXAFMM_MAX(nxmin[d],ix[d] - DP2P - 1);
	for_3d jxmax[d] = EXAFMM_MIN(nxmax[d],ix[d] + DP2P + 1);
	for( ix[2]=ixmin[2]; ix[2]<=ixmax[2]; ix[2]++ ) {
	  for( ix[1]=ixmin[1]; ix[1]<=ixmax[1]; ix[1]++ ) {
	    for( ix[0]=ixmin[0]; ix[0]<=ixmax[0]; ix[0]++ ) {
	      int ixp[3];
	      for_3d ixp[d] = (ix[d] + nunit) % nunit;
	      int ii = getKey(ixp,maxLevel,false);
	      int jx[3];
	      for( jx[2]=jxmin[2]; jx[2]<=jxmax[2]; jx[2]++ ) {
		for( jx[1]=jxmin[1]; jx[1]<=jxmax[1]; jx[1]++ ) {
		  for( jx[0]=jxmin[0]; jx[0]<=jxmax[0]; jx[0]++ ) {
		    int jxp[3];
		    for_3d jxp[d] = (jx[d] + nunit) % nunit;
		    int j = getKey(jxp,maxLevel,false);
		    for_3d jxp[d] = (jx[d] + nunit) / nunit;
#if EXAFMM_SERIAL
		    int rankOffset = 13 * numLeafs;
#else
		    int rankOffset = (jxp[0] + 3 * jxp[1] + 9 * jxp[2]) * numLeafs;
#endif
		    j += rankOffset;
		    rankOffset = 13 * numLeafs;
		    real_t periodic[3] = {0, 0, 0};
		    for_3d jxp[d] = (jx[d] + ixc[d] * nunit + nunitGlob[d]) / nunitGlob[d];
		    for_3d periodic[d] = (jxp[d] - 1) * 2 * RGlob[d];
		    P2P(Leafs[ii+rankOffset][0],Leafs[ii+rankOffset][1],Leafs[j][0],Leafs[j][1],
			Ximin,Ximax,Xjmin,Xjmax,periodic);
		  }
		}
	      }
	    }
	  }
	}
      }
    }

    void P2M() const {
      int ixc[3];
      getGlobIndex(ixc,MPIRANK,maxGlobLevel);
      int nunit = 1 << maxLevel;
      int nunitGlob[3];
      for_3d nunitGlob[d] = nunit * numPartition[maxGlobLevel][d];
      int nxmin[3], nxmax[3];
      for_3d nxmin[d] = -ixc[d] * nunit;
      for_3d nxmax[d] = nunitGlob[d] + nxmin[d] - 1;
      if( numImages != 0 ) {
	for_3d nxmin[d] -= nunitGlob[d];
	for_3d nxmax[d] += nunitGlob[d];
      }
      real_t diameter = 2 * R0 / (1 << maxLevel);
      int levelOffset = ((1 << 3 * maxLevel) - 1) / 7 + 13 * numCells;
#pragma omp parallel for
      for( int i=0; i<numLeafs; i++ ) {
	real_t center[3];
	getCenter(center,i,maxLevel);
	int jx[3] = {0, 0, 0};
	getIndex(jx,i);
	real_t Xjmin[3], Xjmax[3];
	for_3d Xjmin[d] = diameter * jx[d];
	for_3d Xjmax[d] = diameter * (jx[d] + 1);
	int jxmin[3], jxmax[3];
	for_3d jxmin[d] = EXAFMM_MAX(nxmin[d],jx[d] - 1);
	for_3d jxmax[d] = EXAFMM_MIN(nxmax[d],jx[d] + 1);
	for( jx[2]=jxmin[2]; jx[2]<=jxmax[2]; jx[2]++ ) {
	  for( jx[1]=jxmin[1]; jx[1]<=jxmax[1]; jx[1]++ ) {
	    for( jx[0]=jxmin[0]; jx[0]<=jxmax[0]; jx[0]++ ) {
	      int jxp[3];
	      for_3d jxp[d] = (jx[d] + nunit) % nunit;
	      int jj = getKey(jxp,maxLevel,false);
	      for_3d jxp[d] = (jx[d] + nunit) / nunit;
#if EXAFMM_SERIAL
	      int rankOffset = 13 * numLeafs;
#else
	      int rankOffset = (jxp[0] + 3 * jxp[1] + 9 * jxp[2]) * numLeafs;
#endif
	      real_t periodic[3] = {0, 0, 0};
	      for_3d jxp[d] = (jx[d] + ixc[d] * nunit + nunitGlob[d]) / nunitGlob[d];
	      for_3d periodic[d] = (jxp[d] - 1) * 2 * RGlob[d];
	      for( int j=Leafs[jj+rankOffset][0]; j<Leafs[jj+rankOffset][1]; j++ ) {
		real_t weight = 1;
		for_3d weight *= 1 - erfc(ALPHA_M*(Jbodies[j][d] + periodic[d] - Xjmin[d])) / 2;
		for_3d weight *= erfc(ALPHA_M*(Jbodies[j][d] + periodic[d] - Xjmax[d])) / 2;
		real_t dX[3];
		for_3d dX[d] = center[d] - Jbodies[j][d] - periodic[d];
		real_t M[MTERM];
		M[0] = Jbodies[j][3] * weight;
		powerM(M,dX);
		for_m Multipole[i+levelOffset][m] += M[m];
	      }
	    }
	  }
	}
      }
    }

    void M2M() const {
      int rankOffset = 13 * numCells;
      for( int lev=maxLevel; lev>0; lev-- ) {
	int childOffset = ((1 << 3 * lev) - 1) / 7 + rankOffset;
	int parentOffset = ((1 << 3 * (lev - 1)) - 1) / 7 + rankOffset;
	real_t radius = R0 / (1 << lev);
#pragma omp parallel for schedule(static, 8)
	for( int i=0; i<(1 << 3 * lev); i++ ) {
	  int c = i + childOffset;
	  int p = (i >> 3) + parentOffset;
	  int ix[3];
	  ix[0] = 1 - (i & 1) * 2;
	  ix[1] = 1 - ((i / 2) & 1) * 2;
	  ix[2] = 1 - ((i / 4) & 1) * 2;
	  real_t dX[3];
	  for_3d dX[d] = ix[d] * radius;
	  real_t M[MTERM];
	  real_t C[LTERM];
	  C[0] = 1;
	  powerM(C,dX);
	  for_m M[m] = Multipole[c][m];
	  for_m Multipole[p][m] += C[m] * M[0];
	  M2MSum(Multipole[p],C,M);
	}
      }
    }

    void M2L() const {
      int ixc[3];
      int DM2LC = DM2L;
      getGlobIndex(ixc,MPIRANK,maxGlobLevel);
      for( int lev=1; lev<=maxLevel; lev++ ) {
	if (lev==maxLevel) DM2LC = DP2P;
	int levelOffset = ((1 << 3 * lev) - 1) / 7;
	int nunit = 1 << lev;
	int nunitGlob[3];
	for_3d nunitGlob[d] = nunit * numPartition[maxGlobLevel][d];
	int nxmin[3], nxmax[3];
	for_3d nxmin[d] = -ixc[d] * (nunit >> 1);
	for_3d nxmax[d] = (nunitGlob[d] >> 1) + nxmin[d] - 1;
	if( numImages != 0 ) {
	  for_3d nxmin[d] -= (nunitGlob[d] >> 1);
	  for_3d nxmax[d] += (nunitGlob[d] >> 1);
	}
	real_t diameter = 2 * R0 / (1 << lev);
#pragma omp parallel for
	for( int i=0; i<(1 << 3 * lev); i++ ) {
	  real_t L[LTERM];
	  for_l L[l] = 0;
	  int ix[3] = {0, 0, 0};
	  getIndex(ix,i);
	  int jxmin[3];
	  for_3d jxmin[d] = (EXAFMM_MAX(nxmin[d],(ix[d] >> 1) - DM2L) << 1);
	  int jxmax[3];
	  for_3d jxmax[d] = (EXAFMM_MIN(nxmax[d],(ix[d] >> 1) + DM2L) << 1) + 1;
	  int jx[3];
	  for( jx[2]=jxmin[2]; jx[2]<=jxmax[2]; jx[2]++ ) {
	    for( jx[1]=jxmin[1]; jx[1]<=jxmax[1]; jx[1]++ ) {
	      for( jx[0]=jxmin[0]; jx[0]<=jxmax[0]; jx[0]++ ) {
		if(jx[0] < ix[0]-DM2LC || ix[0]+DM2LC < jx[0] ||
		   jx[1] < ix[1]-DM2LC || ix[1]+DM2LC < jx[1] ||
		   jx[2] < ix[2]-DM2LC || ix[2]+DM2LC < jx[2]) {
		  int jxp[3];
		  for_3d jxp[d] = (jx[d] + nunit) % nunit;
		  int j = getKey(jxp,lev);
		  for_3d jxp[d] = (jx[d] + nunit) / nunit;
#if EXAFMM_SERIAL
		  int rankOffset = 13 * numCells;
#else
		  int rankOffset = (jxp[0] + 3 * jxp[1] + 9 * jxp[2]) * numCells;
#endif
		  j += rankOffset;
		  real_t M[MTERM];
		  for_m M[m] = Multipole[j][m];
		  real_t dX[3];
		  for_3d dX[d] = (ix[d] - jx[d]) * diameter;
		  real_t invR2 = 1. / (dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2]);
		  real_t invR  = sqrt(invR2);
		  real_t C[LTERM];
		  getCoef(C,dX,invR2,invR);
		  M2LSum(L,C,M);
		}
	      }
	    }
	  }
	  for_l Local[i+levelOffset][l] += L[l];
	}
      }
    }

    void L2L() const {
      for( int lev=1; lev<=maxLevel; lev++ ) {
	int childOffset = ((1 << 3 * lev) - 1) / 7;
	int parentOffset = ((1 << 3 * (lev - 1)) - 1) / 7;
	real_t radius = R0 / (1 << lev);
#pragma omp parallel for
	for( int i=0; i<(1 << 3 * lev); i++ ) {
	  int c = i + childOffset;
	  int p = (i >> 3) + parentOffset;
	  int ix[3];
	  ix[0] = (i & 1) * 2 - 1;
	  ix[1] = ((i / 2) & 1) * 2 - 1;
	  ix[2] = ((i / 4) & 1) * 2 - 1;
	  real_t dX[3];
	  for_3d dX[d] = ix[d] * radius;
	  real_t C[LTERM];
	  C[0] = 1;
	  powerL(C,dX);
	  for_l Local[c][l] += Local[p][l];
	  for( int l=1; l<LTERM; l++ ) Local[c][0] += C[l] * Local[p][l];
	  L2LSum(Local[c],C,Local[p]);
	}
      }
    }

    void L2P() const {
      int ixc[3];
      getGlobIndex(ixc,MPIRANK,maxGlobLevel);
      int nunit = 1 << maxLevel;
      int nunitGlob[3];
      for_3d nunitGlob[d] = nunit * numPartition[maxGlobLevel][d];
      int nxmin[3], nxmax[3];
      for_3d nxmin[d] = -ixc[d] * nunit;
      for_3d nxmax[d] = nunitGlob[d] + nxmin[d] - 1;
      if( numImages != 0 ) {
	for_3d nxmin[d] -= nunitGlob[d];
	for_3d nxmax[d] += nunitGlob[d];
      }
      real_t diameter = 2 * R0 / (1 << maxLevel);
      int levelOffset = ((1 << 3 * maxLevel) - 1) / 7;
#pragma omp parallel for
      for( int i=0; i<numLeafs; i++ ) {
	real_t center[3];
	getCenter(center,i,maxLevel);
	int ix[3] = {0, 0, 0};
	getIndex(ix,i);
	real_t Ximin[3], Ximax[3];
	for_3d Ximin[d] = diameter * ix[d];
	for_3d Ximax[d] = diameter * (ix[d] + 1);
	int ixmin[3], ixmax[3];
	for_3d ixmin[d] = EXAFMM_MAX(nxmin[d],ix[d] - 1);
	for_3d ixmax[d] = EXAFMM_MIN(nxmax[d],ix[d] + 1);
	for( ix[2]=ixmin[2]; ix[2]<=ixmax[2]; ix[2]++ ) {
	  for( ix[1]=ixmin[1]; ix[1]<=ixmax[1]; ix[1]++ ) {
	    for( ix[0]=ixmin[0]; ix[0]<=ixmax[0]; ix[0]++ ) {
	      int ixp[3];
	      for_3d ixp[d] = (ix[d] + nunit) % nunit;
	      int ii = getKey(ixp,maxLevel,false);
	      for_3d ixp[d] = (ix[d] + nunit) / nunit;
#if EXAFMM_SERIAL
	      int rankOffset = 13 * numLeafs;
#else
	      int rankOffset = (ixp[0] + 3 * ixp[1] + 9 * ixp[2]) * numLeafs;
#endif
	      real_t periodic[3] = {0, 0, 0};
	      for_3d ixp[d] = (ix[d] + ixc[d] * nunit + nunitGlob[d]) / nunitGlob[d];
	      for_3d periodic[d] = (ixp[d] - 1) * 2 * RGlob[d];
	      real_t L[LTERM];
	      for_l L[l] = Local[i+levelOffset][l];
	      for( int j=Leafs[ii+rankOffset][0]; j<Leafs[ii+rankOffset][1]; j++ ) {
		real_t weight = 1;
		for_3d weight *= 1 - erfc(ALPHA_L*(Jbodies[j][d] + periodic[d] - Ximin[d])) / 2;
		for_3d weight *= erfc(ALPHA_L*(Jbodies[j][d] + periodic[d] - Ximax[d])) / 2;
		real_t dX[3];
		for_3d dX[d] = Jbodies[j][d] + periodic[d] - center[d];
		real_t C[LTERM];
		C[0] = weight;
		powerL(C,dX);
		for_4d Ibodies[j][d] += L[d] * weight;
		for( int l=1; l<LTERM; l++ ) Ibodies[j][0] += C[l] * L[l];
		L2PSum(Ibodies[j],C,L);
	      }
	    }
	  }
	}
      }
    }

  public:
    Kernel() : MPISIZE(1), MPIRANK(0) {}
    ~Kernel() {}

    inline int getKey(int *ix, int level, bool levelOffset=true) const {
      int id = 0;
      if( levelOffset ) id = ((1 << 3 * level) - 1) / 7;
      for( int lev=0; lev<level; ++lev ) {
	for_3d id += ix[d] % 2 << (3 * lev + d);
	for_3d ix[d] >>= 1;
      }
      return id;
    }

    inline void getGlobIndex(int *ix, int index, int level) const {
      ix[0] = index % numPartition[level][0];
      ix[1] = index / numPartition[level][0] % numPartition[level][1];
      ix[2] = index / numPartition[level][0] / numPartition[level][1];
    }
  };
}
