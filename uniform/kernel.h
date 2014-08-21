#include "types.h"
#include "core.h"

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

  real X0[3];
  real R0;
  real RGlob[3];
  real (*Ibodies)[4];
  real (*Jbodies)[4];
  real (*Multipole)[MTERM];
  real (*Local)[LTERM];
  real (*globMultipole)[MTERM];
  real (*globLocal)[LTERM];
  int (*Leafs)[2];
  float (*sendJbodies)[4];
  float (*recvJbodies)[4];
  float (*sendMultipole)[MTERM];
  float (*recvMultipole)[MTERM];
  int (*sendLeafs)[2];
  int (*recvLeafs)[2];

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

  void getCenter(real *dist, int index, int level) const {
    real R = R0 / (1 << level);
    int ix[3] = {0, 0, 0};
    getIndex(ix, index);
    for_3d dist[d] = X0[d] - R0 + (2 * ix[d] + 1) * R;
  }

protected:
  inline int getGlobKey(int *ix, int level) const {
    return ix[0] + (ix[1] + ix[2] * numPartition[level][1]) * numPartition[level][0];
  }

  void P2P(int ibegin, int iend, int jbegin, int jend, real *periodic) const {
    int ii;
    for( ii=ibegin; ii<iend-1; ii+=2 ) {
#ifndef SPARC_SIMD
      for( int i=ii; i<=ii+1; i++ ) {
        real Po = 0, Fx = 0, Fy = 0, Fz = 0;
        for( int j=jbegin; j<jend; j++ ) {
          real dist[3];
          for_3d dist[d] = Jbodies[i][d] - Jbodies[j][d] - periodic[d];
          real R2 = dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];
          real invR2 = 1.0 / R2;                                  
          if( R2 == 0 ) invR2 = 0;                                
          real invR = Jbodies[j][3] * sqrt(invR2);
          real invR3 = invR2 * invR;
          Po += invR;
          Fx += dist[0] * invR3;
          Fy += dist[1] * invR3;
          Fz += dist[2] * invR3;
        }
        Ibodies[i][0] += Po;
        Ibodies[i][1] -= Fx;
        Ibodies[i][2] -= Fy;
        Ibodies[i][3] -= Fz;
      }
#else
      __m128d Po = _mm_setzero_pd();
      __m128d Fx = _mm_setzero_pd();
      __m128d Fy = _mm_setzero_pd();
      __m128d Fz = _mm_setzero_pd();
      __m128d zero = _mm_setzero_pd();
      __m128d xi[3], xj[4];
      for_3d xi[d] = _mm_set_pd(Jbodies[ii][d]-periodic[d],Jbodies[ii+1][d]-periodic[d]);
      for( int j=jbegin; j<jend; j++ ) {
        for_4d xj[d] = _mm_set_pd(Jbodies[j][d],Jbodies[j][d]);
        for_3d xj[d] = _mm_sub_pd(xi[d],xj[d]);
        __m128d R2 = _mm_mul_pd(xj[0],xj[0]);
        R2 = _fjsp_madd_v2r8(xj[1],xj[1],R2);
        R2 = _fjsp_madd_v2r8(xj[2],xj[2],R2);
        __m128d invR = _fjsp_rsqrta_v2r8(R2);
        R2 = _mm_cmpneq_pd(R2,zero);
        invR = _mm_and_pd(invR,R2);
        R2 = _mm_mul_pd(invR,invR);
        invR = _mm_mul_pd(invR,xj[3]);
        R2 = _mm_mul_pd(R2,invR);
        Po = _mm_add_pd(Po,invR);
        Fx = _fjsp_madd_v2r8(xj[0],R2,Fx);
        Fy = _fjsp_madd_v2r8(xj[1],R2,Fy);
        Fz = _fjsp_madd_v2r8(xj[2],R2,Fz);
      }
      real Po2[2], Fx2[2], Fy2[2], Fz2[2];
      _mm_store_pd(Po2,Po);
      _mm_store_pd(Fx2,Fx);
      _mm_store_pd(Fy2,Fy);
      _mm_store_pd(Fz2,Fz);
      Ibodies[ii][0] += Po2[1];
      Ibodies[ii][1] -= Fx2[1];
      Ibodies[ii][2] -= Fy2[1];
      Ibodies[ii][3] -= Fz2[1];
      Ibodies[ii+1][0] += Po2[0];
      Ibodies[ii+1][1] -= Fx2[0];
      Ibodies[ii+1][2] -= Fy2[0];
      Ibodies[ii+1][3] -= Fz2[0];
#endif
    }
    for( int i=ii; i<iend; i++ ) {
      real Po = 0, Fx = 0, Fy = 0, Fz = 0;
      for( int j=jbegin; j<jend; j++ ) {
        real dist[3];
        for_3d dist[d] = Jbodies[i][d] - Jbodies[j][d] - periodic[d];
        real R2 = dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];
        real invR2 = 1.0 / R2;
        if( R2 == 0 ) invR2 = 0;
        real invR = Jbodies[j][3] * sqrt(invR2);
        real invR3 = invR2 * invR;
        Po += invR;
        Fx += dist[0] * invR3;
        Fy += dist[1] * invR3;
        Fz += dist[2] * invR3;
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
#pragma omp parallel for
    for( int i=0; i<numLeafs; i++ ) {
      int ix[3] = {0, 0, 0};
      getIndex(ix,i);
      int jxmin[3];
      for_3d jxmin[d] = FMMMAX(nxmin[d],ix[d] - 1);
      int jxmax[3];
      for_3d jxmax[d] = FMMMIN(nxmax[d],ix[d] + 1) + 1;
      int jx[3];
      for( jx[2]=jxmin[2]; jx[2]<jxmax[2]; jx[2]++ ) {
        for( jx[1]=jxmin[1]; jx[1]<jxmax[1]; jx[1]++ ) {
          for( jx[0]=jxmin[0]; jx[0]<jxmax[0]; jx[0]++ ) {
            int jxp[3];
            for_3d jxp[d] = (jx[d] + nunit) % nunit;
            int j = getKey(jxp,maxLevel,false);
            for_3d jxp[d] = (jx[d] + nunit) / nunit;
#if Serial
            int rankOffset = 13 * numLeafs;
#else
            int rankOffset = (jxp[0] + 3 * jxp[1] + 9 * jxp[2]) * numLeafs;
#endif
            j += rankOffset;
            rankOffset = 13 * numLeafs;
            real periodic[3] = {0, 0, 0};
            for_3d jxp[d] = (jx[d] + ixc[d] * nunit + nunitGlob[d]) / nunitGlob[d];
            for_3d periodic[d] = (jxp[d] - 1) * 2 * RGlob[d];
            P2P(Leafs[i+rankOffset][0],Leafs[i+rankOffset][1],Leafs[j][0],Leafs[j][1],periodic);
          }
        }
      }
    }
  }

  void P2M() const {
    int rankOffset = 13 * numLeafs;
    int levelOffset = ((1 << 3 * maxLevel) - 1) / 7 + 13 * numCells;
#pragma omp parallel for
    for( int i=0; i<numLeafs; i++ ) {
      real center[3];
      getCenter(center,i,maxLevel);
      for( int j=Leafs[i+rankOffset][0]; j<Leafs[i+rankOffset][1]; j++ ) {
        real dist[3];
        for_3d dist[d] = center[d] - Jbodies[j][d];
        real M[MTERM];
        M[0] = Jbodies[j][3];
        powerM(M,dist);
        for_m Multipole[i+levelOffset][m] += M[m];
      }
    }
  }

  void M2M() const {
    int rankOffset = 13 * numCells;
    for( int lev=maxLevel; lev>0; lev-- ) {
      int childOffset = ((1 << 3 * lev) - 1) / 7 + rankOffset;
      int parentOffset = ((1 << 3 * (lev - 1)) - 1) / 7 + rankOffset;
      real radius = R0 / (1 << lev);
#pragma omp parallel for schedule(static, 8)
      for( int i=0; i<(1 << 3 * lev); i++ ) {
        int c = i + childOffset;
        int p = (i >> 3) + parentOffset;
        int ix[3];
        ix[0] = 1 - (i & 1) * 2;
        ix[1] = 1 - ((i / 2) & 1) * 2;
        ix[2] = 1 - ((i / 4) & 1) * 2;
        real dist[3];
        for_3d dist[d] = ix[d] * radius;
        real M[MTERM];
        real C[LTERM];
        C[0] = 1;
        powerM(C,dist);
        for_m M[m] = Multipole[c][m];
        for_m Multipole[p][m] += C[m] * M[0];
        M2MSum(Multipole[p],C,M);
      }
    }
  }

  void M2L() const {
    int ixc[3];
    getGlobIndex(ixc,MPIRANK,maxGlobLevel);
    for( int lev=1; lev<=maxLevel; lev++ ) {
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
      real diameter = 2 * R0 / (1 << lev);
#pragma omp parallel for
      for( int i=0; i<(1 << 3 * lev); i++ ) {
        real L[LTERM];
        for_l L[l] = 0;
        int ix[3] = {0, 0, 0};
        getIndex(ix,i);
        int jxmin[3];
        for_3d jxmin[d] =  FMMMAX(nxmin[d],(ix[d] >> 1) - 1)      << 1;
        int jxmax[3];
        for_3d jxmax[d] = (FMMMIN(nxmax[d],(ix[d] >> 1) + 1) + 1) << 1;
        int jx[3];
        for( jx[2]=jxmin[2]; jx[2]<jxmax[2]; jx[2]++ ) {
          for( jx[1]=jxmin[1]; jx[1]<jxmax[1]; jx[1]++ ) {
            for( jx[0]=jxmin[0]; jx[0]<jxmax[0]; jx[0]++ ) {
              if(jx[0] < ix[0]-1 || ix[0]+1 < jx[0] ||
                 jx[1] < ix[1]-1 || ix[1]+1 < jx[1] ||
                 jx[2] < ix[2]-1 || ix[2]+1 < jx[2]) {
                int jxp[3];
                for_3d jxp[d] = (jx[d] + nunit) % nunit;
                int j = getKey(jxp,lev);
                for_3d jxp[d] = (jx[d] + nunit) / nunit;
#if Serial
                int rankOffset = 13 * numCells;
#else
                int rankOffset = (jxp[0] + 3 * jxp[1] + 9 * jxp[2]) * numCells;
#endif
                j += rankOffset;
                real M[MTERM];
                for_m M[m] = Multipole[j][m];
                real dist[3];
                for_3d dist[d] = (ix[d] - jx[d]) * diameter;
                real invR2 = 1. / (dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2]);
                real invR  = sqrt(invR2);
                real C[LTERM];
                getCoef(C,dist,invR2,invR);
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
      real radius = R0 / (1 << lev);
#pragma omp parallel for
      for( int i=0; i<(1 << 3 * lev); i++ ) {
        int c = i + childOffset;
        int p = (i >> 3) + parentOffset;
        int ix[3];
        ix[0] = (i & 1) * 2 - 1;
        ix[1] = ((i / 2) & 1) * 2 - 1;
        ix[2] = ((i / 4) & 1) * 2 - 1;
        real dist[3];
        for_3d dist[d] = ix[d] * radius;
        real C[LTERM];
        C[0] = 1;
        powerL(C,dist);
        for_l Local[c][l] += Local[p][l];
        for( int l=1; l<LTERM; l++ ) Local[c][0] += C[l] * Local[p][l];
        L2LSum(Local[c],C,Local[p]);
      }
    }
  }

  void L2P() const {
    int rankOffset = 13 * numLeafs;
    int levelOffset = ((1 << 3 * maxLevel) - 1) / 7;
#pragma omp parallel for
    for( int i=0; i<numLeafs; i++ ) {
      real center[3];
      getCenter(center,i,maxLevel);
      real L[LTERM];
      for_l L[l] = Local[i+levelOffset][l];
      for( int j=Leafs[i+rankOffset][0]; j<Leafs[i+rankOffset][1]; j++ ) {
        real dist[3];
        for_3d dist[d] = Jbodies[j][d] - center[d];
        real C[LTERM];
        C[0] = 1;
        powerL(C,dist);
        for_4d Ibodies[j][d] += L[d];
        for( int l=1; l<LTERM; l++ ) Ibodies[j][0] += C[l] * L[l];
        L2PSum(Ibodies[j],C,L);
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
