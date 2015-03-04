#include <cassert>
#include <cmath>
#include "types2.h"

class Ewald {
private:
  int numBodies;
  int maxLevel;
  int numLeafs;
  int numWaves;
  int ksize;
  real_t alpha;
  real_t sigma;
  real_t cutoff;
  real_t scale;
  real_t R0;
  real_t X0[3];
  real_t *waveRe;
  real_t *waveIm;
  real_t (*waveK)[3];

private:
  inline int getKey(int *ix, int level) const {
    int id = 0;
    for (int lev=0; lev<level; lev++) {
      for_3 id += ((ix[d] >> lev) & 1) << (3 * lev + d);
    }
    return id;
  }

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

  void initWaves(real_t *waveRe, real_t *waveIm, real_t (*waveK)[3]) {
    numWaves = 0;
    for (int l=0; l<=ksize; l++) {
      int mmin = -ksize;
      if (l==0) mmin = 0;
      for (int m=mmin; m<=ksize; m++) {
	int nmin = -ksize;
	if (l==0 && m==0) nmin=1;
	for (int n=nmin; n<=ksize; n++) {
	  real_t k2 = l * l + m * m + n * n;
	  if (k2 <= ksize * ksize) {
	    waveK[numWaves][0] = l;
	    waveK[numWaves][1] = m;
	    waveK[numWaves][2] = n;
	    waveRe[numWaves] = waveIm[numWaves] = 0;
	    numWaves++;
	  }
	}
      }
    }
    assert(numWaves < 4. / 3 * M_PI * ksize * ksize * ksize);
  }

  void dft(real_t (*Jbodies)[4]) {
#pragma omp parallel for
    for (int i=0; i<numWaves; i++) { 
      waveRe[i] = waveIm[i] = 0;
      for (int j=0; j<numBodies; j++) {
	real_t th = 0;
	for_3 th += waveK[i][d] * Jbodies[j][d] * scale;
	waveRe[i] += Jbodies[j][3] * cos(th);
	waveIm[i] += Jbodies[j][3] * sin(th);
      }
    }
  }

  void idft(real_t (*Ibodies)[4], real_t (*Jbodies)[4]) {
#pragma omp parallel for
    for (int i=0; i<numBodies; i++) {
      for (int j=0; j<numWaves; j++) {
	real_t th = 0;
	for_3 th += waveK[j][d] * Jbodies[i][d] * scale;
	real_t dtmp = waveRe[j] * sin(th) - waveIm[j] * cos(th);
	Ibodies[i][0] += waveRe[j] * cos(th) + waveIm[j] * sin(th);
	for_3 Ibodies[i][d+1] -= dtmp * waveK[j][d] * scale;
      }
    }
  }

  void P2PEwald(int ibegin, int iend, int jbegin, int jend, real_t *Xperiodic,
		real_t (*Ibodies)[4], real_t (*Jbodies)[4]) const {
    for (int i=ibegin; i<iend; i++) {
      real_t Po = 0, Fx = 0, Fy = 0, Fz = 0;
      for (int j=jbegin; j<jend; j++) {
	real_t dX[3];
	for_3 dX[d] = Jbodies[i][d] - Jbodies[j][d] - Xperiodic[d];
	real_t R2 = dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2];
	if (0 < R2 && R2 < cutoff * cutoff) {
	  real_t R2s = R2 * alpha * alpha;
	  real_t Rs = sqrtf(R2s);
	  real_t invRs = 1 / Rs;
	  real_t invR2s = invRs * invRs;
	  real_t invR3s = invR2s * invRs;
	  real_t dtmp = Jbodies[j][3] * (M_2_SQRTPI * exp(-R2s) * invR2s + erfc(Rs) * invR3s);
	  dtmp *= alpha * alpha * alpha;
	  Po += Jbodies[j][3] * erfc(Rs) * invRs * alpha;
	  Fx += dX[0] * dtmp;
	  Fy += dX[1] * dtmp;
	  Fz += dX[2] * dtmp;
	}
      }
      Ibodies[i][0] += Po;
      Ibodies[i][1] -= Fx;
      Ibodies[i][2] -= Fy;
      Ibodies[i][3] -= Fz;
    }
  }
  
public:
  Ewald(int _numBodies, int _maxLevel, real_t cycle,
	int _ksize, real_t _alpha, real_t _sigma, real_t _cutoff) :
    ksize(_ksize), alpha(_alpha), sigma(_sigma), cutoff(_cutoff) {
    numBodies = _numBodies;
    maxLevel = _maxLevel;
    numLeafs = 1 << 3 * maxLevel;
    numWaves = 4. / 3 * M_PI * ksize * ksize * ksize;
    scale = 2 * M_PI / cycle;
    R0 = cycle * .5;
    for_3 X0[d] = R0;
    waveRe = new real_t [numWaves];
    waveIm = new real_t [numWaves];
    waveK = new real_t [numWaves][3]();
  }
  
  ~Ewald() {
    delete[] waveRe;
    delete[] waveIm;
    delete[] waveK;
  }
  
  void wavePart(real_t (*Ibodies2)[4], real_t (*Jbodies)[4]) {
    initWaves(waveRe, waveIm, waveK);
    dft(Jbodies);
    const real_t coef = .25 / M_PI / M_PI / sigma / R0;
    const real_t coef2 = scale * scale / (4 * alpha * alpha);
#pragma omp parallel for
    for (int w=0; w<numWaves; w++) {
      real_t k2 = 0;
      for_3 k2 += waveK[w][d] * waveK[w][d];
      real_t factor = coef * exp(-k2 * coef2) / k2;
      waveRe[w] *= factor;
      waveIm[w] *= factor;
    }
    idft(Ibodies2, Jbodies);
  }

  void realPart(real_t (*Ibodies2)[4], real_t (*Jbodies)[4], int (*Leafs)[2]) {
    int nunit = 1 << maxLevel;
    int nmin = -nunit;
    int nmax = 2 * nunit - 1;
#pragma omp parallel for
    for (int i=0; i<numLeafs; i++) {
      int ix[3] = {0, 0, 0};
      getIndex(ix,i);
      int jxmin[3];
      for_3 jxmin[d] = MAX(nmin, ix[d] - 2);
      int jxmax[3];
      for_3 jxmax[d] = MIN(nmax, ix[d] + 2);
      int jx[3];
      for (jx[2]=jxmin[2]; jx[2]<=jxmax[2]; jx[2]++) {
	for (jx[1]=jxmin[1]; jx[1]<=jxmax[1]; jx[1]++) {
	  for (jx[0]=jxmin[0]; jx[0]<=jxmax[0]; jx[0]++) {
	    int jxp[3];
	    for_3 jxp[d] = (jx[d] + nunit) % nunit;
	    int j = getKey(jxp,maxLevel);
	    real_t Xperiodic[3] = {0, 0, 0};
	    for_3 jxp[d] = (jx[d] + nunit) / nunit;
	    for_3 Xperiodic[d] = (jxp[d] - 1) * 2 * R0;
	    P2PEwald(Leafs[i][0],Leafs[i][1],Leafs[j][0],Leafs[j][1],Xperiodic,
		     Ibodies2,Jbodies);
	  }
	}
      }
    }
    for (int i=0; i<numBodies; i++) {
      Ibodies2[i][0] -= M_2_SQRTPI * Jbodies[i][3] * alpha;
    }
  }
};
