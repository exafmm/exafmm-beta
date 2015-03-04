#include <cassert>
#include <cmath>
#include "types2.h"

class VdW {
private:
  int numBodies;
  int maxLevel;
  int numLeafs;
  int numTypes;
  real_t R0;
  real_t X0[3];
  real_t cuton;
  real_t cutoff;
  real_t *rscale;
  real_t *gscale;
  real_t *fgscale;

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

  void P2PVdW(int ibegin, int iend, int jbegin, int jend, real_t *Xperiodic,
	      real_t (*Ibodies)[4], real_t (*Jbodies)[4]) const {
    for (int i=ibegin; i<iend; i++) {
      int atypei = int(Jbodies[i][3]);
      real_t Po = 0, Fx = 0, Fy = 0, Fz = 0;
      for (int j=jbegin; j<jend; j++) {
	real_t dX[3];
	for_3 dX[d] = Jbodies[i][d] - Jbodies[j][d] - Xperiodic[d];
	real_t R2 = dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2];
	if (R2 != 0) {
          int atypej = int(Jbodies[j][3]);
          real_t rs = rscale[atypei*numTypes+atypej];
          real_t gs = gscale[atypei*numTypes+atypej];
          real_t fgs = fgscale[atypei*numTypes+atypej];
          real_t R2s = R2 * rs;
          real_t invR2 = 1.0 / R2s;
          real_t invR6 = invR2 * invR2 * invR2;
          real_t cuton2 = cuton * cuton;
          real_t cutoff2 = cutoff * cutoff;
          if (R2 < cutoff2) {
            real_t tmp = 0, dtmp = 0;
            if (cuton2 < R2) {
              real_t tmp1 = (cutoff2 - R2) / ((cutoff2-cuton2)*(cutoff2-cuton2)*(cutoff2-cuton2));
              real_t tmp2 = tmp1 * (cutoff2 - R2) * (cutoff2 - 3 * cuton2 + 2 * R2);
              tmp = invR6 * (invR6 - 1) * tmp2;
              dtmp = invR6 * (invR6 - 1) * 12 * (cuton2 - R2) * tmp1
                - 6 * invR6 * (invR6 + (invR6 - 1) * tmp2) * tmp2 / R2;
            } else {
              tmp = invR6 * (invR6 - 1);
              dtmp = invR2 * invR6 * (2 * invR6 - 1);
            }
            dtmp *= fgs;
            Po += gs * tmp;
            Fx += dX[0] * dtmp;
            Fy += dX[1] * dtmp;
            Fz += dX[2] * dtmp;
          }
        }
      }
      Ibodies[i][0] += Po;
      Ibodies[i][1] -= Fx;
      Ibodies[i][2] -= Fy;
      Ibodies[i][3] -= Fz;
    }
  }
  
public:
  VdW(int _numBodies, int _maxLevel, int _numTypes,
      real_t cycle, real_t _cuton, real_t _cutoff,
      real_t *_rscale, real_t *_gscale, real_t *_fgscale) :
    numTypes(_numTypes), cuton(_cuton), cutoff(_cutoff) {
    numBodies = _numBodies;
    maxLevel = _maxLevel;
    numLeafs = 1 << 3 * maxLevel;
    R0 = cycle * .5;
    for_3 X0[d] = R0;
    rscale = new real_t [numTypes*numTypes];
    gscale = new real_t [numTypes*numTypes];
    fgscale = new real_t [numTypes*numTypes];
    for (int i=0; i<numTypes*numTypes; i++) {
      rscale[i] = _rscale[i];
      gscale[i] = _gscale[i];
      fgscale[i] = _fgscale[i];
    }
  }

  ~VdW() {
    delete[] rscale;
    delete[] gscale;
    delete[] fgscale;
  }
  
  void evaluate(real_t (*Ibodies)[4], real_t (*Jbodies)[4], int (*Leafs)[2]) {
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
	    P2PVdW(Leafs[i][0],Leafs[i][1],Leafs[j][0],Leafs[j][1],Xperiodic,
		   Ibodies,Jbodies);
	  }
	}
      }
    }
  }
};
