#include "ewald2.h"
#include "fmm.h"
#include "logger2.h"
#include "vdw.h"

#ifdef SAKURA
#define DIM 3
#define NN 26
#define FN 37
#define CN 152
int common_stencil[DIM*CN] = {0};
int far_stencil[8*DIM*FN] = {0};
int near_stencil[8*DIM*NN] = {0};
void formInteractionStencil(int *common_stencil, int *far_stencil, int *near_stencil);
#endif

static const double Celec = 332.0716;
int numImages;
Fmm * FMM;

int wrap(real bodies[4], real cycle) {
  int iwrap = 0;
  for (int d=0; d<3; d++) {
    if(bodies[d] < -cycle / 2) {
      bodies[d] += cycle;
      iwrap |= 1 << d;
    }
    if(bodies[d] > cycle / 2) {
      bodies[d] -= cycle;
      iwrap |= 1 << d;
    }
  }
  return iwrap;
}

void unwrap(real bodies[4], real cycle, int iwrap) {
  for (int d=0; d<3; d++) {
    if((iwrap >> d) & 1) bodies[d] += (bodies[d] > 0 ? -cycle : cycle);
  }
}

extern "C" void fmm_init_(int & images, double &, int & verbose) {
  numImages = images;
  FMM = new Fmm;
  logger::verbose = verbose;
#ifdef SAKURA
  formInteractionStencil(common_stencil, far_stencil, near_stencil);
#endif
}

extern "C" void fmm_finalize_() {
  delete FMM;
}

extern "C" void fmm_partition_(int &, int *, double *, double *, double *, double &) {}

extern "C" void fmm_coulomb_(int & numBodies, int *,
			     double * x, double * q, double * p, double * f,
			     double & cycle) {
  logger::startTimer("Total FMM");
  const int shift = 29;
  const int mask = ~(0x7U << shift);
  const int ncrit = 100;
  const int maxLevel = numBodies >= ncrit ? 1 + int(log(numBodies / ncrit)/M_LN2/3) : 0;
  const int numNeighbors = 1;
  FMM->allocate(numBodies, maxLevel, numNeighbors, numImages);
  FMM->R0 = cycle * .5;
  for_3 FMM->X0[d] = FMM->R0;
  int iwrap = 0;
  for (int i=0; i<numBodies; i++) {
    for_3 FMM->Jbodies[i][d] = x[3*i+d];
    FMM->Jbodies[i][3] = q[i];
#if NOWRAP
    if(i % 3 == 0) iwrap = wrap(FMM->Jbodies[i], cycle);
    else unwrap(FMM->Jbodies[i], cycle, iwrap);
#else
    iwrap = wrap(FMM->Jbodies[i], cycle);
#endif
    FMM->Index[i] = i | (iwrap << shift);
    for_3 FMM->Jbodies[i][d] += FMM->R0;
    for_4 FMM->Ibodies[i][d] = 0;
  }
  FMM->sortBodies();
  FMM->fillLeafs();
  FMM->P2M();
  FMM->M2M();
#ifndef SAKURA
  FMM->M2L();
#else
  FMM->M2L(common_stencil, far_stencil);
#endif
  FMM->L2L();
  FMM->L2P();
#ifndef SAKURA
  FMM->P2P();
#else
  FMM->P2P(near_stencil);
#endif
  FMM->dipoleCorrection();
  for (int i=0; i<numBodies; i++) {
    int ii = FMM->Index[i] & mask;
    int iwrap = unsigned(FMM->Index[i]) >> shift;
    for_3 FMM->Jbodies[i][d] -= FMM->R0;
    unwrap(FMM->Jbodies[i], cycle, iwrap);
    p[ii] += FMM->Ibodies[i][0] * FMM->Jbodies[i][3] * Celec;
    for_3 f[3*ii+d] += FMM->Ibodies[i][d+1] * FMM->Jbodies[i][3] * Celec;
  }
  FMM->deallocate();
  logger::stopTimer("Total FMM");
}

extern "C" void ewald_coulomb_(int & numBodies, int *, double * x, double * q, double * p, double * f,
			       int & ksize, double & alpha, double & sigma, double & cutoff, double & cycle) {
  logger::startTimer("Total Ewald");
  const int shift = 29;
  const int mask = ~(0x7U << shift);
  const int ncrit = 100;
  const int maxLevel = numBodies >= ncrit ? 1 + int(log(numBodies / ncrit)/M_LN2/3) : 0;
  const int numNeighbors = 1;
  Ewald ewald(numBodies, maxLevel, cycle, ksize, alpha, sigma, cutoff);
  FMM->allocate(numBodies, maxLevel, numNeighbors, numImages);
  FMM->R0 = cycle * .5;
  for_3 FMM->X0[d] = FMM->R0;
  int iwrap = 0;
  for (int i=0; i<numBodies; i++) {
    for_3 FMM->Jbodies[i][d] = x[3*i+d];
    FMM->Jbodies[i][3] = q[i];
#if NOWRAP
    if(i % 3 == 0) iwrap = wrap(FMM->Jbodies[i], cycle);
    else unwrap(FMM->Jbodies[i], cycle, iwrap);
#else
    iwrap = wrap(FMM->Jbodies[i], cycle);
#endif
    FMM->Index[i] = i | (iwrap << shift);
    for_3 FMM->Jbodies[i][d] += FMM->R0;
    for_4 FMM->Ibodies[i][d] = 0;
  }
  FMM->sortBodies();
  FMM->fillLeafs();
  ewald.wavePart(FMM->Ibodies, FMM->Jbodies);
  ewald.realPart(FMM->Ibodies, FMM->Jbodies, FMM->Leafs);
  for (int i=0; i<numBodies; i++) {
    int ii = FMM->Index[i] & mask;
    int iwrap = unsigned(FMM->Index[i]) >> shift;
    for_3 FMM->Jbodies[i][d] -= FMM->R0;
    unwrap(FMM->Jbodies[i], cycle, iwrap);
    p[ii] += FMM->Ibodies[i][0] * FMM->Jbodies[i][3] * Celec;
    for_3 f[3*ii+d] += FMM->Ibodies[i][d+1] * FMM->Jbodies[i][3] * Celec;
  }
  FMM->deallocate();
  logger::stopTimer("Total Ewald");
}

extern "C" void direct_coulomb_(int & numBodies, int *, double * x, double * q, double * p, double * f, double & cycle) {
  logger::startTimer("Direct Coulomb");
  int prange = 0;
  for (int i=0; i<numImages; i++) {
    prange += int(std::pow(3.,i));
  }
  real Xperiodic[3];
  for (int i=0; i<numBodies; i++) {
    real pp = 0, fx = 0, fy = 0, fz = 0;
    for (int ix=-prange; ix<=prange; ix++) {
      for (int iy=-prange; iy<=prange; iy++) {
	for (int iz=-prange; iz<=prange; iz++) {
	  Xperiodic[0] = ix * cycle;
	  Xperiodic[1] = iy * cycle;
	  Xperiodic[2] = iz * cycle;
	  for (int j=0; j<numBodies; j++) {
	    real dX[4];
	    for_3 dX[d] = x[3*i+d] - x[3*j+d] - Xperiodic[d];
	    real R2 = dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2];
	    real invR = 1 / std::sqrt(R2);
	    if (R2 == 0) invR = 0;
	    real invR3 = q[j] * invR * invR * invR;
	    pp += q[j] * invR;
	    fx += dX[0] * invR3;
	    fy += dX[1] * invR3;
	    fz += dX[2] * invR3;
	  }
	}
      }
      p[i] += pp * q[i] * Celec;
      f[3*i+0] -= fx * q[i] * Celec;
      f[3*i+1] -= fy * q[i] * Celec;
      f[3*i+2] -= fz * q[i] * Celec;
    }
  }
  real dipole[3] = {0, 0, 0};
  for (int i=0; i<numBodies; i++) {
    for_3 dipole[d] += x[3*i+d] * q[i];
  }
  real norm = 0;
  for_3 norm += dipole[d] * dipole[d];
  real coef = 4 * M_PI / (3 * cycle * cycle * cycle) * Celec;
  for (int i=0; i<numBodies; i++) {
    p[i] -= coef * norm / numBodies;
    f[3*i+0] -= coef * dipole[0] * q[i];
    f[3*i+1] -= coef * dipole[1] * q[i];
    f[3*i+2] -= coef * dipole[2] * q[i];
  }
  logger::stopTimer("Direct Coulomb");
}

extern "C" void coulomb_exclusion_(int & numBodies, int *,
				   double * x, double * q, double * p, double * f,
				   double & cycle, int * numex, int * natex) {
  logger::startTimer("Coulomb Exclusion");
  for (int i=0, ic=0; i<numBodies; i++) {
    real pp = 0, fx = 0, fy = 0, fz = 0;
    for (int jc=0; jc<numex[i]; jc++, ic++) {
      int j = natex[ic]-1;
      real dX[4];
      for (int d=0; d<3; d++) dX[d] = x[3*i+d] - x[3*j+d];
      wrap(dX, cycle);
      real R2 = dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2];
      real invR = 1 / std::sqrt(R2);
      if (R2 == 0) invR = 0;
      real invR3 = q[j] * invR * invR * invR;
      pp += q[j] * invR;
      fx += dX[0] * invR3;
      fy += dX[1] * invR3;
      fz += dX[2] * invR3;
    }
    p[i] -= pp * q[i] * Celec;
    f[3*i+0] += fx * q[i] * Celec; 
    f[3*i+1] += fy * q[i] * Celec;
    f[3*i+2] += fz * q[i] * Celec;
  }
  logger::stopTimer("Coulomb Exclusion");
}

extern "C" void fmm_vanderwaals_(int & numBodies, int *, int * atype,
				 double * x, double * p, double * f,
				 double & cuton, double & cutoff, double & cycle,
				 int & numTypes, double * rscale, double * gscale, double * fgscale) {
  logger::startTimer("Total VdW");
  const int shift = 29;
  const int mask = ~(0x7U << shift);
  const int ncrit = 100;
  const int maxLevel = numBodies >= ncrit ? 1 + int(log(numBodies / ncrit)/M_LN2/3) : 0;
  const int numNeighbors = 1;
  VdW VDW(numBodies, maxLevel, numTypes, cycle, cuton, cutoff, rscale, gscale, fgscale);
  FMM->allocate(numBodies, maxLevel, numNeighbors, numImages);
  FMM->R0 = cycle * .5;
  for_3 FMM->X0[d] = FMM->R0;
  int iwrap = 0;
  for (int i=0; i<numBodies; i++) {
    for_3 FMM->Jbodies[i][d] = x[3*i+d];
    FMM->Jbodies[i][3] = atype[i] - .5;
#if NOWRAP
    if(i % 3 == 0) iwrap = wrap(FMM->Jbodies[i], cycle);
    else unwrap(FMM->Jbodies[i], cycle, iwrap);
#else
    iwrap = wrap(FMM->Jbodies[i], cycle);
#endif
    FMM->Index[i] = i | (iwrap << shift);
    for_3 FMM->Jbodies[i][d] += FMM->R0;
    for_4 FMM->Ibodies[i][d] = 0;
  }
  FMM->sortBodies();
  FMM->fillLeafs();
  VDW.evaluate(FMM->Ibodies, FMM->Jbodies, FMM->Leafs);
  for (int i=0; i<numBodies; i++) {
    int ii = FMM->Index[i] & mask;
    int iwrap = unsigned(FMM->Index[i]) >> shift;
    for_3 FMM->Jbodies[i][d] -= FMM->R0;
    unwrap(FMM->Jbodies[i], cycle, iwrap);
    p[ii] += FMM->Ibodies[i][0];
    for_3 f[3*ii+d] += FMM->Ibodies[i][d+1];
  }
  FMM->deallocate();
  logger::stopTimer("Total VdW");
}

extern "C" void direct_vanderwaals_(int & numBodies, int *, int * atype,
				    double * x, double * p, double * f,
				    double & cuton, double & cutoff, double & cycle,
				    int & numTypes, double * rscale, double * gscale, double * fgscale) {
  logger::startTimer("Direct VdW");
  for (int i=0; i<numBodies; i++) {
    int atypei = atype[i]-1;
    real pp = 0, fx = 0, fy = 0, fz = 0;
    for (int j=0; j<numBodies; j++) {
      real dX[4];
      for (int d=0; d<3; d++) dX[d] = x[3*i+d] - x[3*j+d];
      wrap(dX, cycle);
      real R2 = dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2];
      if (R2 != 0) {
	int atypej = atype[j]-1;
	real rs = rscale[atypei*numTypes+atypej];
	real gs = gscale[atypei*numTypes+atypej];
	real fgs = fgscale[atypei*numTypes+atypej];
	real R2s = R2 * rs;
	real invR2 = 1.0 / R2s;
	real invR6 = invR2 * invR2 * invR2;
	real cuton2 = cuton * cuton;
	real cutoff2 = cutoff * cutoff;
	if (R2 < cutoff2) {
	  real tmp = 0, dtmp = 0;
	  if (cuton2 < R2) {
	    real tmp1 = (cutoff2 - R2) / ((cutoff2-cuton2)*(cutoff2-cuton2)*(cutoff2-cuton2));
	    real tmp2 = tmp1 * (cutoff2 - R2) * (cutoff2 - 3 * cuton2 + 2 * R2);
	    tmp = invR6 * (invR6 - 1) * tmp2;
	    dtmp = invR6 * (invR6 - 1) * 12 * (cuton2 - R2) * tmp1
	      - 6 * invR6 * (invR6 + (invR6 - 1) * tmp2) * tmp2 / R2;
	  } else {
	    tmp = invR6 * (invR6 - 1);
	    dtmp = invR2 * invR6 * (2 * invR6 - 1);
	  }
	  dtmp *= fgs;
	  pp += gs * tmp;
	  fx += dX[0] * dtmp;
	  fy += dX[1] * dtmp;
	  fz += dX[2] * dtmp;
	}
      }
    }
    p[i] += pp;
    f[3*i+0] -= fx;
    f[3*i+1] -= fy;
    f[3*i+2] -= fz;
  }
  logger::stopTimer("Direct VdW");
}

extern "C" void vanderwaals_exclusion_(int & numBodies, int *, int * atype,
				       double * x, double * p, double * f,
				       double & cuton, double & cutoff, double & cycle,
				       int & numTypes, double * rscale, double * gscale,
				       double * fgscale, int * numex, int * natex) {
  logger::startTimer("VdW Exclusion");
  for (int i=0, ic=0; i<numBodies; i++) {
    int atypei = atype[i]-1;
    for (int jc=0; jc<numex[i]; jc++, ic++) {
      int j = natex[ic]-1;
      real dX[4];
      for (int d=0; d<3; d++) dX[d] = x[3*i+d] - x[3*j+d];
      wrap(dX, cycle);
      real R2 = dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2];
      if (R2 != 0) {
	int atypej = atype[j]-1;
	real rs = rscale[atypei*numTypes+atypej];
	real gs = gscale[atypei*numTypes+atypej];
	real fgs = fgscale[atypei*numTypes+atypej];
	real R2s = R2 * rs;
	real invR2 = 1.0 / R2s;
	real invR6 = invR2 * invR2 * invR2;
	real cuton2 = cuton * cuton;
	real cutoff2 = cutoff * cutoff;
	if (R2 < cutoff2) {
	  real tmp = 0, dtmp = 0;
	  if (cuton2 < R2) {
	    real tmp1 = (cutoff2 - R2) / ((cutoff2-cuton2)*(cutoff2-cuton2)*(cutoff2-cuton2));
	    real tmp2 = tmp1 * (cutoff2 - R2) * (cutoff2 - 3 * cuton2 + 2 * R2);
	    tmp = invR6 * (invR6 - 1) * tmp2;
	    dtmp = invR6 * (invR6 - 1) * 12 * (cuton2 - R2) * tmp1
	      - 6 * invR6 * (invR6 + (invR6 - 1) * tmp2) * tmp2 / R2;
	  } else {
	    tmp = invR6 * (invR6 - 1);
	    dtmp = invR2 * invR6 * (2 * invR6 - 1);
	  }
	  dtmp *= fgs;
	  p[i] -= gs * tmp;
	  f[3*i+0] += dX[0] * dtmp;
	  f[3*i+1] += dX[1] * dtmp;
	  f[3*i+2] += dX[2] * dtmp;
	}
      }
    }
  }
  logger::stopTimer("VdW Exclusion");
}
