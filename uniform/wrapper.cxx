#include "parallelfmm.h"

static const double Celec = 332.0716;

ParallelFMM * FMM;

int wrap(double x[], double cycle) {
  int iwrap = 0;
  for (int d=0; d<3; d++) {
    if (x[d] < -cycle / 2) {
      x[d] += cycle;
      iwrap |= 1 << d;
    }
    if (x[d] > cycle / 2) {
      x[d] -= cycle;
      iwrap |= 1 << d;
    }
  }
  return iwrap;
}

void unwrap(double x[], double cycle, int iwrap) {
  for (int d=0; d<3; d++) {
    if ((iwrap >> d) & 1) x[d] += (x[d] > 0 ? -cycle : cycle);
  }
}

extern "C" void fmm_init_(int & N, int & images, int & verbose) {
  const int ncrit = 100;
  const int maxLevel = N >= ncrit ? 1 + int(log(N / ncrit)/M_LN2/3) : 0;
  const int numImages = images;
  FMM = new ParallelFMM;
  FMM->allocate(N, maxLevel, numImages);
  logger::verbose = true;
  FMM->printNow = false;
}

extern "C" void fmm_finalize_() {
  FMM->deallocate();
  delete FMM;
}

extern "C" void fmm_partition_(int & nglobal, int * icpumap, double * x, double * q,
			       double * xold, double & cycle) {
  const int gatherLevel = 1;
  logger::startTimer("Partition");
  FMM->partitioner(gatherLevel);
  logger::stopTimer("Partition");

  logger::startTimer("Copy in");
  const int shift = 29;
  const int mask = ~(0x7U << shift);
  int nlocal = 0;
  for (int i=0; i<nglobal; i++) {
    if (icpumap[i] == 1) nlocal++;
  }
  assert(FMM->numBodies == nlocal);
  int ix[3] = {0, 0, 0};
  FMM->R0 = 0.5 * cycle / FMM->numPartition[FMM->maxGlobLevel][0];
  for_3d FMM->RGlob[d] = FMM->R0 * FMM->numPartition[FMM->maxGlobLevel][d];
  FMM->getGlobIndex(ix,FMM->MPIRANK,FMM->maxGlobLevel);
  for_3d FMM->X0[d] = 2 * FMM->R0 * (ix[d] + .5);
  for (int i=0,b=0; i<nglobal; i++) {
    if (icpumap[i] == 1) {
      for_3d FMM->Jbodies[b][d] = x[3*i+d];
      FMM->Jbodies[b][3] = q[i];
      for_3d FMM->Ibodies[b][d] = xold[3*i+d];
      int iwrap = wrap(FMM->Jbodies[b], cycle);
      FMM->Index[b] = i | (iwrap << shift);
      for_3d FMM->Jbodies[b][d] += FMM->RGlob[d];
      b++;
    }
  }
  logger::stopTimer("Copy in");
 
  logger::startTimer("Partition Comm.");
  FMM->partitionComm();
  logger::stopTimer("Partition Comm.");

  logger::startTimer("Copy out"); 
  for (int i=0; i<nglobal; i++) {
    icpumap[i] = 0;
  }
  for (int b=0; b<nlocal; b++) {
    for_3d FMM->Jbodies[b][d] -= FMM->RGlob[d];
    int iwrap = unsigned(FMM->Index[b]) >> shift;
    unwrap(FMM->Jbodies[b], cycle, iwrap);
    int i = FMM->Index[b] & mask;
    for_3d x[3*i+d] = FMM->Jbodies[b][d];
    q[i] = FMM->Jbodies[b][3];
    for_3d xold[3*i+d] = FMM->Ibodies[b][d];
    icpumap[i] = 1;
  }
  logger::stopTimer("Copy out"); 
}

extern "C" void fmm_coulomb_(int & nglobal, int * icpumap,
			     double * x, double * q, double * p, double * f,
			     double & cycle) {
  const int shift = 29;
  const int mask = ~(0x7U << shift);
  int nlocal = 0;
  for (int i=0; i<nglobal; i++) {
    if (icpumap[i] == 1) {
      nlocal++;
    } else {
      icpumap[i] = 0;
    }
  }
  logger::startTimer("Total FMM");
  for (int i=0,b=0; i<nglobal; i++) {
    if (icpumap[i] == 1) {
      for_3d FMM->Jbodies[b][d] = x[3*i+d];
      FMM->Jbodies[b][3] = q[i];
      for_4d FMM->Ibodies[b][d] = 0;
      int iwrap = wrap(FMM->Jbodies[b], cycle);
      FMM->Index[b] = i | (iwrap << shift);
      for_3d FMM->Jbodies[b][d] += FMM->RGlob[d];
      b++;
    }
  }

  logger::startTimer("Grow tree");
  FMM->sortBodies();
  FMM->buildTree();
  logger::stopTimer("Grow tree");

  logger::startTimer("Upward pass");
  FMM->upwardPass();
  logger::stopTimer("Upward pass");

  FMM->periodicM2L();

  logger::startTimer("Downward pass");
  FMM->globL2L();
  logger::stopTimer("Downward pass");

  FMM->downwardPass();

  logger::stopTimer("Total FMM");
  logger::printTitle("Total runtime");
  logger::printTime("Total FMM");

  for (int b=0; b<nlocal; b++) {
    int i = FMM->Index[b] & mask;
    p[i] += FMM->Ibodies[b][0] * FMM->Jbodies[b][3] * Celec;
    for_3d f[3*i+d] += FMM->Ibodies[b][d+1] * FMM->Jbodies[b][3] * Celec;
    for_3d FMM->Jbodies[b][d] -= FMM->RGlob[d];
    int iwrap = unsigned(FMM->Index[b]) >> shift;
    unwrap(FMM->Jbodies[b], cycle, iwrap);
    for_3d x[3*i+d] = FMM->Jbodies[b][d];
    q[i] = FMM->Jbodies[b][3];
  }
}

extern "C" void ewald_coulomb_(int & nglobal, int * icpumap, double * x, double * q, double * p, double * f,
			       int & ksize, double & alpha, double & sigma, double & cutoff, double & cycle) {
}

extern "C" void direct_coulomb_(int & nglobal, int * icpumap, double * x, double * q, double * p, double * f, double & cycle) {
  logger::startTimer("Direct Coulomb");
  int images = FMM->numImages;
  int prange = 0;
  for (int i=0; i<images; i++) {
    prange += int(std::pow(3.,i));
  }
  real Xperiodic[3];
  for (int i=0; i<nglobal; i++) {
    if (icpumap[i] == 1) {
      real pp = 0, fx = 0, fy = 0, fz = 0;
      for (int ix=-prange; ix<=prange; ix++) {
	for (int iy=-prange; iy<=prange; iy++) {
	  for (int iz=-prange; iz<=prange; iz++) {
	    Xperiodic[0] = ix * cycle;
	    Xperiodic[1] = iy * cycle;
	    Xperiodic[2] = iz * cycle;
	    for (int j=0; j<nglobal; j++) {
	      real dX[3];
	      for (int d=0; d<3; d++) dX[d] = x[3*i+d] - x[3*j+d] - Xperiodic[d];
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
      }
      p[i] += pp * q[i] * Celec;
      f[3*i+0] -= fx * q[i] * Celec;
      f[3*i+1] -= fy * q[i] * Celec;
      f[3*i+2] -= fz * q[i] * Celec;
    }
  }
  real dipole[3] = {0, 0, 0};
  for (int i=0; i<nglobal; i++) {
    for (int d=0; d<3; d++) dipole[d] += x[3*i+d] * q[i];
  }
  real norm = 0;
  for (int d=0; d<3; d++) {
    norm += dipole[d] * dipole[d];
  }
  real coef = 4 * M_PI / (3 * cycle * cycle * cycle) * Celec;
  for (int i=0; i<nglobal; i++) {
    if (icpumap[i] == 1) {
      p[i] -= coef * norm / nglobal;
      f[3*i+0] -= coef * dipole[0] * q[i];
      f[3*i+1] -= coef * dipole[1] * q[i];
      f[3*i+2] -= coef * dipole[2] * q[i];
    }
  }
  logger::stopTimer("Direct Coulomb");
}

extern "C" void coulomb_exclusion_(int & nglobal, int * icpumap,
				   double * x, double * q, double * p, double * f,
				   double & cycle, int * numex, int * natex) {
  logger::startTimer("Coulomb Exclusion");
  for (int i=0, ic=0; i<nglobal; i++) {
    if (icpumap[i] == 1) {
      real pp = 0, fx = 0, fy = 0, fz = 0;
      for (int jc=0; jc<numex[i]; jc++, ic++) {
	int j = natex[ic]-1;
	real dX[3];
	for (int d=0; d<3; d++) dX[d] = x[3*i+d] - x[3*j+d];
        wrap(dX, cycle);
	real R2 = dX[0] * dX[0] + dX[1] * dX[1] + dX[2] + dX[2];
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
    } else {
      ic += numex[i];
    }
  }
  logger::stopTimer("Coulomb Exclusion");
}

extern "C" void fmm_vanderwaals_(int & nglobal, int * icpumap, int * atype,
				 double * x, double * p, double * f,
				 double & cuton, double & cutoff, double & cycle,
				 int & numTypes, double * rscale, double * gscale, double * fgscale) {
}

extern "C" void direct_vanderwaals_(int & nglobal, int * icpumap, int * atype,
				    double * x, double * p, double * f,
				    double & cuton, double & cutoff, double & cycle,
				    int & numTypes, double * rscale, double * gscale, double * fgscale) {
  logger::startTimer("Direct VdW");
  for (int i=0; i<nglobal; i++) {
    if (icpumap[i] == 1) {
      int atypei = atype[i]-1;
      real pp = 0, fx = 0, fy = 0, fz = 0;
      for (int j=0; j<nglobal; j++) {
	real dX[3];
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
  }
  logger::stopTimer("Direct VdW");
}

extern "C" void vanderwaals_exclusion_(int & nglobal, int * icpumap, int * atype,
				       double * x, double * p, double * f,
				       double & cuton, double & cutoff, double & cycle,
				       int & numTypes, double * rscale, double * gscale,
				       double * fgscale, int * numex, int * natex) {
  logger::startTimer("VdW Exclusion");
  for (int i=0, ic=0; i<nglobal; i++) {
    if (icpumap[i] == 1) {
      int atypei = atype[i]-1;
      for (int jc=0; jc<numex[i]; jc++, ic++) {
	int j = natex[ic]-1;
	real dX[3];
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
    } else {
      ic += numex[i];
    }
  }
  logger::stopTimer("VdW Exclusion");
}
