#include <cassert>
#include <cstdlib>
#include <cstdio>
#include "kernels2.h"

#ifdef SAKURA
void decomposeSpacePermute(int N, float * Y, float * X, uint * keys, 
			   uint *permutation_vectors, int maxlev);
#endif

class Fmm : public Kernel {
private:
  inline void getIndex(int i, int *ix, real_t diameter) {
#if NOWRAP
    i = (i / 3) * 3;
#endif
    for_3 ix[d] = int((Jbodies[i][d] + R0 - X0[d]) / diameter);
  }

  void sort(real_t (*bodies)[4], real_t (*bodies2)[4], int *Index, int *Index2, int *key) const {
    int Imax = key[0];
    int Imin = key[0];
    for( int i=0; i<numBodies; i++ ) {
      Imax = MAX(Imax,key[i]);
      Imin = MIN(Imin,key[i]);
    }
    int numBucket = Imax - Imin + 1;
    int *bucket = new int [numBucket];
    for (int i=0; i<numBucket; i++) bucket[i] = 0;
    for (int i=0; i<numBodies; i++) bucket[key[i]-Imin]++;
    for (int i=1; i<numBucket; i++) bucket[i] += bucket[i-1];
    for (int i=numBodies-1; i>=0; i--) {
      bucket[key[i]-Imin]--;
      int inew = bucket[key[i]-Imin];
      Index2[inew] = Index[i];
      for_4 bodies2[inew][d] = bodies[i][d];
    }
    delete[] bucket;
  }

public:
  void allocate(int NumBodies, int MaxLevel, int NumNeighbors, int NumImages) {
    numBodies = NumBodies;
    maxLevel = MaxLevel;
    numNeighbors = NumNeighbors;
    numImages = NumImages;
    numCells = ((1 << 3 * (maxLevel + 1)) - 1) / 7;
    numLeafs = 1 << 3 * maxLevel;
    Ibodies = new real_t [numBodies][4]();
    Ibodies2 = new real_t [numBodies][4]();
    Jbodies = new real_t [numBodies][4]();
    Multipole = new real_t [numCells][MTERM]();
    Local = new real_t [numCells][LTERM]();
    Leafs = new int [numLeafs][2]();
    Index = new int [numBodies];
    Index2 = new int [numBodies];
    for (int i=0; i<numCells; i++) {
      for_m Multipole[i][m] = 0;
      for_l Local[i][l] = 0;
    }
    for (int i=0; i<numLeafs; i++) {
      Leafs[i][0] = Leafs[i][1] = 0;
    }
  }

  void deallocate() {
    delete[] Ibodies;
    delete[] Ibodies2;
    delete[] Jbodies;
    delete[] Multipole;
    delete[] Local;
    delete[] Leafs;
    delete[] Index;
    delete[] Index2;
  }

  void initBodies(real_t cycle) {
    int ix[3] = {0, 0, 0};
    R0 = cycle * .5;
    for_3 X0[d] = R0;
    srand48(0);
    real_t average = 0;
    for (int i=0; i<numBodies; i++) {
      Jbodies[i][0] = 2 * R0 * (drand48() + ix[0]);
      Jbodies[i][1] = 2 * R0 * (drand48() + ix[1]);
      Jbodies[i][2] = 2 * R0 * (drand48() + ix[2]);
      Jbodies[i][3] = (drand48() - .5) / numBodies;
      average += Jbodies[i][3];
    }
    average /= numBodies;
    for (int i=0; i<numBodies; i++) {
      Jbodies[i][3] -= average;
    }
  }

  void encode(real_t (*Jbodies)[4], int *key, int N, int offset, real_t diameter, real_t R0, real_t (X0)[3], int maxLevel){
    int ix[3] = {0, 0, 0};
    for(int i=0; i<N; i++){
      for_3 ix[d] = int((Jbodies[i+offset][d] + R0 - X0[d]) / diameter);
      key[i+offset] = getKey(ix,maxLevel);
    }
  }

#ifndef SAKURA
  void sortBodies() {
    int *key = new int [numBodies];
    real_t diameter = 2 * R0 / (1 << maxLevel);
    int ix[3] = {0, 0, 0};
    for (int i=0; i<numBodies; i++) {
      getIndex(i,ix,diameter);
      key[i] = getKey(ix,maxLevel);
    }
    sort(Jbodies, Ibodies, Index, Index2, key);
    for (int i=0; i<numBodies; i++) {
      Index[i] = Index2[i];
      for_4 Jbodies[i][d] = Ibodies[i][d];
      for_4 Ibodies[i][d] = 0;
    }
    delete[] key;
  }
#else
  void sortBodies() {
    uint *key = new uint [numBodies];
    uint *permutation_vector = new uint [numBodies];
    float *X = (float*)Jbodies;
    float *Y = (float*)Ibodies;
    real_t diameter = 2 * R0 / (1 << maxLevel);
    int ix[3] = {0, 0, 0};
    for (int i=0; i<numBodies; i++) {
      getIndex(i,ix,diameter);
      key[i] = getKey(ix,maxLevel);
    }
    decomposeSpacePermute(numBodies, Y, X, key, permutation_vector, maxLevel);
    for (int i=0; i<numBodies; i++) {
      for_4 Jbodies[i][d] = Ibodies[i][d];
      for_4 Ibodies[i][d] = 0;
    }

    delete[] key;
    delete[] permutation_vector;
  }
#endif

  void fillLeafs() {
    real_t diameter = 2 * R0 / (1 << maxLevel);
    int ix[3] = {0, 0, 0};
    getIndex(0,ix,diameter);
    int ileaf = getKey(ix, maxLevel, false);
    Leafs[ileaf][0] = 0;
    for (int i=0; i<numBodies; i++) {
      getIndex(i,ix,diameter);
      int inew = getKey(ix, maxLevel, false);
      if (ileaf != inew) {
        Leafs[ileaf][1] = Leafs[inew][0] = i;
        ileaf = inew;
      }
    }
    Leafs[ileaf][1] = numBodies;
    for (int i=0; i<numLeafs; i++) {
      if (Leafs[i][1] == Leafs[i][0]) printf("Warning: Cell %d is empty.\n",i);
    }
  }

  void dipoleCorrection() {
    real_t dipole[3] = {0, 0, 0};
    for (int i=0; i<numBodies; i++) {
      for_3 dipole[d] += (Jbodies[i][d] - X0[d]) * Jbodies[i][3];
    }
    real_t norm = dipole[0] * dipole[0] + dipole[1] * dipole[1] + dipole[2] * dipole[2];
    real_t cycle = 2 * R0;
    real_t coef = 4 * M_PI / (3 * cycle * cycle * cycle);
    for (int i=0; i<numBodies; i++) {
      Ibodies[i][0] -= coef * norm / numBodies / Jbodies[i][3];
      for_3 Ibodies[i][d+1] -= coef * dipole[d];
    }
  }

  void direct() {
    real_t Ibodies3[4], Jbodies2[4], dX[3];
    int range = (pow(3,numImages) - 1) / 2;
    for (int i=0; i<100; i++) {
      for_4 Ibodies3[d] = 0;
      for_4 Jbodies2[d] = Jbodies[i][d];
      int jx[3];
      for (jx[2]=-range; jx[2]<=range; jx[2]++) {
	for (jx[1]=-range; jx[1]<=range; jx[1]++) {
	  for (jx[0]=-range; jx[0]<=range; jx[0]++) {	
	    for (int j=0; j<numBodies; j++) {
	      for_3 dX[d] = Jbodies2[d] - Jbodies[j][d] - jx[d] * 2 * R0;
	      real_t R2 = dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2];
	      real_t invR2 = R2 == 0 ? 0 : 1.0 / R2;
	      real_t invR = Jbodies[j][3] * sqrtf(invR2);
	      for_3 dX[d] *= invR2 * invR;
	      Ibodies3[0] += invR;
	      Ibodies3[1] -= dX[0];
	      Ibodies3[2] -= dX[1];
	      Ibodies3[3] -= dX[2];
	    }
	  }
	}
      }
      for_4 Ibodies2[i][d] = Ibodies3[d];
    }
  }

  void verify(int numTargets, real_t & potDif, real_t & potNrm, real_t & accDif, real_t & accNrm) {
    real_t potSum = 0, potSum2 = 0;
    for (int i=0; i<numTargets; i++) {
      potSum += Ibodies[i][0] * Jbodies[i][3];
      potSum2 += Ibodies2[i][0] * Jbodies[i][3];
    }
    potDif = (potSum - potSum2) * (potSum - potSum2);
    potNrm = potSum2 * potSum2;
    for (int i=0; i<numTargets; i++) {
      for_3 accDif += (Ibodies[i][d+1] - Ibodies2[i][d+1]) * (Ibodies[i][d+1] - Ibodies2[i][d+1]);
      for_3 accNrm += (Ibodies2[i][d+1] * Ibodies2[i][d+1]);
    }
  }
};
