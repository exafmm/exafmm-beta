#include <cassert>
#include <cstdlib>
#include <cstdio>
#include "kernels.h"

class Fmm : public Kernel {
private:
  void sort(real (*bodies)[4], real (*bodies2)[4], int *key) const {
    int Imax = key[0];
    int Imin = key[0];
    for( int i=0; i<numBodies; i++ ) {
      Imax = MAX(Imax,key[i]);
      Imin = MIN(Imin,key[i]);
    }
    int numBucket = Imax - Imin + 1;
    int *bucket = new int [numBucket];
    for( int i=0; i<numBucket; i++ ) bucket[i] = 0;
    for( int i=0; i<numBodies; i++ ) bucket[key[i]-Imin]++;
    for( int i=1; i<numBucket; i++ ) bucket[i] += bucket[i-1];
    for( int i=numBodies-1; i>=0; --i ) {
      bucket[key[i]-Imin]--;
      int inew = bucket[key[i]-Imin];
      for_4d bodies2[inew][d] = bodies[i][d];
    }
    delete[] bucket;
  }

public:
  void allocate(int NumBodies, int MaxLevel, int NumNeighbors) {
    numBodies = NumBodies;
    maxLevel = MaxLevel;
    numNeighbors = NumNeighbors;
    numCells = ((1 << 3 * (maxLevel + 1)) - 1) / 7;
    numLeafs = 1 << 3 * maxLevel;
    Ibodies = new real [numBodies][4]();
    Jbodies = new real [numBodies][4]();
    Multipole = new real [numCells][MTERM]();
    Local = new real [numCells][LTERM]();
    Leafs = new int [numLeafs][2]();
    for( int i=0; i<numCells; i++ ) {
      for_m Multipole[i][m] = 0;
      for_l Local[i][l] = 0;
    }
    for( int i=0; i<numLeafs; i++ ) {
      Leafs[i][0] = Leafs[i][1] = 0;
    }
  }

  void deallocate() {
    delete[] Ibodies;
    delete[] Jbodies;
    delete[] Multipole;
    delete[] Local;
    delete[] Leafs;
  }

  void initBodies() {
    int ix[3] = {0, 0, 0};
    R0 = M_PI;
    for_3d X0[d] = R0;
    srand48(0);
    real average = 0;
    for( int i=0; i<numBodies; i++ ) {
      Jbodies[i][0] = 2 * R0 * (drand48() + ix[0]);
      Jbodies[i][1] = 2 * R0 * (drand48() + ix[1]);
      Jbodies[i][2] = 2 * R0 * (drand48() + ix[2]);
      Jbodies[i][3] = (drand48() - .5) / numBodies;
      average += Jbodies[i][3];
    }
    average /= numBodies;
    for( int i=0; i<numBodies; i++ ) {
      Jbodies[i][3] -= average;
    }
  }

  void sortBodies() const {
    int *key = new int [numBodies];
    real diameter = 2 * R0 / (1 << maxLevel);
    int ix[3] = {0, 0, 0};
    for( int i=0; i<numBodies; i++ ) {
      for_3d ix[d] = int((Jbodies[i][d] + R0 - X0[d]) / diameter);
      key[i] = getKey(ix,maxLevel);
    }
    sort(Jbodies,Ibodies,key);
    for( int i=0; i<numBodies; i++ ) {
      for_4d Jbodies[i][d] = Ibodies[i][d];
      for_4d Ibodies[i][d] = 0;
    }
    delete[] key;
  }

  void fillLeafs() const {
    real diameter = 2 * R0 / (1 << maxLevel);
    int ix[3] = {0, 0, 0};
    for_3d ix[d] = int((Jbodies[0][d] + R0 - X0[d]) / diameter);
    int ileaf = getKey(ix,maxLevel,false);
    Leafs[ileaf][0] = 0;
    for( int i=0; i<numBodies; i++ ) {
      for_3d ix[d] = int((Jbodies[i][d] + R0 - X0[d]) / diameter);
      int inew = getKey(ix,maxLevel,false);
      if( ileaf != inew ) {
        Leafs[ileaf][1] = Leafs[inew][0] = i;
        ileaf = inew;
      }
    }
    Leafs[ileaf][1] = numBodies;
    for( int i=0; i<numLeafs; i++ ) {
      if( Leafs[i][1] == Leafs[i][0] ) printf("Warning: Cell %d is empty.\n",i);
    }
  }

  void verify(real & potDif, real & potNrm, real & accDif, real & accNrm) {
    for (int i=0; i<100; i++) {
      real Ibodies2[4] = {0, 0, 0, 0};
      real Jbodies2[4], dX[3];
      for_4d Jbodies2[d] = Jbodies[i][d];
      for (int j=0; j<numBodies; j++) {
	for_3d dX[d] = Jbodies2[d] - Jbodies[j][d];
	real R2 = dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2];
	real invR2 = R2 == 0 ? 0 : 1.0 / R2;
	real invR = Jbodies[j][3] * sqrtf(invR2);
	for_3d dX[d] *= invR2 * invR;
	Ibodies2[0] += invR;
	Ibodies2[1] -= dX[0];
	Ibodies2[2] -= dX[1];
	Ibodies2[3] -= dX[2];
      }
      potDif += (Ibodies[i][0] - Ibodies2[0]) * (Ibodies[i][0] - Ibodies2[0]);
      potNrm += Ibodies2[0] * Ibodies2[0];
      for_3d accDif += (Ibodies[i][d+1] - Ibodies2[d+1]) * (Ibodies[i][d+1] - Ibodies2[d+1]);
      for_3d accNrm += (Ibodies2[d+1] * Ibodies2[d+1]);
    }
  }
};
