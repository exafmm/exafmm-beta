#include <cmath>

typedef double real;
const int PP = 6;
const int MTERM = PP*(PP+1)*(PP+2)/6;
const int LTERM = (PP+1)*(PP+2)*(PP+3)/6;

#include "core.h"

#define for_3d for( int d=0; d<3; d++ )
#define for_4d for( int d=0; d<4; d++ )
#define for_m for( int m=0; m<MTERM; m++ )
#define for_l for( int l=0; l<LTERM; l++ )
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

#define DIM 3

#define NN 26
#define FN 37
#define CN 152

int executeStencilNN(int *nn_codes, int *node_code, int *stencil, int periodic, int lv);
int executeStencilFN(int *fn_codes, int *node_code, int *fn_stencil, int *cn_stencil, int periodic, int lv);

class Kernel {
public:
  int maxLevel;
  int numBodies;
  int numCells;
  int numLeafs;
  int numNeighbors;

  real X0[3];
  real R0;
  real (*Ibodies)[4];
  real (*Jbodies)[4];
  real (*Multipole)[MTERM];
  real (*Local)[LTERM];
  int (*Leafs)[2];

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
  inline int getKey(int *ix, int level, bool levelOffset=true) const {
    int id = 0;
    if( levelOffset ) id = ((1 << 3 * level) - 1) / 7;
    for( int lev=0; lev<level; ++lev ) {
      for_3d id += (ix[d] >> lev) % 2 << (3 * lev + d);
    }
    return id;
  }

public:
  void P2PSum(int ibegin, int iend, int jbegin, int jend) const {
    for( int i=ibegin; i<iend; i++ ) {
      real Po = 0, Fx = 0, Fy = 0, Fz = 0;
      for( int j=jbegin; j<jend; j++ ) {
	real dist[3];
	for_3d dist[d] = Jbodies[i][d] - Jbodies[j][d];
	real R2 = dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];
	real invR2 = R2 == 0 ? 0 : 1.0 / R2;
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
    int nunit = 1 << maxLevel;
#pragma omp parallel for
    for( int i=0; i<numLeafs; i++ ) {
      int ix[3] = {0, 0, 0};
      getIndex(ix,i);
      int jxmin[3];
      for_3d jxmin[d] = MAX(0,         ix[d] - numNeighbors);
      int jxmax[3];
      for_3d jxmax[d] = MIN(nunit - 1, ix[d] + numNeighbors);
      int jx[3];

      for( jx[2]=jxmin[2]; jx[2]<=jxmax[2]; jx[2]++ ) {
        for( jx[1]=jxmin[1]; jx[1]<=jxmax[1]; jx[1]++ ) {
          for( jx[0]=jxmin[0]; jx[0]<=jxmax[0]; jx[0]++ ) {
            int j = getKey(jx,maxLevel,false);

	    if(ix[0]==7 && ix[1]==7 && ix[2]==7){
	      printf("%d, %d, %d\n", jx[0], jx[1], jx[2]);
	    }

            P2PSum(Leafs[i][0],Leafs[i][1],Leafs[j][0],Leafs[j][1]);
          }
        }
      }
    }
  }


  void P2P(int *nn_stencil) const {
    //int nunit = 1 << maxLevel;
#pragma omp parallel for
    for( int i=0; i<numLeafs; i++ ) {
      int ix[3] = {0, 0, 0};
      getIndex(ix,i);
      //int jxmin[3];
      //for_3d jxmin[d] = MAX(0,         ix[d] - numNeighbors);
      //int jxmax[3];
      //for_3d jxmax[d] = MIN(nunit - 1, ix[d] + numNeighbors);
      int jx[3];

      int node_code[DIM];
      node_code[0] = ix[2]; node_code[1] = ix[1]; node_code[2] = ix[0];

      int nn_codes[DIM*NN];
      int nn_count = executeStencilNN(nn_codes, node_code, nn_stencil, 0, maxLevel);

      for(int k=0; k<nn_count; k++){
        jx[2] = nn_codes[k*DIM + 0];
        jx[1] = nn_codes[k*DIM + 1];
        jx[0] = nn_codes[k*DIM + 2];


	int j = getKey(jx,maxLevel,false);                                                                                                                           
	P2PSum(Leafs[i][0],Leafs[i][1],Leafs[j][0],Leafs[j][1]);     
      }

      /* P2P interaction within the same box */
      jx[0] = ix[0]; jx[1] = ix[1]; jx[2] = ix[2]; 
      int j = getKey(jx,maxLevel,false);
      P2PSum(Leafs[i][0],Leafs[i][1],Leafs[j][0],Leafs[j][1]);

      /*
      for( jx[2]=jxmin[2]; jx[2]<=jxmax[2]; jx[2]++ ) {
        for( jx[1]=jxmin[1]; jx[1]<=jxmax[1]; jx[1]++ ) {
          for( jx[0]=jxmin[0]; jx[0]<=jxmax[0]; jx[0]++ ) {
            int j = getKey(jx,maxLevel,false);
            P2PSum(Leafs[i][0],Leafs[i][1],Leafs[j][0],Leafs[j][1]);
          }
        }
      }
      */

    }
  }




  void P2M() const {
    int levelOffset = ((1 << 3 * maxLevel) - 1) / 7;
#pragma omp parallel for
    for( int i=0; i<numLeafs; i++ ) {
      real center[3];
      getCenter(center,i,maxLevel);
      for( int j=Leafs[i][0]; j<Leafs[i][1]; j++ ) {
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
    for( int lev=maxLevel; lev>0; lev-- ) {
      int childOffset = ((1 << 3 * lev) - 1) / 7;
      int parentOffset = ((1 << 3 * (lev - 1)) - 1) / 7;
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
        real C[LTERM];
        C[0] = 1;
        powerM(C,dist);
        for_m Multipole[p][m] += C[m] * Multipole[c][0];
        M2MSum(Multipole[p],C,Multipole[c]);
      }
    }
  }

  void M2L() const {
    for( int lev=1; lev<=maxLevel; lev++ ) {
      int levelOffset = ((1 << 3 * lev) - 1) / 7;
      int nunit = 1 << lev;
      real diameter = 2 * R0 / (1 << lev);
#pragma omp parallel for
      for( int i=0; i<(1 << 3 * lev); i++ ) {
        real L[LTERM];
        for_l L[l] = 0;
        int ix[3] = {0, 0, 0};
        getIndex(ix,i);
        int jxmin[3];
        for_3d jxmin[d] =  MAX(0,                (ix[d] >> 1) - numNeighbors) << 1;
        int jxmax[3];
        for_3d jxmax[d] = (MIN((nunit >> 1) - 1, (ix[d] >> 1) + numNeighbors) << 1) + 1;
        int jx[3];
        for( jx[2]=jxmin[2]; jx[2]<=jxmax[2]; jx[2]++ ) {
          for( jx[1]=jxmin[1]; jx[1]<=jxmax[1]; jx[1]++ ) {
            for( jx[0]=jxmin[0]; jx[0]<=jxmax[0]; jx[0]++ ) {
              if(jx[0] < ix[0]-numNeighbors || ix[0]+numNeighbors < jx[0] ||
                 jx[1] < ix[1]-numNeighbors || ix[1]+numNeighbors < jx[1] ||
                 jx[2] < ix[2]-numNeighbors || ix[2]+numNeighbors < jx[2]) {

		/*
		if(ix[0]==0 && ix[1]==0 && ix[2]==0){
		  printf("%d, %d, %d\n", jx[0], jx[1], jx[2]);
		}
		*/

                int j = getKey(jx,lev);
                real dist[3];
                for_3d dist[d] = (ix[d] - jx[d]) * diameter;
                real invR2 = 1. / (dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2]);
                real invR  = sqrt(invR2);
                real C[LTERM];
                getCoef(C,dist,invR2,invR);
                M2LSum(L,C,Multipole[j]);
              }
            }
          }
        }
        for_l Local[i+levelOffset][l] += L[l];
      }
    }
  }

  void M2L(int *common_stencil, int *far_stencil) const {
    for( int lev=1; lev<=maxLevel; lev++ ) {
      int levelOffset = ((1 << 3 * lev) - 1) / 7;
      //int nunit = 1 << lev;
      real diameter = 2 * R0 / (1 << lev);
#pragma omp parallel for
      for( int i=0; i<(1 << 3 * lev); i++ ) {
        real L[LTERM];
        for_l L[l] = 0;
        int ix[3] = {0, 0, 0};
        getIndex(ix,i);
        //int jxmin[3];
        //for_3d jxmin[d] =  MAX(0,                (ix[d] >> 1) - numNeighbors) << 1;
        //int jxmax[3];
        //for_3d jxmax[d] = (MIN((nunit >> 1) - 1, (ix[d] >> 1) + numNeighbors) << 1) + 1;
        int jx[3];

	int node_code[DIM];
	node_code[0] = ix[2]; node_code[1] = ix[1]; node_code[2] = ix[0];

	int fn_codes[DIM*(FN+CN)];
	int fn_count = executeStencilFN(fn_codes, node_code, far_stencil, common_stencil, 0, lev);


	for(int k=0; k<fn_count; k++){
	  jx[2] = fn_codes[k*DIM + 0];
	  jx[1] = fn_codes[k*DIM + 1];
	  jx[0] = fn_codes[k*DIM + 2];

	  /*
	  if(ix[0]==0 && ix[1]==0 && ix[2]==0){
	    printf("%d, %d, %d\n", jx[0], jx[1], jx[2]);
	  }
	  */

	  int j = getKey(jx,lev);                                                                                                                                  
	  real dist[3];                                                                                                                                            
	  for_3d dist[d] = (ix[d] - jx[d]) * diameter;                                                                                                             
	  real invR2 = 1. / (dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2]);                                                                           
	  real invR  = sqrt(invR2);                                                                                                                                
	  real C[LTERM];                                                                                                                                           
	  getCoef(C,dist,invR2,invR);                                                                                                                              
	  M2LSum(L,C,Multipole[j]);        

	}

	/*
        for( jx[2]=jxmin[2]; jx[2]<=jxmax[2]; jx[2]++ ) {
          for( jx[1]=jxmin[1]; jx[1]<=jxmax[1]; jx[1]++ ) {
            for( jx[0]=jxmin[0]; jx[0]<=jxmax[0]; jx[0]++ ) {
              if(jx[0] < ix[0]-numNeighbors || ix[0]+numNeighbors < jx[0] ||
                 jx[1] < ix[1]-numNeighbors || ix[1]+numNeighbors < jx[1] ||
                 jx[2] < ix[2]-numNeighbors || ix[2]+numNeighbors < jx[2]) {
                int j = getKey(jx,lev);
                real dist[3];
                for_3d dist[d] = (ix[d] - jx[d]) * diameter;
                real invR2 = 1. / (dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2]);
                real invR  = sqrt(invR2);
                real C[LTERM];
                getCoef(C,dist,invR2,invR);
                M2LSum(L,C,Multipole[j]);
              }
            }
          }
        }
	*/

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
    int levelOffset = ((1 << 3 * maxLevel) - 1) / 7;
#pragma omp parallel for
    for( int i=0; i<numLeafs; i++ ) {
      real center[3];
      getCenter(center,i,maxLevel);
      real L[LTERM];
      for_l L[l] = Local[i+levelOffset][l];
      for( int j=Leafs[i][0]; j<Leafs[i][1]; j++ ) {
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

};
