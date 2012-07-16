#ifndef bottomup_h
#define bottomup_h
#include "topdown.h"

class BottomUp : public TopDown {
private:
  int getMaxLevel(Bodies &bodies) {
    const long N = bodies.size();
    int level;
    level = N >= NCRIT ? 1 + int(log(N / NCRIT)/M_LN2/3) : 0;
    return level;
  }

  inline void getIndex(Bodies &bodies) {
    float d = 2 * R0 / (1 << MAXLEVEL);
#if MTHREADS
#else
#pragma omp parallel for num_threads(OMP_NUM_THREADS)
#endif
    for( uint b=0; b<bodies.size(); b++ ) {
      B_iter B = bodies.begin() + b;
      int ix = int((B->X[0] + R0 - X0[0]) / d);
      int iy = int((B->X[1] + R0 - X0[1]) / d);
      int iz = int((B->X[2] + R0 - X0[2]) / d);
      int id = 0;
      for( int l=0; l!=MAXLEVEL; ++l ) {
        id += (ix & 1) << (3 * l);
        id += (iy & 1) << (3 * l + 1);
        id += (iz & 1) << (3 * l + 2);
        ix >>= 1;
        iy >>= 1;
        iz >>= 1;
      }
      B->ICELL = id;
    }
  }

  void radixSort(Bodies &bodies, int **index) {
    const int n = bodies.size();
    const int bitStride = 8;
    const int stride = 1 << bitStride;
    const int mask = stride - 1;
    int (*bucket2D)[stride] = new int [OMP_NUM_THREADS][stride]();
    int aMaxPerThread[OMP_NUM_THREADS] = {0};
#if MTHREADS
#else
#pragma omp parallel for num_threads(OMP_NUM_THREADS)
#endif
    for( int b=0; b<n; b++ ) {
      int i = bodies[b].ICELL;
      index[0][b] = i;
      index[2][b] = b;
      index[4][b] = i;
      if( i > aMaxPerThread[omp_get_thread_num()] )
        aMaxPerThread[omp_get_thread_num()] = i;
    }
    int aMax = 0;
    for( int i=0; i<OMP_NUM_THREADS; i++ )
      if( aMaxPerThread[i] > aMax ) aMax = aMaxPerThread[i];
    while( aMax > 0 ) {
      int bucket2[stride] = {0};
      for( int t=0; t<OMP_NUM_THREADS; t++ )
        for( int i=0; i<stride; i++ )
          bucket2D[t][i] = 0;
#if MTHREADS
#else
#pragma omp parallel for num_threads(OMP_NUM_THREADS)
#endif
      for( int i=0; i<n; i++ )
        bucket2D[omp_get_thread_num()][index[0][i] & mask]++;
      for( int t=0; t<OMP_NUM_THREADS; t++ )
        for( int i=0; i<stride; i++ )
          bucket2[i] += bucket2D[t][i];
      for( int i=1; i<stride; i++ )
        bucket2[i] += bucket2[i-1];
      for( int i=n-1; i>=0; i-- )
        index[3][i] = --bucket2[index[0][i] & mask];
#if MTHREADS
#else
#pragma omp parallel for num_threads(OMP_NUM_THREADS)
#endif
      for( int i=0; i<n; i++ )
        index[1][index[3][i]] = index[2][i];
#if MTHREADS
#else
#pragma omp parallel for num_threads(OMP_NUM_THREADS)
#endif
      for( int i=0; i<n; i++ )
        index[2][i] = index[1][i];
#if MTHREADS
#else
#pragma omp parallel for num_threads(OMP_NUM_THREADS)
#endif
      for( int i=0; i<n; i++ )
        index[1][index[3][i]] = index[0][i];
#if MTHREADS
#else
#pragma omp parallel for num_threads(OMP_NUM_THREADS)
#endif
      for( int i=0; i<n; i++ )
        index[0][i] = index[1][i] >> bitStride;
      aMax >>= bitStride;
    }
    Bodies buffer = bodies;
#if MTHREADS
#else
#pragma omp parallel for num_threads(OMP_NUM_THREADS)
#endif
    for( int b=0; b<n; b++ )
      bodies[b] = buffer[index[2][b]];
    delete[] bucket2D;
  }

  inline void initCell(Cell &cell, int child, B_iter LEAF, real diameter) {
    cell.NCHILD = 0;
    cell.NCLEAF = 0;
    cell.NDLEAF = 0;
    cell.CHILD  = child;
    cell.LEAF   = LEAF;
    int ix = int((LEAF->X[0] + R0 - X0[0]) / diameter);
    int iy = int((LEAF->X[1] + R0 - X0[1]) / diameter);
    int iz = int((LEAF->X[2] + R0 - X0[2]) / diameter);
    cell.X[0]   = diameter * (ix + .5) + X0[0] - R0;
    cell.X[1]   = diameter * (iy + .5) + X0[1] - R0;
    cell.X[2]   = diameter * (iz + .5) + X0[2] - R0;
    cell.R      = diameter * .5;
  }

  void buildBottom(Bodies &bodies, Cells &cells) {
    int I = -1;
    C_iter C;
    cells.clear();
    cells.reserve(1 << (3 * MAXLEVEL));
    float d = 2 * R0 / (1 << MAXLEVEL);
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int IC = B->ICELL;
      if( IC != I ) {
        Cell cell;
        initCell(cell,0,B,d);
        cells.push_back(cell);
        C = cells.end()-1;
        I = IC;
      }
      C->NCLEAF++;
      C->NDLEAF++;
    }
  }

  void twigs2cells(Cells &cells) {
    int begin = 0, end = cells.size();
    float d = 2 * R0 / (1 << MAXLEVEL);
    for( int l=0; l!=MAXLEVEL; ++l ) {
      int div = (8 << (3 * l));
      int I = -1;
      int p = end - 1;
      d *= 2;
      for( int c=begin; c!=end; ++c ) {
        B_iter B = cells[c].LEAF;
        int IC = B->ICELL / div;
        if( IC != I ) {
          Cell cell;
          initCell(cell,c,cells[c].LEAF,d);
          cells.push_back(cell);
          p++;
          I = IC;
        }
        cells[p].NCHILD++;
        cells[p].NDLEAF += cells[c].NDLEAF;
        cells[c].PARENT = p;
      }
      begin = end;
      end = cells.size();
    }
  }

protected:
  void setDomain(Bodies &bodies) {
    startTimer("Set domain");
    MAXLEVEL = getMaxLevel(bodies);
    vect xmin, xmax;
    X0 = 0;
    xmax = xmin = bodies.begin()->X;
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      for( int d=0; d!=3; ++d ) {
        if     (B->X[d] < xmin[d]) xmin[d] = B->X[d];
        else if(B->X[d] > xmax[d]) xmax[d] = B->X[d];
      }
      X0 += B->X;
    }
    if( IMAGES != 0 ) {
      X0 = 0;
      R0 = M_PI;
    } else {
      X0 /= bodies.size();
      for( int d=0; d!=3; ++d ) {
        X0[d] = int(X0[d]+.5);
        R0 = std::max(xmax[d] - X0[d], R0);
        R0 = std::max(X0[d] - xmin[d], R0);
      }
      R0 *= 1.000001;
    }
    stopTimer("Set domain",printNow);
  }

  void buildTree(Bodies &bodies, Cells &cells) {
    startTimer("Morton index");
    getIndex(bodies);
    stopTimer("Morton index",printNow);
    startTimer("Sort bodies");
    int **index = new int* [5];
    for( int i=0; i<5; i++ ) index[i] = new int [bodies.size()];
    radixSort(bodies,index);
    for( int i=0; i<5; i++ ) delete[] index[i];
    delete[] index;
    stopTimer("Sort bodies",printNow);
    startTimer("Build bottom");
    buildBottom(bodies,cells);
    stopTimer("Build bottom",printNow);
  }

  void linkTree(Cells &cells) {
    startTimer("Link tree");
    twigs2cells(cells);
    stopTimer("Link tree",printNow);
  }

  void upwardPass(Cells &cells) {
    startTimer("Upward pass");
    setRootCell(cells);
    for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {
      C->M = 0;
      C->L = 0;
    }
    for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {
      real Rmax = 0;
      setCenter(C);
      P2M(C,Rmax);
      M2M(C,Rmax);
    }
#if Cartesian
    for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {
      for( int i=1; i<MTERM; ++i ) C->M[i] /= C->M[0];
    }
#endif
    setRcrit(cells);
    stopTimer("Upward pass",printNow);
  }

  void downwardPass(Cells &cells) const {
    C_iter C1 = cells.end() - 1;
    L2P(C1);
    for( C_iter C=C1-1; C!=cells.begin()-1; --C ) {
      L2L(C);
      L2P(C);
    }
  }
};

#endif
