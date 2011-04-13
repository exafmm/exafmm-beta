#ifndef vortex_h
#define vortex_h
#include "fft.h"

class Vortex : public FastFourierTransform {
private:
  float dx;
  float *r, *x;
  float *dxdt, *dydt, *dzdt;

  void rbf(Bodies &bodies, int d) {
    const int itmax = 5;
    const float tol = 1e-4;
    Cells cells;

    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int i = B-bodies.begin();
      B->SRC[0] = x[i] = B->TRG[d+1] * dx * dx * dx;
      B->TRG[0] = 0;
    }
    setKernel("Gaussian");
    octsection(bodies);
    bottomup(bodies,cells);
    commBodies(cells);
    Bodies jbodies = bodies;
    Cells jcells = cells;
    commCells(jbodies,jcells);
    downward(cells,jcells,1);
    unpartition(bodies);
    std::sort(bodies.begin(),bodies.end());

    float resRecv, resSend = 0;
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int i = B-bodies.begin();
      r[i] = B->TRG[d+1] - B->TRG[0];
      B->SRC[0] = r[i];
      B->TRG[0] = 0;
      resSend += r[i] * r[i];
    }
    MPI_Allreduce(&resSend,&resRecv,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
    float res0 = resRecv;
    int it = 0;
    while( sqrt(res0) > tol && sqrt(resRecv / res0) > tol && it < itmax ) {
      print("iteration : ",0);
      print(it,0);
      print(", residual  : ",0);
      print(sqrt(resRecv / res0));
      cells.clear();
      jcells.clear();
      octsection(bodies);
      bottomup(bodies,cells);
      commBodies(cells);
      jbodies = bodies;
      jcells = cells;
      commCells(jbodies,jcells);
      downward(cells,jcells,1);
      unpartition(bodies);
      std::sort(bodies.begin(),bodies.end());
      float pApRecv, pApSend = 0;
      for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
        pApSend += B->SRC[0] * B->TRG[0];
      }
      MPI_Allreduce(&pApSend,&pApRecv,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
      float alpha = resRecv / pApRecv;
      float resOld = resRecv;
      resSend = 0;
      for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
        int i = B-bodies.begin();
        x[i] += alpha * B->SRC[0];
        r[i] -= alpha * B->TRG[0];
        resSend += r[i] * r[i];
      }
      MPI_Allreduce(&resSend,&resRecv,1,MPI_FLOAT,MPI_SUM,MPI_COMM_WORLD);
      float beta = resRecv / resOld;
      for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
        int i = B-bodies.begin();
        B->SRC[0] = r[i] + beta * B->SRC[0];
        B->TRG[0] = 0;
      }
      it++;
    }
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int i = B-bodies.begin();
      B->SRC[d] = x[i];
    }
  }

public:
  Vortex(int numGrid1D) : FastFourierTransform(numGrid1D) {
    dx   = 2 * M_PI / numGrid1D;
    r    = new float [numBodies];
    x    = new float [numBodies];
    dxdt = new float [numBodies];
    dydt = new float [numBodies];
    dzdt = new float [numBodies];
  }

  ~Vortex() {
    delete[] r;
    delete[] x;
    delete[] dxdt;
    delete[] dydt;
    delete[] dzdt;
  }

  void readData(Bodies &bodies) {                               // Initialize source values
#if 1
    char fname[256];
    sprintf(fname,"../../isotropic/spectral/initialu%4.4d",MPIRANK);
    std::ifstream fid(fname,std::ios::in);
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      fid >> B->SRC[0];
    }
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      fid >> B->SRC[1];
    }
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      fid >> B->SRC[2];
    }
#else
    std::ifstream fid("../../isotropic/spectral/initialu",std::ios::in);
    float dummy;
    for( int rank=0; rank!=MPISIZE; ++rank ) {
      if( rank == MPIRANK ) {
        for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) fid >> B->SRC[0];
      } else {
        for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) fid >> dummy;
      }
    }
    for( int rank=0; rank!=MPISIZE; ++rank ) {
      if( rank == MPIRANK ) {
        for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) fid >> B->SRC[1];
      } else {
        for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) fid >> dummy;
      }
    }
    for( int rank=0; rank!=MPISIZE; ++rank ) {
      if( rank == MPIRANK ) {
        for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) fid >> B->SRC[2];
      } else {
        for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) fid >> dummy;
      }
    }
#endif
    fid.close();

    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int i = B-bodies.begin();
      int ix = (i + numBodies * MPIRANK) / nx / nx;
      int iy = (i + numBodies * MPIRANK) / nx % nx;
      int iz = (i + numBodies * MPIRANK) % nx;
      B->IBODY = i;                                             //  Tag body with initial index
      B->IPROC = MPIRANK;                                       //  Tag body with initial MPI rank
      B->X[0] = (ix + .5) * dx - M_PI;                          //  Initialize x position
      B->X[1] = (iy + .5) * dx - M_PI;                          //  Initialize y position
      B->X[2] = (iz + .5) * dx - M_PI;                          //  Initialize z position
      B->SRC[3] = dx;                                           //  Initialize core radius
      realRecv[i] = B->SRC[1];
    }
    zDerivative();
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      B->TRG[1] = realSend[B-bodies.begin()];
      realRecv[B-bodies.begin()] = B->SRC[2];
    }
    yDerivative();
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      B->TRG[1] -= realSend[B-bodies.begin()];
      realRecv[B-bodies.begin()] = B->SRC[2];
    }
    xDerivative();
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      B->TRG[2] = realSend[B-bodies.begin()];
      realRecv[B-bodies.begin()] = B->SRC[0];
    }
    zDerivative();
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      B->TRG[2] -= realSend[B-bodies.begin()];
      realRecv[B-bodies.begin()] = B->SRC[0];
    }
    yDerivative();
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      B->TRG[3] = realSend[B-bodies.begin()];
      realRecv[B-bodies.begin()] = B->SRC[1];
    }
    xDerivative();
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      B->TRG[3] -= realSend[B-bodies.begin()];
    }
    setGlobDomain(bodies);
    rbf(bodies,2);
    rbf(bodies,1);
    rbf(bodies,0);
  }

  void initialError(Bodies &bodies) {
    Cells cells;
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      B->TRG = 0;
    }
    setKernel("BiotSavart");
    octsection(bodies);
    bottomup(bodies,cells);
    commBodies(cells);
    Bodies jbodies = bodies;
    Cells jcells = cells;
    commCells(jbodies,jcells);
    downward(cells,jcells,1);
    unpartition(bodies);
    std::sort(bodies.begin(),bodies.end());

    float u, v, w;
    double diff = 0, norm = 0;
#if 1
    char fname[256];
    sprintf(fname,"../../isotropic/spectral/initialu%4.4d",MPIRANK);
    std::ifstream fid(fname,std::ios::in);
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      fid >> u;
      diff += (B->TRG[0] - u) * (B->TRG[0] - u);
      norm += u * u;
    }
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      fid >> v;
      diff += (B->TRG[1] - v) * (B->TRG[1] - v);
      norm += v * v;
    }
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      fid >> w;
      diff += (B->TRG[2] - w) * (B->TRG[2] - w);
      norm += w * w;
    }
#else
    std::ifstream fid("../../isotropic/spectral/initialu",std::ios::in);
    float dummy;
    for( int rank=0; rank!=MPISIZE; ++rank ) {
      if( rank == MPIRANK ) {
        for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
          fid >> u;
          diff += (B->TRG[0] - u) * (B->TRG[0] - u);
          norm += u * u;
        }
      } else {
        for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) fid >> dummy;
      }
    }
    for( int rank=0; rank!=MPISIZE; ++rank ) {
      if( rank == MPIRANK ) {
        for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
          fid >> v;
          diff += (B->TRG[1] - v) * (B->TRG[1] - v);
          norm += v * v;
        }
      } else {
        for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) fid >> dummy;
      }
    }
    for( int rank=0; rank!=MPISIZE; ++rank ) {
      if( rank == MPIRANK ) {
        for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
          fid >> w;
          diff += (B->TRG[2] - w) * (B->TRG[2] - w);
          norm += w * w;
        }
      } else {
        for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) fid >> dummy;
      }
    }
#endif
    fid.close();
    print("Error     : ",0);
    print(sqrt(diff/norm));
  }

  void statistics(Bodies &bodies, bool fft=true) {
    if( fft ) {
      initSpectrum();
      for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) realRecv[B-bodies.begin()] = B->TRG[0];
      forwardFFT();
      addSpectrum();
      for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) realRecv[B-bodies.begin()] = B->TRG[1];
      forwardFFT();
      addSpectrum();
      for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) realRecv[B-bodies.begin()] = B->TRG[2];
      forwardFFT();
      addSpectrum();
      writeSpectrum();
    }
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int i = B-bodies.begin();
      dxdt[i] = B->TRG[0];
      dydt[i] = B->TRG[1];
      dzdt[i] = B->TRG[2];
    }
  }

  void convect(Bodies &bodies, float nu, float dt) {
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int i = B-bodies.begin();
      B->X[0]   += dxdt[i] * dt;
      B->X[1]   += dydt[i] * dt;
      B->X[2]   += dzdt[i] * dt;
      B->SRC[0] += B->TRG[0] * dt;
      B->SRC[1] += B->TRG[1] * dt;
      B->SRC[2] += B->TRG[2] * dt;
      B->SRC[3] += nu / B->SRC[3] * dt;
    }
  }

  void reinitialize(Bodies &bodies) {
    Cells cells, jcells;
    Bodies bodies2 = bodies, jbodies = bodies;
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int i = B-bodies.begin();
      int ix = (i + numBodies * MPIRANK) / nx / nx;
      int iy = (i + numBodies * MPIRANK) / nx % nx;
      int iz = (i + numBodies * MPIRANK) % nx;
      B->X[0] = (ix + .5) * dx - M_PI;
      B->X[1] = (iy + .5) * dx - M_PI;
      B->X[2] = (iz + .5) * dx - M_PI;
      B->TRG = 0;
    }

    setKernel("Gaussian");
    octsection(bodies);
    octsection(jbodies);
    bottomup(bodies,cells);
    bottomup(jbodies,jcells);
    commBodies(jcells);
    commCells(jbodies,jcells);
    downward(cells,jcells,1);
    unpartition(bodies);
    std::sort(bodies.begin(),bodies.end());
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      B->TRG[1] = B->TRG[0];
      B->TRG[0] = 0;
    }
    for( B_iter B=bodies2.begin(); B!=bodies2.end(); ++B ) {
      B->SRC[0] = B->SRC[1];
    }

    cells.clear();
    jcells.clear();
    jbodies = bodies2;
    octsection(bodies);
    octsection(jbodies);
    bottomup(bodies,cells);
    bottomup(jbodies,jcells);
    commBodies(jcells);
    commCells(jbodies,jcells);
    downward(cells,jcells,1);
    unpartition(bodies);
    std::sort(bodies.begin(),bodies.end());
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      B->TRG[2] = B->TRG[0];
      B->TRG[0] = 0;
    }
    for( B_iter B=bodies2.begin(); B!=bodies2.end(); ++B ) {
      B->SRC[0] = B->SRC[2];
    }

    cells.clear();
    jcells.clear();
    jbodies = bodies2;
    octsection(bodies);
    octsection(jbodies);
    bottomup(bodies,cells);
    bottomup(jbodies,jcells);
    commBodies(jcells);
    commCells(jbodies,jcells);
    downward(cells,jcells,1);
    unpartition(bodies);
    std::sort(bodies.begin(),bodies.end());
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      B->TRG[3] = B->TRG[0];
      B->SRC[3] = dx;
    }

    rbf(bodies,2);
    rbf(bodies,1);
    rbf(bodies,0);
  }
};

#endif
