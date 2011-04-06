#ifndef vortex_h
#define vortex_h
#include "fft.h"
#include "let.h"

class Vortex : public LocalEssentialTree, public FastFourierTransform {
private:
  float dx;
  float *r, *x;
  float *dxdt, *dydt, *dzdt;

  void rbf(Bodies &bodies, int d) {
    const int itmax = 5;
    const float tol = 1e-4;
    Cells cells, jcells;

    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int i = B-bodies.begin();
      B->SRC[0] = x[i] = B->TRG[d+1] * dx * dx * dx;
      B->TRG[0] = 0;
    }
    setKernel("Gaussian");
    setDomain(bodies);
    bottomup(bodies,cells);
    jcells = cells;
    downward(cells,jcells,1);
    std::sort(bodies.begin(),bodies.end());

    float res = 0;
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int i = B-bodies.begin();
      r[i] = B->TRG[d+1] - B->TRG[0];
      B->SRC[0] = r[i];
      B->TRG[0] = 0;
      res += r[i] * r[i];
    }
    float res0 = res;
    int it = 0;
    while( sqrt(res0) > tol && sqrt(res / res0) > tol && it < itmax ) {
      std::cout << "iteration : " << it << ", residual : " << sqrt(res / res0) << std::endl;
      cells.clear();
      bottomup(bodies,cells);
      jcells = cells;
      downward(cells,jcells,1);
      std::sort(bodies.begin(),bodies.end());
      float pAp = 0;
      for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
        pAp += B->SRC[0] * B->TRG[0];
      }
      float alpha = res / pAp;
      float resOld = res;
      res = 0;
      for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
        int i = B-bodies.begin();
        x[i] += alpha * B->SRC[0];
        r[i] -= alpha * B->TRG[0];
        res += r[i] * r[i];
      }
      float beta = res / resOld;
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
    std::ifstream fid("../../isotropic/spectral/initialu",std::ios::in);
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      fid >> B->SRC[0];
    }
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      fid >> B->SRC[1];
    }
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      fid >> B->SRC[2];
    }
    fid.close();

    int *plus1 = new int [nx];
    for( int i=0; i!=nx-1; ++i ) {
      plus1[i] = i+1;
    }
    plus1[nx-1] = 0;

    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      int i = B-bodies.begin();
      int ix = i / nx / nx;
      int iy = i / nx % nx;
      int iz = i % nx;
      Body B0 = bodies[ix        * nx * nx + iy        * nx + iz       ];
      Body B1 = bodies[ix        * nx * nx + iy        * nx + plus1[iz]];
      Body B2 = bodies[ix        * nx * nx + plus1[iy] * nx + iz       ];
      Body B3 = bodies[ix        * nx * nx + plus1[iy] * nx + plus1[iz]];
      Body B4 = bodies[plus1[ix] * nx * nx + iy        * nx + iz       ];
      Body B5 = bodies[plus1[ix] * nx * nx + iy        * nx + plus1[iz]];
      Body B6 = bodies[plus1[ix] * nx * nx + plus1[iy] * nx + iz       ];
      Body B7 = bodies[plus1[ix] * nx * nx + plus1[iy] * nx + plus1[iz]];
      B->IBODY = i;                                             //  Tag body with initial index
      B->IPROC = RANK;                                          //  Tag body with initial MPI rank
      B->X[0] = (ix + .5) * dx - M_PI;                          //  Initialize x position
      B->X[1] = (iy + .5) * dx - M_PI;                          //  Initialize y position
      B->X[2] = (iz + .5) * dx - M_PI;                          //  Initialize z position
      float uy = B2.SRC[0] + B3.SRC[0] + B6.SRC[0] + B7.SRC[0]
               - B0.SRC[0] - B1.SRC[0] - B4.SRC[0] - B5.SRC[0];
      float uz = B1.SRC[0] + B3.SRC[0] + B5.SRC[0] + B7.SRC[0]
               - B0.SRC[0] - B2.SRC[0] - B4.SRC[0] - B6.SRC[0];
      float vx = B4.SRC[1] + B5.SRC[1] + B6.SRC[1] + B7.SRC[1]
               - B0.SRC[1] - B1.SRC[1] - B2.SRC[1] - B3.SRC[1];
      float vz = B1.SRC[1] + B3.SRC[1] + B5.SRC[1] + B7.SRC[1]
               - B0.SRC[1] - B2.SRC[1] - B4.SRC[1] - B6.SRC[1];
      float wx = B4.SRC[2] + B5.SRC[2] + B6.SRC[2] + B7.SRC[2]
               - B0.SRC[2] - B1.SRC[2] - B2.SRC[2] - B3.SRC[2];
      float wy = B2.SRC[2] + B3.SRC[2] + B6.SRC[2] + B7.SRC[2]
               - B0.SRC[2] - B1.SRC[2] - B4.SRC[2] - B5.SRC[2];
      B->TRG[1] = (vz - wy) / 4 / dx;                           //  Initialize x vorticity
      B->TRG[2] = (wx - uz) / 4 / dx;                           //  Initialize y vorticity
      B->TRG[3] = (uy - vx) / 4 / dx;                           //  Initialize z vorticity
      B->SRC[3] = dx;                                           //  Initialize core radius
    }                                                           // End loop over bodies
    rbf(bodies,2);
    rbf(bodies,1);
    rbf(bodies,0);
    delete[] plus1;
  }

  void initialError(Bodies &bodies) {
    Bodies jbodies = bodies;
    Cells cells, jcells;
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int i = B-bodies.begin();
      int ix = i / nx / nx;
      int iy = i / nx % nx;
      int iz = i % nx;
      B->X[0] = ix * dx - M_PI;
      B->X[1] = iy * dx - M_PI;
      B->X[2] = iz * dx - M_PI;
      B->TRG = 0;
    }
    setKernel("BiotSavart");
    setDomain(bodies);
    bottomup(bodies,cells);
    bottomup(jbodies,jcells);
    downward(cells,jcells,1);
    std::sort(bodies.begin(),bodies.end());
    std::sort(jbodies.begin(),jbodies.end());

    float u, v, w;
    double diff = 0, norm = 0;
    std::ifstream fid("../../isotropic/spectral/initialu",std::ios::in);
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
    fid.close();
    std::cout << "Error : " << std::sqrt(diff/norm) << std::endl;
    bodies = jbodies;
  }

  void statistics(Bodies &bodies, bool fft=true) {
    Bodies jbodies = bodies;
    Cells cells, jcells;
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int i = B-bodies.begin();
      int ix = i / nx / nx;
      int iy = i / nx % nx;
      int iz = i % nx;
      B->X[0] = ix * dx - M_PI;
      B->X[1] = iy * dx - M_PI;
      B->X[2] = iz * dx - M_PI;
      B->TRG = 0;
    }
    setKernel("BiotSavart");
    setDomain(bodies);
    bottomup(bodies,cells);
    bottomup(jbodies,jcells);
    downward(cells,jcells,1);
    std::sort(bodies.begin(),bodies.end());
    std::sort(jbodies.begin(),jbodies.end());

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

    bodies = jbodies;
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      B->TRG = 0;
    }
    cells.clear();
    bottomup(bodies,cells);
    jcells = cells;
    downward(cells,jcells,1);
    std::sort(bodies.begin(),bodies.end());
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
    Bodies jbodies = bodies;
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int i = B-bodies.begin();
      int ix = i / nx / nx;
      int iy = i / nx % nx;
      int iz = i % nx;
      B->X[0] = (ix + .5) * dx - M_PI;
      B->X[1] = (iy + .5) * dx - M_PI;
      B->X[2] = (iz + .5) * dx - M_PI;
      B->TRG = 0;
    }

    setKernel("Gaussian");
    setDomain(bodies);
    bottomup(bodies,cells);
    bottomup(jbodies,jcells);
    downward(cells,jcells,1);
    std::sort(bodies.begin(),bodies.end());
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      B->TRG[1] = B->TRG[0];
      B->TRG[0] = 0;
    }
    for( B_iter B=jbodies.begin(); B!=jbodies.end(); ++B ) {
      B->SRC[0] = B->SRC[1];
    }

    cells.clear();
    jcells.clear();
    bottomup(bodies,cells);
    bottomup(jbodies,jcells);
    downward(cells,jcells,1);
    std::sort(bodies.begin(),bodies.end());
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      B->TRG[2] = B->TRG[0];
      B->TRG[0] = 0;
    }
    for( B_iter B=jbodies.begin(); B!=jbodies.end(); ++B ) {
      B->SRC[0] = B->SRC[2];
    }

    cells.clear();
    jcells.clear();
    bottomup(bodies,cells);
    bottomup(jbodies,jcells);
    downward(cells,jcells,1);
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
