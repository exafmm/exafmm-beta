#ifndef vortex_h
#define vortex_h
#include "fft.h"

class Vortex : public FastFourierTransform {
private:
  float dx;
  float *r, *x;
  float *dxdt, *dydt, *dzdt;

  void rbf(Bodies &bodies, Cells &cells, int d) {
    const int itmax = 5;
    const float tol = 1e-4;

    setKernel("Gaussian");
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int i = B-bodies.begin();
      B->SRC[0] = x[i] = B->TRG[d+1] * dx * dx * dx;
      B->TRG[0] = 0;
    }
    commBodies(cells);
    Bodies jbodies = bodies;
    Cells jcells;
    bodies2cells(jbodies,jcells);
    downward(cells,jcells,1);

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
      print("iteration     : ",0);
      print(it,0);
      print(", residual  : ",0);
      print(sqrt(resRecv / res0));
      evalP2M(cells);
      evalM2M(cells);
      updateBodies();
      jbodies = bodies;
      jcells.clear();
      bodies2cells(jbodies,jcells);
      downward(cells,jcells,1);
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

  void readData(Bodies &bodies, Cells &cells) {                 // Initialize source values
#if 1
    std::ifstream fid("../../hioki/3d/isotropic/initialuc",std::ios::in|std::ios::binary);
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int byte;
      int i = B-bodies.begin();
      int ix = (i + numBodies * MPIRANK) / nx / nx;
      int iy = (i + numBodies * MPIRANK) / nx % nx;
      int iz = (i + numBodies * MPIRANK) % nx;
      i = iz * nx * nx + iy * nx + ix;
      fid.read((char*)&byte,sizeof(int));
      fid.read((char*)&bodies[i].SRC[0],byte);
      fid.read((char*)&byte,sizeof(int));
    }
#elif 0
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
#if 1
    int *l1 = new int [nx];
    for( int i=0; i!=nx-1; ++i ) l1[i] = i+1;
    l1[nx-1] = 0;
#endif
    
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
#if 1
      int i1 =    ix  * nx * nx +    iy  * nx +    iz;
      int i2 =    ix  * nx * nx +    iy  * nx + l1[iz];
      int i3 =    ix  * nx * nx + l1[iy] * nx +    iz;
      int i4 =    ix  * nx * nx + l1[iy] * nx + l1[iz];
      int i5 = l1[ix] * nx * nx +    iy  * nx +    iz;
      int i6 = l1[ix] * nx * nx +    iy  * nx + l1[iz];
      int i7 = l1[ix] * nx * nx + l1[iy] * nx +    iz;
      int i8 = l1[ix] * nx * nx + l1[iy] * nx + l1[iz];
//      float ux = (bodies[i5].SRC[0] + bodies[i6].SRC[0] + bodies[i7].SRC[0] + bodies[i8].SRC[0]
//                - bodies[i1].SRC[0] - bodies[i2].SRC[0] - bodies[i3].SRC[0] - bodies[i4].SRC[0]) / 4 / dx;
      float uy = (bodies[i3].SRC[0] + bodies[i4].SRC[0] + bodies[i7].SRC[0] + bodies[i8].SRC[0]
                - bodies[i1].SRC[0] - bodies[i2].SRC[0] - bodies[i5].SRC[0] - bodies[i6].SRC[0]) / 4 / dx;
      float uz = (bodies[i2].SRC[0] + bodies[i4].SRC[0] + bodies[i6].SRC[0] + bodies[i8].SRC[0]
                - bodies[i1].SRC[0] - bodies[i3].SRC[0] - bodies[i5].SRC[0] - bodies[i7].SRC[0]) / 4 / dx;
      float vx = (bodies[i5].SRC[1] + bodies[i6].SRC[1] + bodies[i7].SRC[1] + bodies[i8].SRC[1]
                - bodies[i1].SRC[1] - bodies[i2].SRC[1] - bodies[i3].SRC[1] - bodies[i4].SRC[1]) / 4 / dx;
//      float vy = (bodies[i3].SRC[1] + bodies[i4].SRC[1] + bodies[i7].SRC[1] + bodies[i8].SRC[1]
//                - bodies[i1].SRC[1] - bodies[i2].SRC[1] - bodies[i5].SRC[1] - bodies[i6].SRC[1]) / 4 / dx;
      float vz = (bodies[i2].SRC[1] + bodies[i4].SRC[1] + bodies[i6].SRC[1] + bodies[i8].SRC[1]
                - bodies[i1].SRC[1] - bodies[i3].SRC[1] - bodies[i5].SRC[1] - bodies[i7].SRC[1]) / 4 / dx;
      float wx = (bodies[i5].SRC[2] + bodies[i6].SRC[2] + bodies[i7].SRC[2] + bodies[i8].SRC[2]
                - bodies[i1].SRC[2] - bodies[i2].SRC[2] - bodies[i3].SRC[2] - bodies[i4].SRC[2]) / 4 / dx;
      float wy = (bodies[i3].SRC[2] + bodies[i4].SRC[2] + bodies[i7].SRC[2] + bodies[i8].SRC[2]
                - bodies[i1].SRC[2] - bodies[i2].SRC[2] - bodies[i5].SRC[2] - bodies[i6].SRC[2]) / 4 / dx;
//      float wz = (bodies[i2].SRC[2] + bodies[i4].SRC[2] + bodies[i6].SRC[2] + bodies[i8].SRC[2]
//                - bodies[i1].SRC[2] - bodies[i3].SRC[2] - bodies[i5].SRC[2] - bodies[i7].SRC[2]) / 4 / dx;
      B->TRG[1] = vz - wy;
      B->TRG[2] = wx - uz;
      B->TRG[3] = uy - vx;
#endif
    }
#if 0
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
#endif
    setGlobDomain(bodies);
    octsection(bodies);
    bottomup(bodies,cells);
    rbf(bodies,cells,2);
    rbf(bodies,cells,1);
    rbf(bodies,cells,0);
    unpartition(bodies);
    std::sort(bodies.begin(),bodies.end());
  }

  void gridVelocity(Bodies &bodies, Cells &cells) {
    Bodies jbodies = bodies;
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int i = B-bodies.begin();
      int ix = (i + numBodies * MPIRANK) / nx / nx;
      int iy = (i + numBodies * MPIRANK) / nx % nx;
      int iz = (i + numBodies * MPIRANK) % nx;
      B->IBODY = i;                                             //  Tag body with initial index
      B->IPROC = MPIRANK;                                       //  Tag body with initial MPI rank
      B->X[0] = ix * dx - M_PI + 1e-5;
      B->X[1] = iy * dx - M_PI + 1e-5;
      B->X[2] = iz * dx - M_PI + 1e-5;
    }

    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      B->TRG = 0;
    }
    setKernel("BiotSavart");
    cells.clear();
    Cells jcells;
    octsection(bodies);
    octsection(jbodies);
    bottomup(bodies,cells);
    bottomup(jbodies,jcells);
    commBodies(jcells);
    commCells(jbodies,jcells);
    downward(cells,jcells,1);
    unpartition(bodies);
    std::sort(bodies.begin(),bodies.end());
  }

  void initialError(Bodies &bodies) {
    float u, v, w;
    double diff = 0, norm = 0;
#if 1
    std::ifstream fid("../../hioki/3d/isotropic/initialuc",std::ios::in|std::ios::binary);
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int byte;
      int i = B-bodies.begin();
      int ix = (i + numBodies * MPIRANK) / nx / nx;
      int iy = (i + numBodies * MPIRANK) / nx % nx;
      int iz = (i + numBodies * MPIRANK) % nx;
      i = iz * nx * nx + iy * nx + ix;
      fid.read((char*)&byte,sizeof(int));
      fid.read((char*)&u,sizeof(float));
      fid.read((char*)&v,sizeof(float));
      fid.read((char*)&w,sizeof(float));
      fid.read((char*)&byte,sizeof(int));
      diff += (bodies[i].TRG[0] - u) * (bodies[i].TRG[0] - u)
            + (bodies[i].TRG[1] - v) * (bodies[i].TRG[1] - v)
            + (bodies[i].TRG[2] - w) * (bodies[i].TRG[2] - w);
      norm += u * u + v * v + w * w;
    }
    fid.close();
#elif 0
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
    print("Error         : ",0);
    print(sqrt(diff/norm));
  }

  void statistics(Bodies &bodies, float nu, float dt) {
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

    float umax = 0;
    double uu = 0, vv = 0, ww = 0, uw = 0, ds = 0;
    double ux2 = 0, ux3 = 0, ux4 = 0;
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      umax = std::max(umax,std::abs(B->TRG[0]));
      umax = std::max(umax,std::abs(B->TRG[1]));
      umax = std::max(umax,std::abs(B->TRG[2]));
      uu += B->TRG[0] * B->TRG[0];
      vv += B->TRG[1] * B->TRG[1];
      ww += B->TRG[2] * B->TRG[2];
      uw += B->TRG[0] * B->TRG[2];
      realRecv[B-bodies.begin()] = B->TRG[0];
    }
    xDerivative();
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      float ux = realSend[B-bodies.begin()];
      ds += ux * ux;
      ux2 += ux * ux;
      ux3 += ux * ux * ux;
      ux4 += ux * ux * ux * ux;
      realRecv[B-bodies.begin()] = B->TRG[0];
    }
    yDerivative();
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      float uy = realSend[B-bodies.begin()];
      ds += uy * uy;
      realRecv[B-bodies.begin()] = B->TRG[0];
    }
    zDerivative();
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      float uz = realSend[B-bodies.begin()];
      ds += uz * uz;
      realRecv[B-bodies.begin()] = B->TRG[1];
    }
    xDerivative();
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      float vx = realSend[B-bodies.begin()];
      ds += vx * vx;
      realRecv[B-bodies.begin()] = B->TRG[1];
    }
    yDerivative();
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      float vy = realSend[B-bodies.begin()];
      ds += vy * vy;
      realRecv[B-bodies.begin()] = B->TRG[1];
    }
    zDerivative();
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      float vz = realSend[B-bodies.begin()];
      ds += vz * vz;
      realRecv[B-bodies.begin()] = B->TRG[2];
    }
    xDerivative();
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      float wx = realSend[B-bodies.begin()];
      ds += wx * wx;
      realRecv[B-bodies.begin()] = B->TRG[2];
    }
    yDerivative();
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      float wy = realSend[B-bodies.begin()];
      ds += wy * wy;
      realRecv[B-bodies.begin()] = B->TRG[2];
    }
    zDerivative();
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      float wz = realSend[B-bodies.begin()];
      ds += wz * wz;
    }
    float ek = uu + vv + ww;
    uu /= ek;
    vv /= ek;
    ww /= ek;
    uw /= ek;
    ek /= 2 * bodies.size();
    ds *= nu / bodies.size() / 4;
    ux2 /= bodies.size();
    ux3 /= bodies.size();
    ux4 /= bodies.size();
    float sk = ux3 / std::pow(ux2,1.5);
    float fl = ux4 / ux2 / ux2;
    float ret = ek * ek / nu / ds;
    float rel = sqrt(20 * ret / 3);
    float cfl = umax * dt / dx;
    if( MPIRANK == 0 ) {
      std::ofstream fid("statistics.dat",std::ios::out | std::ios::app);
      fid << ek << std::endl;
      fid << ds << std::endl;
      fid << sk << std::endl;
      fid << fl << std::endl;
      fid << cfl << std::endl;
      fid << ret << std::endl;
      fid << rel << std::endl;
      fid << uu << std::endl;
      fid << vv << std::endl;
      fid << ww << std::endl;
      fid << uw << std::endl;
      fid.close();
      std::cout << "-------------------------------\n";
      std::cout << "| energy      : " << ek << std::endl;
      std::cout << "| dissipation : " << ds << std::endl;
      std::cout << "| skewness    : " << sk << std::endl;
      std::cout << "| flatness    : " << fl << std::endl;
      std::cout << "| CFL         : " << cfl << std::endl;
      std::cout << "| Re T        : " << ret << std::endl;
      std::cout << "| Re lambda   : " << rel << std::endl;
      std::cout << "| uu          : " << uu << std::endl;
      std::cout << "| vv          : " << vv << std::endl;
      std::cout << "| ww          : " << ww << std::endl;
      std::cout << "| uw          : " << uw << std::endl;
      std::cout << "-------------------------------\n";
    }
  }

  void BiotSavart(Bodies &bodies, Cells &cells) {
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      B->TRG = 0;
    }
    setKernel("BiotSavart");
    cells.clear();
    octsection(bodies);
    bottomup(bodies,cells);
    commBodies(cells);
    Bodies jbodies = bodies;
    Cells jcells = cells;
    commCells(jbodies,jcells);
    downward(cells,jcells,1);
    unpartition(bodies);
    std::sort(bodies.begin(),bodies.end());
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int i = B-bodies.begin();
      dxdt[i] = B->TRG[0];
      dydt[i] = B->TRG[1];
      dzdt[i] = B->TRG[2];
    }
  }

  void Stretching(Bodies &bodies, Cells &cells) {
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      B->TRG = 0;
    }
    setKernel("Stretching");
    cells.clear();
    octsection(bodies);
    bottomup(bodies,cells);
    commBodies(cells);
    Bodies jbodies = bodies;
    Cells jcells = cells;
    commCells(jbodies,jcells);
    downward(cells,jcells,1);
    unpartition(bodies);
    std::sort(bodies.begin(),bodies.end());
  }

  void update(Bodies &bodies, float nu, float dt) {
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

  void reinitialize(Bodies &bodies, Bodies &bodies2) {
    Cells cells, jcells;
    Bodies jbodies = bodies;
    bodies2 = bodies;
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

    cells.clear();
    octsection(bodies);
    bottomup(bodies,cells);
    rbf(bodies,cells,2);
    rbf(bodies,cells,1);
    rbf(bodies,cells,0);
    unpartition(bodies);
    std::sort(bodies.begin(),bodies.end());
  }
};

#endif
