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
    const float tol = 1e-3;

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
      print(sqrt(resRecv / res0),0);
      print("\n",0);
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
    dxdt = new float [2*numBodies];
    dydt = new float [2*numBodies];
    dzdt = new float [2*numBodies];
  }

  ~Vortex() {
    delete[] r;
    delete[] x;
    delete[] dxdt;
    delete[] dydt;
    delete[] dzdt;
  }

  void readData(Bodies &bodies, Bodies &bodies2, Cells &cells) {// Initialize source values
    char fname[256];
#if 0
    sprintf(fname,"/work0/t2g-ppc-all/11ITA070/initial%4.4d",nx);
    std::ifstream fid(fname,std::ios::in|std::ios::binary);
    int byte;
    float dummy[3];
    for( int rank=0; rank!=MPISIZE; ++rank ) {
      if( rank == MPIRANK ) {
        for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
          fid.read((char*)&byte,sizeof(int));
          fid.read((char*)&bodies[B-bodies.begin()].SRC[0],byte);
          fid.read((char*)&byte,sizeof(int));
        }
      } else {
        for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
          fid.read((char*)&byte,sizeof(int));
          fid.read((char*)dummy,byte);
          fid.read((char*)&byte,sizeof(int));
        }
      }
    }
#else
    sprintf(fname,"/work0/t2g-ppc-all/11ITA070/initial%5.5d-%5.5d",nx,MPIRANK);
    std::ifstream fid(fname,std::ios::in|std::ios::binary);
    int byte;
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      fid.read((char*)&byte,sizeof(int));
      fid.read((char*)&bodies[B-bodies.begin()].SRC[0],byte);
      fid.read((char*)&byte,sizeof(int));
    }
#endif
    fid.close();
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int i = B-bodies.begin();
      int iGlob = i + numBodies * MPIRANK;
      int ix = iGlob % nx;
      int iy = iGlob / nx % nx;
      int iz = iGlob / nx / nx;
      B->IBODY = iGlob;                                         //  Tag body with initial index
      B->IPROC = MPIRANK;                                       //  Tag body with initial MPI rank
      B->X[0] = (ix + .5) * dx - M_PI;                          //  Initialize x position
      B->X[1] = (iy + .5) * dx - M_PI;                          //  Initialize y position
      B->X[2] = (iz + .5) * dx - M_PI;                          //  Initialize z position
      B->SRC[3] = dx;                                           //  Initialize core radius
      realRecv[i] = B->SRC[1];
    }

    for( B_iter B=bodies2.begin(); B!=bodies2.end(); ++B ) {
      int i = B-bodies2.begin();
      int iGlob = i + numBodies * MPIRANK;
      int ix = iGlob % nx;
      int iy = iGlob / nx % nx;
      int iz = iGlob / nx / nx;
      B->IBODY = iGlob;                                         //  Tag body with initial index
      B->IPROC = MPIRANK;                                       //  Tag body with initial MPI rank
      B->X[0] = ix * dx - M_PI + 1e-5;
      B->X[1] = iy * dx - M_PI + 1e-5;
      B->X[2] = iz * dx - M_PI + 1e-5;
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
    octsection(bodies);
    octsection(bodies2);
    bottomup(bodies,cells);
    rbf(bodies,cells,2);
    rbf(bodies,cells,1);
    rbf(bodies,cells,0);
  }

  void gridVelocity(Bodies &bodies, Bodies &bodies2, Cells &cells) {
    Bodies jbodies = bodies;
    for( B_iter B=bodies2.begin(); B!=bodies2.end(); ++B ) B->TRG = 0;
    setKernel("BiotSavart");
    cells.clear();
    Cells jcells;
    bottomup(bodies2,cells);
    bottomup(jbodies,jcells);
    commBodies(jcells);
    commCells(jbodies,jcells);
    downward(cells,jcells,1);
  }

  void initialError(Bodies bodies) {
    char fname[256];
    int byte;
    float dummy[3];
    float u, v, w;
    float diffSend = 0, normSend = 0, diffRecv, normRecv;
    unpartition(bodies);
    std::sort(bodies.begin(),bodies.end());
#if 0
    sprintf(fname,"/work0/t2g-ppc-all/11ITA070/initial%4.4d",nx);
    std::ifstream fid(fname,std::ios::in|std::ios::binary);
    for( int rank=0; rank!=MPISIZE; ++rank ) {
      if( rank == MPIRANK ) {
        for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
          int i = B-bodies.begin();
          fid.read((char*)&byte,sizeof(int));
          fid.read((char*)&u,sizeof(float));
          fid.read((char*)&v,sizeof(float));
          fid.read((char*)&w,sizeof(float));
          fid.read((char*)&byte,sizeof(int));
          diffSend += (bodies[i].TRG[0] - u) * (bodies[i].TRG[0] - u)
                    + (bodies[i].TRG[1] - v) * (bodies[i].TRG[1] - v)
                    + (bodies[i].TRG[2] - w) * (bodies[i].TRG[2] - w);
          normSend += u * u + v * v + w * w;
        }
      } else {
        for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
          fid.read((char*)&byte,sizeof(int));
          fid.read((char*)dummy,byte);
          fid.read((char*)&byte,sizeof(int));
        }
      }
    }
#else
    sprintf(fname,"/work0/t2g-ppc-all/11ITA070/initial%5.5d-%5.5d",nx,MPIRANK);
    std::ifstream fid(fname,std::ios::in|std::ios::binary);
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int i = B-bodies.begin();
      fid.read((char*)&byte,sizeof(int));
      fid.read((char*)&u,sizeof(float));
      fid.read((char*)&v,sizeof(float));
      fid.read((char*)&w,sizeof(float));
      fid.read((char*)&byte,sizeof(int));
      diffSend += (bodies[i].TRG[0] - u) * (bodies[i].TRG[0] - u)
                + (bodies[i].TRG[1] - v) * (bodies[i].TRG[1] - v)
                + (bodies[i].TRG[2] - w) * (bodies[i].TRG[2] - w);
      normSend += u * u + v * v + w * w;
    }
#endif
    fid.close();
    MPI_Reduce(&diffSend,&diffRecv,1,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&normSend,&normRecv,1,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
    print("Error         : ",0);
    print(sqrt(diffRecv/normRecv),0);
    print("\n",0);
  }

  void statistics(Bodies bodies, float nu, float dt) {
    unpartition(bodies);
    std::sort(bodies.begin(),bodies.end());
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
    float uu = 0, vv = 0, ww = 0, uw = 0, ds = 0;
    float ux2 = 0, ux3 = 0, ux4 = 0;
    float statSend[9], statRecv[9];
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
    statSend[0] = umax;
    statSend[1] = uu;
    statSend[2] = vv;
    statSend[3] = ww;
    statSend[4] = uw;
    statSend[5] = ds;
    statSend[6] = ux2;
    statSend[7] = ux3;
    statSend[8] = ux4;
    MPI_Reduce(&statSend,&statRecv,9,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
    if( MPIRANK == 0 ) {
      umax = statRecv[0] / MPISIZE;
      uu   = statRecv[1];
      vv   = statRecv[2];
      ww   = statRecv[3];
      uw   = statRecv[4];
      ds   = statRecv[5];
      ux2  = statRecv[6];
      ux3  = statRecv[7];
      ux4  = statRecv[8];
      float ek = uu + vv + ww;
      uu /= ek;
      vv /= ek;
      ww /= ek;
      uw /= ek;
      ek /= 2 * numGlobal;
      ds *= nu / numGlobal / 4;
      ux2 /= numGlobal;
      ux3 /= numGlobal;
      ux4 /= numGlobal;
      float sk = ux3 / std::pow(ux2,1.5);
      float fl = ux4 / ux2 / ux2;
      float ret = ek * ek / nu / ds;
      float rel = sqrt(20 * ret / 3);
      float cfl = umax * dt / dx;
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
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) B->TRG = 0;
    setKernel("BiotSavart");
    cells.clear();
    octsection(bodies);
    bottomup(bodies,cells);
    commBodies(cells);
    Bodies jbodies = bodies;
    Cells jcells = cells;
    commCells(jbodies,jcells);
    downward(cells,jcells,1);
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int i = B-bodies.begin();
      dxdt[i] = B->TRG[0];
      dydt[i] = B->TRG[1];
      dzdt[i] = B->TRG[2];
    }
  }

  void Stretching(Bodies &bodies, Cells &cells) {
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) B->TRG = 0;
    setKernel("Stretching");
    evalP2M(cells);
    evalM2M(cells);
    updateBodies();
    Bodies jbodies = bodies;
    Cells jcells = cells;
    commCells(jbodies,jcells);
    downward(cells,jcells,1);
  }

  void update(Bodies &bodies, float nu, float dt) {
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int i = B-bodies.begin();
      B->X[0] += dxdt[i] * dt;
      B->X[1] += dydt[i] * dt;
      B->X[2] += dzdt[i] * dt;
      for( int d=0; d!=3; ++d ) {
        if( B->X[d] < -M_PI ) {
          B->X[d] += 2 * M_PI;
        } else if( M_PI < B->X[d] ) {
          B->X[d] -= 2 * M_PI;
        }
      }
      B->SRC[0] += B->TRG[0] * dt;
      B->SRC[1] += B->TRG[1] * dt;
      B->SRC[2] += B->TRG[2] * dt;
      B->SRC[3] += nu / B->SRC[3] * dt;
    }
  }

  void reinitialize(Bodies &bodies, Cells &cells) {
    Bodies jbodies = bodies;
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      int ix = B->IBODY % nx;
      int iy = B->IBODY / nx % nx;
      int iz = B->IBODY / nx / nx;
      B->X[0] = (ix + .5) * dx - M_PI;
      B->X[1] = (iy + .5) * dx - M_PI;
      B->X[2] = (iz + .5) * dx - M_PI;
      B->TRG = 0;
    }

    setKernel("Gaussian");
    cells.clear();
    Cells jcells;
    octsection(bodies);
    bottomup(bodies,cells);
    bottomup(jbodies,jcells);
    commBodies(jcells);
    commCells(jbodies,jcells);
    int numCells = jcells.size();
    downward(cells,jcells,1);
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      B->TRG[1] = B->TRG[0];
      B->TRG[0] = 0;
    }

    jcells.resize(numCells);
    for( B_iter B=jbodies.begin(); B!=jbodies.end(); ++B ) {
      B->SRC[0] = B->SRC[1];
    }
    downward(cells,jcells,1);
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      B->TRG[2] = B->TRG[0];
      B->TRG[0] = 0;
    }

    jcells.resize(numCells);
    for( B_iter B=jbodies.begin(); B!=jbodies.end(); ++B ) {
      B->SRC[0] = B->SRC[2];
    }
    downward(cells,jcells,1);
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      B->TRG[3] = B->TRG[0];
      B->SRC[3] = dx;
    }

    rbf(bodies,cells,2);
    rbf(bodies,cells,1);
    rbf(bodies,cells,0);
  }
};

#endif
