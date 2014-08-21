#include "serialfmm.h"

class ParallelFMM : public SerialFMM {
private:
  int EXTERNAL;
  MPI_Request *requests;
#if PRINT_COMM
  std::ofstream fid;
#endif

  void gatherMultipoles() {
    int i = getGlobKey(IX[gatherLevel],gatherLevel) + globLevelOffset[gatherLevel];
    for_m sendMultipole[0][m] = globMultipole[i][m];
    int numGather = numPartition[gatherLevel][0] * numPartition[gatherLevel][1] * numPartition[gatherLevel][2];
    assert( numGather <= numSendCells ); // resize recvMultipole to avoid this
    int rank;
    MPI_Comm_rank(MPI_COMM_LOCAL,&rank);
    if( rank == 0 ) {
      MPI_Allgather(sendMultipole[0],MTERM,MPI_FLOAT,
                    recvMultipole[0],MTERM,MPI_FLOAT,MPI_COMM_GLOBAL);
    }
    MPI_Bcast(recvMultipole[0],numGather*MTERM,MPI_FLOAT,0,MPI_COMM_LOCAL);
    for( int c=0; c<numGather; c++ ) {
      for_m globMultipole[c+globLevelOffset[gatherLevel]][m] = recvMultipole[c][m];
    }
  }

public:
  ParallelFMM() {
    int argc(0);
    char **argv;
    MPI_Initialized(&EXTERNAL);
    if(!EXTERNAL) MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&MPISIZE);
    MPI_Comm_rank(MPI_COMM_WORLD,&MPIRANK);
    printNow = MPIRANK == 0;
    requests = new MPI_Request [104];
#if PRINT_COMM
    char fname[256];
    sprintf(fname,"time%4.4d.dat",MPIRANK);
    fid.open(fname);
#endif
  }
  ~ParallelFMM() {
#if PRINT_COMM
    fid.close();
#endif
    delete[] requests;
    if(!EXTERNAL) MPI_Finalize();
  }

  void P2PSend() {
    MPI_Status stats[52];
    int rankOffset = 13 * numLeafs;
    int ixc[3];
    getGlobIndex(ixc,MPIRANK,maxGlobLevel);
    int nunit[3];
    for_3d nunit[d] = numPartition[maxGlobLevel][d];
    int ileaf = 0;
    int iforward = 0;
    int ix[3];
    float commBytes = 0;
    for( ix[2]=-1; ix[2]<=1; ix[2]++ ) {
      for( ix[1]=-1; ix[1]<=1; ix[1]++ ) {
        for( ix[0]=-1; ix[0]<=1; ix[0]++ ) {
          if( ix[0] != 0 || ix[1] != 0 || ix[2] != 0 ) {
            assert( ileaf == leafsDispl[iforward] );
            int ibody = bodiesDispl[iforward];
            int nxmin[3] = {(1 << maxLevel) - 1, 0, 0};
            int nxmax[3] = {1 << maxLevel, 1 << maxLevel, 1};
            int jx[3];
            for( jx[2]=nxmin[ix[2]+1]; jx[2]<nxmax[ix[2]+1]; jx[2]++ ) {
              for( jx[1]=nxmin[ix[1]+1]; jx[1]<nxmax[ix[1]+1]; jx[1]++ ) {
                for( jx[0]=nxmin[ix[0]+1]; jx[0]<nxmax[ix[0]+1]; jx[0]++, ileaf++ ) {
                  int jxp[3] = {jx[0], jx[1], jx[2]};
                  int j = getKey(jxp,maxLevel,false) + rankOffset;
                  sendLeafs[ileaf][0] = ibody;
                  for( int jbody=Leafs[j][0]; jbody<Leafs[j][1]; ibody++, jbody++ ) {
                    for_4d sendJbodies[ibody][d] = Jbodies[jbody][d];
                  }
                  sendLeafs[ileaf][1] = ibody;
                }
              }
            }
            if(iforward != 25 ) {
               if( ibody > bodiesDispl[iforward+1] ) std::cout << "ibody: " << ibody << " bodiesDispl: " << bodiesDispl[iforward+1] << " @rank: " << MPIRANK << std::endl;
            }
            int ixp[3];
            for_3d ixp[d] = (ixc[d] - ix[d] + nunit[d]) % nunit[d];
            int sendRank = getGlobKey(ixp,maxGlobLevel);
            for_3d ixp[d] = (ixc[d] + ix[d] + nunit[d]) % nunit[d];
            int recvRank = getGlobKey(ixp,maxGlobLevel);
            int sendDispl = leafsDispl[iforward];
            int sendCount = leafsCount[iforward];
            commBytes += sendCount * 2 * 4;
            MPI_Isend(sendLeafs[sendDispl],sendCount*2,MPI_INT,
                      sendRank,iforward,MPI_COMM_WORLD,&requests[iforward]);
            int recvDispl = leafsDispl[iforward];
            int recvCount = leafsCount[iforward];
            MPI_Irecv(recvLeafs[recvDispl],recvCount*2,MPI_INT,
                      recvRank,iforward,MPI_COMM_WORLD,&requests[iforward+52]);
            sendDispl = bodiesDispl[iforward];
            sendCount = bodiesCount[iforward];
            commBytes += sendCount * 4 * 4;
            MPI_Isend(sendJbodies[sendDispl],sendCount*4,MPI_FLOAT,
                      sendRank,iforward+26,MPI_COMM_WORLD,&requests[iforward+26]);
            recvDispl = bodiesDispl[iforward];
            recvCount = bodiesCount[iforward];
            MPI_Irecv(recvJbodies[recvDispl],recvCount*4,MPI_FLOAT,
                      recvRank,iforward+26,MPI_COMM_WORLD,&requests[iforward+78]);
            iforward++;
          }
        }
      }
    }
#if PRINT_COMM
    int cells = (pow((1 << maxLevel) + 2,3) - (1 << (3 * maxLevel)));
    float theoBytes = cells * (2 + numBodies / numLeafs * 8) * 4;
    fid << "level : " << maxGlobLevel+maxLevel << " P2P comm (theoretical) : " << std::setw(8) << theoBytes << ",   (actual) : " << std::setw(8) << commBytes << " Bytes" << std::endl;
#endif
    MPI_Waitall(52,requests,stats);
  }

  void P2PRecv() const {
    MPI_Status stats[52];
    MPI_Waitall(52,&requests[52],stats);
    int ileaf = 0;
    int iforward = 0;
    int ix[3];
    for( ix[2]=-1; ix[2]<=1; ix[2]++ ) {
      for( ix[1]=-1; ix[1]<=1; ix[1]++ ) {
        for( ix[0]=-1; ix[0]<=1; ix[0]++ ) {
          if( ix[0] != 0 || ix[1] != 0 || ix[2] != 0 ) {
            assert( ileaf == leafsDispl[iforward] );
            int rankIndex = (ix[0] + 1) + 3 * (ix[1] + 1) + 9 * (ix[2] + 1);
            int rankOffset = rankIndex * numLeafs;
            int ibody = numBodies + bodiesDispl[iforward];
            int nxmin[3] = {(1 << maxLevel) - 1, 0, 0};
            int nxmax[3] = {1 << maxLevel, 1 << maxLevel, 1};
            int jx[3];
            for( jx[2]=nxmin[ix[2]+1]; jx[2]<nxmax[ix[2]+1]; jx[2]++ ) {
              for( jx[1]=nxmin[ix[1]+1]; jx[1]<nxmax[ix[1]+1]; jx[1]++ ) {
                for( jx[0]=nxmin[ix[0]+1]; jx[0]<nxmax[ix[0]+1]; jx[0]++, ileaf++ ) {
                  int jxp[3] = {jx[0], jx[1], jx[2]};
                  int j = getKey(jxp,maxLevel,false) + rankOffset;
                  Leafs[j][0] = ibody;
                  for( int jbody=recvLeafs[ileaf][0]; jbody<recvLeafs[ileaf][1]; ibody++, jbody++ ) {
                    for_4d Jbodies[ibody][d] = recvJbodies[jbody][d];
                  }
                  Leafs[j][1] = ibody;
                }
              }
            }
            iforward++;
          }
        }
      }
    }
  }

  void M2LSend(int lev) {
    MPI_Status stats[26];
    int rankOffset = 13 * numCells;
    int ixc[3];
    getGlobIndex(ixc,MPIRANK,maxGlobLevel);
    int nunit[3];
    for_3d nunit[d] = numPartition[maxGlobLevel][d];
    int nxmin[3] = {(1 << lev) - 2, 0, 0};
    int nxmax[3] = {1 << lev, 1 << lev, 2};
    int i = 0;
    int iforward = 0;
    int ix[3];
    float commBytes = 0;
    for( ix[2]=-1; ix[2]<=1; ix[2]++ ) {
      for( ix[1]=-1; ix[1]<=1; ix[1]++ ) {
        for( ix[0]=-1; ix[0]<=1; ix[0]++ ) {
          if( ix[0] != 0 || ix[1] != 0 || ix[2] != 0 ) {
            int jx[3];
            for( jx[2]=nxmin[ix[2]+1]; jx[2]<nxmax[ix[2]+1]; jx[2]++ ) {
              for( jx[1]=nxmin[ix[1]+1]; jx[1]<nxmax[ix[1]+1]; jx[1]++ ) {
                for( jx[0]=nxmin[ix[0]+1]; jx[0]<nxmax[ix[0]+1]; jx[0]++, i++ ) {
                  int jxp[3] = {jx[0], jx[1], jx[2]};
                  int j = getKey(jxp,lev) + rankOffset;
                  for_m sendMultipole[i][m] = Multipole[j][m];
                  commBytes += MTERM * 4;
                }
              }
            }
            int ixp[3];
            for_3d ixp[d] = (ixc[d] - ix[d] + nunit[d]) % nunit[d];
            int sendRank = getGlobKey(ixp,maxGlobLevel);
            for_3d ixp[d] = (ixc[d] + ix[d] + nunit[d]) % nunit[d];
            int recvRank = getGlobKey(ixp,maxGlobLevel);
            int sendDispl = multipoleDispl[lev][iforward];
            int sendCount = multipoleCount[lev][iforward];
            MPI_Isend(sendMultipole[sendDispl],sendCount*MTERM,MPI_FLOAT,
                      sendRank,iforward,MPI_COMM_WORLD,&requests[iforward]);
            int recvDispl = multipoleDispl[lev][iforward];
            int recvCount = multipoleCount[lev][iforward];
            MPI_Irecv(recvMultipole[recvDispl],recvCount*MTERM,MPI_FLOAT,
                      recvRank,iforward,MPI_COMM_WORLD,&requests[iforward+26]);
            iforward++;
          }
        }
      }
    }
#if PRINT_COMM
    int cells = (pow((1 << lev) + 4,3) - (1 << (3 * lev)));
    float theoBytes = cells * MTERM * 4;
    fid << "level : " << maxGlobLevel+lev << " M2L comm (theoretical) : " << std::setw(8) << theoBytes << ",   (actual) : " << std::setw(8) << commBytes << " Bytes" << std::endl;
#endif
    MPI_Waitall(26,requests,stats);
  }

  void M2LRecv(int lev) const {
    MPI_Status stats[26];
    int nxmin[3] = {(1 << lev) - 2, 0, 0};
    int nxmax[3] = {1 << lev, 1 << lev, 2};
    for( int iforward=0; iforward<26; iforward++ ) {
      int irequest;
      MPI_Waitany(26,&requests[26],&irequest,stats);
      int rankIndex = irequest < 13 ? irequest : irequest+1;
      int ix[3] = {rankIndex % 3, rankIndex / 3 % 3, rankIndex / 9};
      for_3d ix[d]--;
      int i = multipoleDispl[lev][irequest];
      int rankOffset = rankIndex * numCells;
      int jx[3];
      for( jx[2]=nxmin[ix[2]+1]; jx[2]<nxmax[ix[2]+1]; jx[2]++ ) {
        for( jx[1]=nxmin[ix[1]+1]; jx[1]<nxmax[ix[1]+1]; jx[1]++ ) {
          for( jx[0]=nxmin[ix[0]+1]; jx[0]<nxmax[ix[0]+1]; jx[0]++, i++ ) {
            int jxp[3] = {jx[0], jx[1], jx[2]};
            int j = getKey(jxp,lev) + rankOffset;
            for_m Multipole[j][m] = recvMultipole[i][m];
          }
        }
      }
    }
  }

  void rootGather() {
#pragma omp parallel for
    for(int i=0;i<numGlobCells;i++){
      for_m globMultipole[i][m] = 0;
    }
#pragma omp parallel for
    for( int lev=0; lev<=maxGlobLevel; lev++ ) {
      for_l globLocal[lev][l] = 0;
    }
  }

  void globM2MSend(int level) {
    int numChild[3];
    for_3d numChild[d] = numPartition[level][d] / numPartition[level-1][d];
    int numStride[3];
    for_3d numStride[d] = numPartition[maxGlobLevel][d] / numPartition[level][d];
    int ix[3];
    for_3d ix[d] = IX[level][d];
    int ixoff[3];
    for_3d ixoff[d] = IX[maxGlobLevel][d] % numStride[d];
    int jxoff[3];
    for_3d jxoff[d] = (IX[level][d] / numChild[d]) * numChild[d];
    int i = getGlobKey(ix,level) + globLevelOffset[level];
    for_m sendMultipole[0][m] = globMultipole[i][m];
    int iforward = 0;
    int numComm = numChild[0] * numChild[1] * numChild[2] - 1;
    MPI_Status *stats = new MPI_Status[numComm];
    float commBytes = 0;
    int jx[3];
    for( jx[2]=jxoff[2]; jx[2]<jxoff[2]+numChild[2]; jx[2]++ ) {
      for( jx[1]=jxoff[1]; jx[1]<jxoff[1]+numChild[1]; jx[1]++ ) {
        for( jx[0]=jxoff[0]; jx[0]<jxoff[0]+numChild[0]; jx[0]++ ) {
          if( ix[0] != jx[0] || ix[1] != jx[1] || ix[2] != jx[2] ) {
            int jxp[3];
            for_3d jxp[d] = ixoff[d] + jx[d] * numStride[d];
            int commRank = getGlobKey(jxp,maxGlobLevel);
            commBytes += MTERM * 4;
            MPI_Isend(sendMultipole[0],MTERM,MPI_FLOAT,
                      commRank,0,MPI_COMM_WORLD,&requests[iforward]);
            MPI_Irecv(recvMultipole[iforward],MTERM,MPI_FLOAT,
                      commRank,0,MPI_COMM_WORLD,&requests[iforward+numComm]);
            iforward++;
          }
        }
      }
    }
#if PRINT_COMM
    float theoBytes = numComm * MTERM * 4;
    fid << "level : " << level << " M2M comm (theoretical) : " << std::setw(8) << theoBytes << ",   (actual) : " << std::setw(8) << commBytes << " Bytes" << std::endl;
#endif
    MPI_Waitall(numComm,requests,stats);
  }

  void globM2MRecv(int level) const {
    int numChild[3];
    for_3d numChild[d] = numPartition[level][d] / numPartition[level-1][d];
    int ix[3];
    for_3d ix[d] = IX[level][d];
    int jxoff[3];
    for_3d jxoff[d] = (ix[d] / numChild[d]) * numChild[d];
    int iforward = 0;
    int numComm = numChild[0] * numChild[1] * numChild[2] - 1;
    MPI_Status *stats = new MPI_Status[numComm];
    MPI_Waitall(numComm,&requests[numComm],stats);
    int jx[3];
    for( jx[2]=jxoff[2]; jx[2]<jxoff[2]+numChild[2]; jx[2]++ ) {
      for( jx[1]=jxoff[1]; jx[1]<jxoff[1]+numChild[1]; jx[1]++ ) {
        for( jx[0]=jxoff[0]; jx[0]<jxoff[0]+numChild[0]; jx[0]++ ) {
          if( ix[0] != jx[0] || ix[1] != jx[1] || ix[2] != jx[2] ) {
            int j = getGlobKey(jx,level) + globLevelOffset[level];
            for_m globMultipole[j][m] = recvMultipole[iforward][m];
            iforward++;
          }
        }
      }
    }
  }

  void globM2M() {
    int rankOffset = 13 * numCells;
    int i = MPIRANK + globLevelOffset[maxGlobLevel];
    for_m globMultipole[i][m] = Multipole[rankOffset][m];
    for( int lev=maxGlobLevel; lev>gatherLevel; lev-- ) {
      logger::startTimer("Comm LET cells");
      double tic = getTime();
      globM2MSend(lev);
      globM2MRecv(lev);
      double toc = getTime();
      if( printNow ) printf("M2M Comm: %lf @ lev: %d\n",toc-tic,lev);
      logger::stopTimer("Comm LET cells");
      logger::startTimer("Upward pass");
      tic = getTime();
      int numChild[3];
      for_3d numChild[d] = numPartition[lev][d] / numPartition[lev-1][d];
      int jxoff[3];
      for_3d jxoff[d] = (IX[lev][d] / numChild[d]) * numChild[d];
      int childOffset = globLevelOffset[lev];
      int parentOffset = globLevelOffset[lev-1];
      real diameter[3];
      for_3d diameter[d] = 2 * RGlob[d] / numPartition[lev][d];
      int jx[3];
      for( jx[2]=jxoff[2]; jx[2]<jxoff[2]+numChild[2]; jx[2]++ ) {
        for( jx[1]=jxoff[1]; jx[1]<jxoff[1]+numChild[1]; jx[1]++ ) {
          for( jx[0]=jxoff[0]; jx[0]<jxoff[0]+numChild[0]; jx[0]++ ) {
            int ix[3];
            for_3d ix[d] = jx[d] / numChild[d];
            int c = getGlobKey(jx,lev) + childOffset;
            int p = getGlobKey(ix,lev-1) + parentOffset;
            real dist[3];
            for_3d dist[d] = (ix[d] + .5) * numChild[d] * diameter[d] - (jx[d] + .5) * diameter[d];
            real M[MTERM];
            real C[LTERM];
            C[0] = 1;
	    powerM(C,dist);
            for_m M[m] = globMultipole[c][m];
            for_m globMultipole[p][m] += C[m] * M[0];
	    M2MSum(globMultipole[p],C,M);
          }
        }
      }
      toc = getTime();
      if( printNow ) printf("M2M Glob: %lf @ lev: %d\n",toc-tic,lev);
      logger::stopTimer("Upward pass");
    }
    logger::startTimer("Comm LET cells");
    double tic = getTime();
    gatherMultipoles();
    double toc = getTime();
    if( printNow ) printf("M2M Comm: %lf @ lev: %d\n",toc-tic,gatherLevel);
    logger::stopTimer("Comm LET cells");
    logger::startTimer("Upward pass");
    for( int lev=gatherLevel; lev>0; lev-- ) {
      tic = getTime();
      int numChild[3];
      for_3d numChild[d] = numPartition[lev][d] / numPartition[lev-1][d];
      int childOffset = globLevelOffset[lev];
      int parentOffset = globLevelOffset[lev-1];
      real diameter[3];
      for_3d diameter[d] = 2 * RGlob[d] / numPartition[lev][d];
      int jx[3];
      for( jx[2]=0; jx[2]<numPartition[lev][2]; jx[2]++ ) {
        for( jx[1]=0; jx[1]<numPartition[lev][1]; jx[1]++ ) {
          for( jx[0]=0; jx[0]<numPartition[lev][0]; jx[0]++ ) {
            int ix[3];
            for_3d ix[d] = jx[d] / numChild[d];
            int c = getGlobKey(jx,lev) + childOffset;
            int p = getGlobKey(ix,lev-1) + parentOffset;
            real dist[3];
            for_3d dist[d] = (ix[d] + .5) * numChild[d] * diameter[d] - (jx[d] + .5) * diameter[d];
            real M[MTERM];
            real C[LTERM];
            C[0] = 1;
            powerM(C,dist);
            for_m M[m] = globMultipole[c][m];
            for_m globMultipole[p][m] += C[m] * M[0];
            M2MSum(globMultipole[p],C,M);
          }
        }
      }
      toc = getTime();
      if( printNow ) printf("M2M Glob: %lf @ lev: %d\n",toc-tic,lev);
    }
    logger::stopTimer("Upward pass");
  }

  void globM2LSend(int level) {
    MPI_Status stats[26];
    int numChild[3];
    for_3d numChild[d] = numPartition[level][d] / numPartition[level-1][d];
    int numStride[3];
    for_3d numStride[d] = numPartition[maxGlobLevel][d] / numPartition[level-1][d];
    int ixc[3];
    for_3d ixc[d] = IX[level-1][d];
    int ixoff[3];
    for_3d ixoff[d] = IX[maxGlobLevel][d] % numStride[d];
    int numGroup = numChild[0] * numChild[1] * numChild[2];
    float commBytes = 0;
    int i = 0;
    int iforward = 0;
    int ix[3];
    for( ix[2]=-1; ix[2]<=1; ix[2]++ ) {
      for( ix[1]=-1; ix[1]<=1; ix[1]++ ) {
        for( ix[0]=-1; ix[0]<=1; ix[0]++ ) {
          if( ix[0] != 0 || ix[1] != 0 || ix[2] != 0 ) {
            int jx[3];
            for( jx[2]=ixc[2]*numChild[2]; jx[2]<(ixc[2]+1)*numChild[2]; jx[2]++ ) {
              for( jx[1]=ixc[1]*numChild[1]; jx[1]<(ixc[1]+1)*numChild[1]; jx[1]++ ) {
                for( jx[0]=ixc[0]*numChild[0]; jx[0]<(ixc[0]+1)*numChild[0]; jx[0]++, i++ ) {
                  int j = getGlobKey(jx,level) + globLevelOffset[level];
                  for_m sendMultipole[i][m] = globMultipole[j][m];
                }
              }
            }
            int ixp[3];
            for_3d ixp[d] = (ixc[d] + ix[d] + numPartition[level-1][d]) % numPartition[level-1][d];
            for_3d ixp[d] = ixoff[d] + ixp[d] * numStride[d];
            int sendRank = getGlobKey(ixp,maxGlobLevel);
            commBytes += numGroup * MTERM * 4;
            MPI_Isend(sendMultipole[iforward*numGroup],numGroup*MTERM,MPI_FLOAT,
                      sendRank,iforward,MPI_COMM_WORLD,&requests[iforward]);
            for_3d ixp[d] = (ixc[d] - ix[d] + numPartition[level-1][d]) % numPartition[level-1][d];
            for_3d ixp[d] = ixoff[d] + ixp[d] * numStride[d];
            int recvRank = getGlobKey(ixp,maxGlobLevel);
            MPI_Irecv(recvMultipole[iforward*numGroup],numGroup*MTERM,MPI_FLOAT,
                      recvRank,iforward,MPI_COMM_WORLD,&requests[iforward+26]);
            iforward++;
          }
        }
      }
    }
#if PRINT_COMM
    float theoBytes = 26 * numGroup * MTERM * 4;
    fid << "level : " << level << " M2L comm (theoretical) : " << std::setw(8) << theoBytes << ",   (actual) : " << std::setw(8) << commBytes << " Bytes" << std::endl;
#endif
    MPI_Waitall(26,requests,stats);
  }

  void globM2LRecv(int level) const {
    MPI_Status stats[26];
    MPI_Waitall(26,&requests[26],stats);
    int numChild[3];
    for_3d numChild[d] = numPartition[level][d] / numPartition[level-1][d];
    int ixc[3];
    for_3d ixc[d] = IX[level-1][d];
    int i = 0;
    int iforward = 0;
    int ix[3];
    for( ix[2]=-1; ix[2]<=1; ix[2]++ ) {
      for( ix[1]=-1; ix[1]<=1; ix[1]++ ) {
        for( ix[0]=-1; ix[0]<=1; ix[0]++ ) {
          if( ix[0] != 0 || ix[1] != 0 || ix[2] != 0 ) {
            int ixp[3];
            for_3d ixp[d] = (ixc[d] - ix[d] + numPartition[level-1][d]) % numPartition[level-1][d];
            int jx[3];
            for( jx[2]=ixp[2]*numChild[2]; jx[2]<(ixp[2]+1)*numChild[2]; jx[2]++ ) {
              for( jx[1]=ixp[1]*numChild[1]; jx[1]<(ixp[1]+1)*numChild[1]; jx[1]++ ) {
                for( jx[0]=ixp[0]*numChild[0]; jx[0]<(ixp[0]+1)*numChild[0]; jx[0]++, i++ ) {
                  int j = getGlobKey(jx,level) + globLevelOffset[level];
                  for_m globMultipole[j][m] = recvMultipole[i][m];
                }
              }
            }
            iforward++;
          }
        }
      }
    }
  }

  void globM2L(std::ofstream &fid2) {
    for( int lev=maxGlobLevel; lev>0; lev-- ) {
      MPI_Barrier(MPI_COMM_WORLD);
      logger::startTimer("Comm LET cells");
      double tic = getTime();
      if( lev > gatherLevel ) {
        globM2LSend(lev);
        globM2LRecv(lev);
      }
      double toc = getTime();
      if( lev > 1 ) fid2 << toc-tic << std::endl;
      if( printNow ) printf("M2L Comm: %lf @ lev: %d\n",toc-tic,lev);
      logger::stopTimer("Comm LET cells");
      logger::startTimer("Traverse");
      tic = getTime();
      int nxmin[3] = {0, 0, 0};
      int nxmax[3] = {numPartition[lev-1][0]-1,numPartition[lev-1][1]-1,numPartition[lev-1][2]-1};
      int nunit[3] = {numPartition[lev][0],numPartition[lev][1],numPartition[lev][2]};
      real diameter[3];
      for_3d diameter[d] = 2 * RGlob[d] / numPartition[lev][d];
      if( numImages != 0 ) {
        for_3d nxmin[d] = -nxmax[d] - 1;
        for_3d nxmax[d] = 2 * nxmax[d] + 1;
      }
      real L[LTERM];
      for_l L[l] = 0;
      int ix[3];
      for_3d ix[d] = IX[lev][d];
      int ixp[3];
      for_3d ixp[d] = IX[lev-1][d];
      int jxmin[3];
      for_3d jxmin[d] =  FMMMAX(nxmin[d], ixp[d] - 1)      * numPartition[lev][d] / numPartition[lev-1][d];
      int jxmax[3];
      for_3d jxmax[d] = (FMMMIN(nxmax[d], ixp[d] + 1) + 1) * numPartition[lev][d] / numPartition[lev-1][d];
      int jx[3];
      for( jx[2]=jxmin[2]; jx[2]<jxmax[2]; jx[2]++ ) {
        for( jx[1]=jxmin[1]; jx[1]<jxmax[1]; jx[1]++ ) {
          for( jx[0]=jxmin[0]; jx[0]<jxmax[0]; jx[0]++ ) {
            if(jx[0] < ix[0]-1 || ix[0]+1 < jx[0] ||
               jx[1] < ix[1]-1 || ix[1]+1 < jx[1] ||
               jx[2] < ix[2]-1 || ix[2]+1 < jx[2]) {
              int jxp[3];
              for_3d jxp[d] = (jx[d] + nunit[d]) % nunit[d];
              int j = getGlobKey(jxp,lev) + globLevelOffset[lev];
              real M[MTERM];
              for_m M[m] = globMultipole[j][m];
              real dist[3];
              for_3d dist[d] = (ix[d] - jx[d]) * diameter[d];
              real invR2 = 1. / (dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2]);
              real invR  = sqrt(invR2);
              real C[LTERM];
	      getCoef(C,dist,invR2,invR);
	      M2LSum(L,C,M);
            }
          }
        }
      }
      for_l globLocal[lev][l] += L[l];
      toc = getTime();
      if( printNow ) printf("M2L Glob: %lf @ lev: %d\n",toc-tic,lev);
      logger::stopTimer("Traverse");
    }
  }

  void globL2L() const {
    for( int lev=1; lev<=maxGlobLevel; lev++ ) {
      real diameter[3];
      for_3d diameter[d] = 2 * RGlob[d] / numPartition[lev][d];
      real dist[3];
      for_3d dist[d] = (IX[lev][d] + .5) * diameter[d] - (IX[lev-1][d] + .5) * 2 * diameter[d];
      real C[LTERM];
      C[0] = 1;
      powerL(C,dist);
      for_l globLocal[lev][l] += globLocal[lev-1][l];
      for( int l=1; l<LTERM; l++ ) globLocal[lev][0] += C[l] * globLocal[lev-1][l];
      L2LSum(globLocal[lev],C,globLocal[lev-1]);
    }
    for_l Local[0][l] += globLocal[maxGlobLevel][l];
  }

  void globDirect() {
    const int numTarget = 100;
    MPI_Status stats[2];
    real (*Ibodies2)[4] = new real [numTarget][4];
    real (*Jbodies2)[4] = new real [numBodies][4];
    for( int i=0; i<numTarget; i++ ) {
      for_4d Ibodies2[i][d] = 0;
    }
    for( int i=0; i<numBodies; i++ ) {
      for_4d Jbodies2[i][d] = Jbodies[i][d];
    }
    const int sendRank = (MPIRANK + 1          ) % MPISIZE;
    const int recvRank = (MPIRANK - 1 + MPISIZE) % MPISIZE;
    for( int irank=0; irank<MPISIZE; irank++ ) {
      for( int i=0; i<numBodies; i++ ) {
        for_4d sendJbodies[i][d] = Jbodies2[i][d];
      }
      int newBodies = 0;
      MPI_Isend(&numBodies,1,MPI_INT,sendRank,
                1,MPI_COMM_WORLD,&requests[0]);
      MPI_Irecv(&newBodies,1,MPI_INT,recvRank,
                1,MPI_COMM_WORLD,&requests[1]);
      MPI_Waitall(2,requests,stats);

      MPI_Isend(sendJbodies[0],numBodies*4,MPI_FLOAT,sendRank,
                1,MPI_COMM_WORLD,&requests[0]);
      MPI_Irecv(recvJbodies[0],newBodies*4,MPI_FLOAT,recvRank,
                1,MPI_COMM_WORLD,&requests[1]);
      int prange = numImages == 0 ? 0 : pow(3,numImages - 1);
#pragma omp parallel for
      for( int i=0; i<numTarget; i++ ) {
        real bodies[4] = {0, 0, 0, 0};
        int jx[3];
        for( jx[2]=-prange; jx[2]<=prange; jx[2]++ ) {
          for( jx[1]=-prange; jx[1]<=prange; jx[1]++ ) {
            for( jx[0]=-prange; jx[0]<=prange; jx[0]++ ) {
              for( int j=0; j<numBodies; j++ ) {
                real dist[3];
                for_3d dist[d] = Jbodies[i][d] - Jbodies2[j][d] - jx[d] * 2 * RGlob[d];
                real R2 = dist[0] * dist[0] + dist[1] * dist[1] + dist[2] * dist[2];
                real invR2 = 1.0 / R2;
                if( R2 == 0 ) invR2 = 0;
                real invR = Jbodies2[j][3] * sqrt(invR2);
                real invR3 = invR2 * invR;
                bodies[0] += invR;
                for_3d bodies[d+1] -= dist[d] * invR3;
              }
            }
          }
        }
        for_4d Ibodies2[i][d] += bodies[d];
      }
      MPI_Waitall(2,requests,stats);
      numBodies = newBodies;
      for( int i=0; i<numBodies; i++ ) {
        for_4d Jbodies2[i][d] = recvJbodies[i][d];
      }
    }
    real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
    for( int i=0; i<numTarget; i++ ) {
      diff1 += (Ibodies[i][0] - Ibodies2[i][0]) * (Ibodies[i][0] - Ibodies2[i][0]);
      norm1 += Ibodies2[i][0] * Ibodies2[i][0];
      for_3d diff2 += (Ibodies[i][d+1] - Ibodies2[i][d+1]) * (Ibodies[i][d+1] - Ibodies2[i][d+1]);
      for_3d norm2 += Ibodies2[i][d+1] * Ibodies2[i][d+1];
    }
    real diff3 = 0, norm3 = 0, diff4 = 0, norm4 = 0;
    MPI_Reduce(&diff1, &diff3, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&norm1, &norm3, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&diff2, &diff4, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&norm2, &norm4, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if( MPIRANK == 0 ) printf("Err Pot : %lf\n",sqrt(diff3/norm3));
    if( MPIRANK == 0 ) printf("Err Forc: %lf\n",sqrt(diff4/norm4));
    delete[] Ibodies2;
    delete[] Jbodies2;
  }
};
