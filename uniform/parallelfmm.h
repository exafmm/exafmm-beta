#include "serialfmm.h"

namespace exafmm {
  class ParallelFMM : public SerialFMM {
  private:
    int EXTERNAL;
    std::ofstream fid;
    MPI_Request *requests;

    template<typename T>
    void print(T data) {
      for (int irank=0; irank<MPISIZE; irank++ ) {                // Loop over ranks
	MPI_Barrier(MPI_COMM_WORLD);                              //  Sync processes
	usleep(100);                                              //  Wait 100 milliseconds
	if (MPIRANK == irank) std::cout << data << " ";           //  If it's my turn print "data"
      }                                                           // End loop over ranks
      MPI_Barrier(MPI_COMM_WORLD);                                // Sync processes
      usleep(100);                                                // Wait 100 milliseconds
      if (MPIRANK == MPISIZE-1) std::cout << std::endl;           // New line
    }

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
#ifdef EXAFMM_IJHPCA
      char fname[256];
      sprintf(fname,"time%5.5d.dat",MPIRANK);
      fid.open(fname);
#endif
    }
    ~ParallelFMM() {
      delete[] requests;
      if(!EXTERNAL) MPI_Finalize();
#ifdef EXAFMM_IJHPCA
      fid.close();
#endif
    }

    void partitionComm() {
      int ix[3];
      for( int i=0; i<MPISIZE; i++ ) sendBodiesCount[i] = 0;
      assert(numBodies % 3 == 0);
      for( int i=0; i<numBodies; i++ ) {
	setGlobIndex(i,ix);
	int sendRank = getGlobKey(ix,maxGlobLevel);
	Rank[i] = sendRank;
	sendBodiesCount[sendRank] += 4;
      }
      for( int i=0; i<MPISIZE; i++ ) assert(sendBodiesCount[i] % 12 == 0);
      MPI_Alltoall(sendBodiesCount,1,MPI_INT,recvBodiesCount,1,MPI_INT,MPI_COMM_WORLD);
      sendBodiesDispl[0] = recvBodiesDispl[0] = 0;
      for( int i=1; i<MPISIZE; i++ ) {
	sendBodiesDispl[i] = sendBodiesDispl[i-1] + sendBodiesCount[i-1];
	recvBodiesDispl[i] = recvBodiesDispl[i-1] + recvBodiesCount[i-1];
      }
      sort(Jbodies,sendJbodies,Index,sendIndex,Rank);
      MPI_Alltoallv(sendJbodies[0], sendBodiesCount, sendBodiesDispl, MPI_FLOAT,
		    recvJbodies[0], recvBodiesCount, recvBodiesDispl, MPI_FLOAT,
		    MPI_COMM_WORLD);
      int newBodies = (recvBodiesDispl[MPISIZE-1] + recvBodiesCount[MPISIZE-1]) / 4;
      for( int i=0; i<newBodies; i++ ) {
	for_4d Jbodies[i][d] = recvJbodies[i][d];
      }
      sort(Ibodies,sendJbodies,Index,sendIndex,Rank);
      MPI_Alltoallv(sendJbodies[0], sendBodiesCount, sendBodiesDispl, MPI_FLOAT,
		    recvJbodies[0], recvBodiesCount, recvBodiesDispl, MPI_FLOAT,
		    MPI_COMM_WORLD);
      for( int i=0; i<MPISIZE; i++ ) {
	sendBodiesCount[i] /= 4;
	sendBodiesDispl[i] /= 4;
	recvBodiesCount[i] /= 4;
	recvBodiesDispl[i] /= 4;
      }
      MPI_Alltoallv(sendIndex, sendBodiesCount, sendBodiesDispl, MPI_INT,
		    recvIndex, recvBodiesCount, recvBodiesDispl, MPI_INT,
		    MPI_COMM_WORLD);
      numBodies = newBodies;
      for( int i=0; i<numBodies; i++ ) {
	Index[i] = recvIndex[i];
	for_4d Ibodies[i][d] = recvJbodies[i][d];
      }
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
	      assert(0<=sendRank && sendRank<MPISIZE);
	      assert(0<=recvRank && recvRank<MPISIZE);
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
      MPI_Waitall(52,requests,stats);
    }

    void P2PRecv() {
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
#ifdef EXAFMM_IJHPCA
      MPI_Barrier(MPI_COMM_WORLD);
      logger::startTimer("M2L Comm");
#endif
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
      MPI_Waitall(26,requests,stats);
    }

    void M2LRecv(int lev) {
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
#ifdef EXAFMM_IJHPCA
      double time = logger::stopTimer("M2L Comm", 0);
      fid << time << std::endl;
      logger::resetTimer("M2L Comm");
#endif
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
      MPI_Status stats[8];
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
      MPI_Waitall(numComm,requests,stats);
    }

    void globM2MRecv(int level) {
      MPI_Status stats[8];
      int numChild[3];
      for_3d numChild[d] = numPartition[level][d] / numPartition[level-1][d];
      int ix[3];
      for_3d ix[d] = IX[level][d];
      int jxoff[3];
      for_3d jxoff[d] = (ix[d] / numChild[d]) * numChild[d];
      int iforward = 0;
      int numComm = numChild[0] * numChild[1] * numChild[2] - 1;
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
	globM2MSend(lev);
	globM2MRecv(lev);
	logger::stopTimer("Comm LET cells", 0);
	logger::startTimer("Upward pass");
	int numChild[3];
	for_3d numChild[d] = numPartition[lev][d] / numPartition[lev-1][d];
	int jxoff[3];
	for_3d jxoff[d] = (IX[lev][d] / numChild[d]) * numChild[d];
	int childOffset = globLevelOffset[lev];
	int parentOffset = globLevelOffset[lev-1];
	real_t diameter[3];
	for_3d diameter[d] = 2 * RGlob[d] / numPartition[lev][d];
	int jx[3];
	for( jx[2]=jxoff[2]; jx[2]<jxoff[2]+numChild[2]; jx[2]++ ) {
	  for( jx[1]=jxoff[1]; jx[1]<jxoff[1]+numChild[1]; jx[1]++ ) {
	    for( jx[0]=jxoff[0]; jx[0]<jxoff[0]+numChild[0]; jx[0]++ ) {
	      int ix[3];
	      for_3d ix[d] = jx[d] / numChild[d];
	      int c = getGlobKey(jx,lev) + childOffset;
	      int p = getGlobKey(ix,lev-1) + parentOffset;
	      real_t dX[3];
	      for_3d dX[d] = (ix[d] + .5) * numChild[d] * diameter[d] - (jx[d] + .5) * diameter[d];
	      real_t M[MTERM];
	      real_t C[LTERM];
	      C[0] = 1;
	      powerM(C,dX);
	      for_m M[m] = globMultipole[c][m];
	      for_m globMultipole[p][m] += C[m] * M[0];
	      M2MSum(globMultipole[p],C,M);
	    }
	  }
	}
	logger::stopTimer("Upward pass", 0);
      }
      logger::startTimer("Comm LET cells");
      gatherMultipoles();
      logger::stopTimer("Comm LET cells", 0);
      logger::startTimer("Upward pass");
      for( int lev=gatherLevel; lev>0; lev-- ) {
	int numChild[3];
	for_3d numChild[d] = numPartition[lev][d] / numPartition[lev-1][d];
	int childOffset = globLevelOffset[lev];
	int parentOffset = globLevelOffset[lev-1];
	real_t diameter[3];
	for_3d diameter[d] = 2 * RGlob[d] / numPartition[lev][d];
	int jx[3];
	for( jx[2]=0; jx[2]<numPartition[lev][2]; jx[2]++ ) {
	  for( jx[1]=0; jx[1]<numPartition[lev][1]; jx[1]++ ) {
	    for( jx[0]=0; jx[0]<numPartition[lev][0]; jx[0]++ ) {
	      int ix[3];
	      for_3d ix[d] = jx[d] / numChild[d];
	      int c = getGlobKey(jx,lev) + childOffset;
	      int p = getGlobKey(ix,lev-1) + parentOffset;
	      real_t dX[3];
	      for_3d dX[d] = (ix[d] + .5) * numChild[d] * diameter[d] - (jx[d] + .5) * diameter[d];
	      real_t M[MTERM];
	      real_t C[LTERM];
	      C[0] = 1;
	      powerM(C,dX);
	      for_m M[m] = globMultipole[c][m];
	      for_m globMultipole[p][m] += C[m] * M[0];
	      M2MSum(globMultipole[p],C,M);
	    }
	  }
	}
      }
      logger::stopTimer("Upward pass", 0);
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
      MPI_Waitall(26,requests,stats);
    }

    void globM2LRecv(int level) {
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

    void globM2L() {
      for( int lev=maxGlobLevel; lev>0; lev-- ) {
	MPI_Barrier(MPI_COMM_WORLD);
	logger::startTimer("Comm LET cells");
	if( lev > gatherLevel ) {
#ifdef EXAFMM_IJHPCA
	  logger::startTimer("M2L Comm");
	  globM2LSend(lev);
	  globM2LRecv(lev);
	  double time = logger::stopTimer("M2L Comm", 0);
	  fid << time << std::endl;
	  logger::resetTimer("M2L Comm");
#endif
	}
	logger::stopTimer("Comm LET cells");
	logger::startTimer("Traverse");
	int nxmin[3] = {0, 0, 0};
	int nxmax[3] = {numPartition[lev-1][0]-1,numPartition[lev-1][1]-1,numPartition[lev-1][2]-1};
	int nunit[3] = {numPartition[lev][0],numPartition[lev][1],numPartition[lev][2]};
	real_t diameter[3];
	for_3d diameter[d] = 2 * RGlob[d] / numPartition[lev][d];
	if( numImages != 0 ) {
	  for_3d nxmin[d] = -nxmax[d] - 1;
	  for_3d nxmax[d] = 2 * nxmax[d] + 1;
	}
	real_t L[LTERM];
	for_l L[l] = 0;
	int ix[3];
	for_3d ix[d] = IX[lev][d];
	int ixp[3];
	for_3d ixp[d] = IX[lev-1][d];
	int jxmin[3];
	for_3d jxmin[d] =  EXAFMM_MAX(nxmin[d], ixp[d] - 1)      * numPartition[lev][d] / numPartition[lev-1][d];
	int jxmax[3];
	for_3d jxmax[d] = (EXAFMM_MIN(nxmax[d], ixp[d] + 1) + 1) * numPartition[lev][d] / numPartition[lev-1][d];
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
		real_t M[MTERM];
		for_m M[m] = globMultipole[j][m];
		real_t dX[3];
		for_3d dX[d] = (ix[d] - jx[d]) * diameter[d];
		real_t invR2 = 1. / (dX[0] * dX[0] + dX[1] * dX[1] + dX[2] * dX[2]);
		real_t invR  = sqrt(invR2);
		real_t C[LTERM];
		getCoef(C,dX,invR2,invR);
		M2LSum(L,C,M);
	      }
	    }
	  }
	}
	for_l globLocal[lev][l] += L[l];
	logger::stopTimer("Traverse", 0);
      }
    }

    void globL2L() {
      for( int lev=1; lev<=maxGlobLevel; lev++ ) {
	real_t diameter[3];
	for_3d diameter[d] = 2 * RGlob[d] / numPartition[lev][d];
	real_t dX[3];
	for_3d dX[d] = (IX[lev][d] + .5) * diameter[d] - (IX[lev-1][d] + .5) * 2 * diameter[d];
	real_t C[LTERM];
	C[0] = 1;
	powerL(C,dX);
	for_l globLocal[lev][l] += globLocal[lev-1][l];
	for( int l=1; l<LTERM; l++ ) globLocal[lev][0] += C[l] * globLocal[lev-1][l];
	L2LSum(globLocal[lev],C,globLocal[lev-1]);
      }
      for_l Local[0][l] += globLocal[maxGlobLevel][l];
    }
  };
}
