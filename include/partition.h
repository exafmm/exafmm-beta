#ifndef partition_h
#define partition_h
#include "mympi.h"
#include "types.h"
#include "sort.h"

class Partition : public MyMPI, virtual public Sort {
private:
  Bodies &bodies;                                               // Bodies to be partitioned
  int numBodies;                                                // Initial local number of bodies
protected:
  int LEVEL;                                                    // Level of the MPI process binary tree
  std::vector<vect> XMIN;                                       // Minimum position vector of bodies
  std::vector<vect> XMAX;                                       // Maximum position vector of bodies
  int oldnprocs;                                                // Old number of processes in the group
  int nprocs[2];                                                // Number of processes in the two split groups
  int offset[2];                                                // Offset of body in the two split groups
  int  color[3];                                                // Color of Gather, Scatter, and Alltoall communicators
  int    key[3];                                                // Key of Gather, Scatter, and Alltoall communicators
  MPI_Comm MPI_COMM[3];                                         // Communicators for Gather, Scatter, and Alltoall
public:
  Partition(Bodies &b) : bodies(b) {                            // Constructor
    numBodies = bodies.size();                                  // Initialize local number of bodies
    LEVEL = int(log(SIZE) / M_LN2 - 1e-5) + 1;                  // Level of the process binary tree
    XMIN.resize(LEVEL+1);                                       // Minimum position vector at each level
    XMAX.resize(LEVEL+1);                                       // Maximum position vector at each level
  }
  ~Partition() {}                                               // Destructor

  void setGlobDomain(real &R0, vect &X0) {                      // Set bounds of domain to be partitioned
    B_iter B = bodies.begin();                                  // Reset body iterator
    XMIN[0] = XMAX[0] = B->pos;                                 // Initialize xmin,xmax
    int const MPI_TYPE = getType(XMIN[0][0]);                   // Get MPI data type
    for( B=bodies.begin(); B!=bodies.end(); ++B ) {             // Loop over all bodies
      for( int d=0; d!=3; ++d ) {                               //  Loop over each dimension
        if     (B->pos[d] < XMIN[0][d]) XMIN[0][d] = B->pos[d]; //   Determine xmin
        else if(B->pos[d] > XMAX[0][d]) XMAX[0][d] = B->pos[d]; //   Determine xmax
      }                                                         //  End loop over each dimension
    }                                                           // End loop over all bodies
    vect buffer;                                                // Define recv buffer
    MPI_Reduce(XMAX[0],buffer,3,MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD);// Reduce global maximum
    XMAX[0] = buffer;                                           // Get data from buffer
    MPI_Bcast(XMAX[0],3,MPI_TYPE,0,MPI_COMM_WORLD);             // Broadcast global maximum
    MPI_Reduce(XMIN[0],buffer,3,MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD);// Reduce global minimum
    XMIN[0] = buffer;                                           // Get data from buffer
    MPI_Bcast(XMIN[0],3,MPI_TYPE,0,MPI_COMM_WORLD);             // Broadcast global minimum
    R0 = 0;                                                     // Initialize root radius
    for( int d=0; d!=3; ++d ) {                                 // Loop over each dimension
      X0[d] = (XMAX[0][d] + XMIN[0][d]) / 2;                    //  Calculate center of domain
      X0[d] = int(X0[d]+.5);                                    //  Shift center to nearest integer
      R0 = std::max(XMAX[0][d] - X0[d], R0);                    //  Calculate max distance from center
      R0 = std::max(X0[d] - XMIN[0][d], R0);                    //  Calculate max distance from center
    }                                                           // End loop over each dimension
    R0 = pow(2.,int(1. + log(R0) / M_LN2));                     // Add some leeway to root radius
  }

  void splitDomain(int iSplit, int l, int d) {
    real X = iSplit * (XMAX[l][d] - XMIN[l][d]) / numBodies + XMIN[l][d];
    XMAX[l+1] = XMAX[l];
    XMIN[l+1] = XMIN[l];
    if( color[0] % 2 == 0 )
      XMAX[l+1][d] = X;
    else
      XMIN[l+1][d] = X;
  }

  void binBodies(Bigints &Ibody, int d) {                       // Turn positions into indices of bins
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over all bodies
      Ibody[B-bodies.begin()] = bigint((B->pos[d] - XMIN[0][d]) // Bin body positions into integers
        / (XMAX[0][d] - XMIN[0][d]) * numBodies);
    }                                                           // End loop over all bodies
  }

  void shiftBodies(Bigints &Ibody, Bodies &buffer) {            // Send bodies to next rank (wrap around)
    int newSize;                                                // New number of bodies
    int oldSize = bodies.size();                                // Current number of bodies
    int const MPI_TYPE = getType(Ibody[0]);                     // Get MPI data type
    int const bytes = sizeof(bodies[0]);                        // Byte size of body structure
    int const isend = (RANK + 1       ) % SIZE;                 // Send to next rank (wrap around)
    int const irecv = (RANK - 1 + SIZE) % SIZE;                 // Receive from previous rank (wrap around)
    MPI_Request sreq,rreq;                                      // Send, receive request handles

    MPI_Isend(&oldSize,1,MPI_INT,irecv,0,MPI_COMM_WORLD,&sreq); // Send current number of bodies
    MPI_Irecv(&newSize,1,MPI_INT,isend,0,MPI_COMM_WORLD,&rreq); // Receive new number of bodies
    MPI_Wait(&sreq,MPI_STATUS_IGNORE);                          // Wait for send to complete
    MPI_Wait(&rreq,MPI_STATUS_IGNORE);                          // Wait for receive to complete

    ibuffer.resize(newSize);                                    // Resize buffer to new number of index
    MPI_Isend(&Ibody[0],  oldSize,MPI_TYPE,irecv,               // Send index to next rank
              1,MPI_COMM_WORLD,&sreq);
    MPI_Irecv(&ibuffer[0],newSize,MPI_TYPE,isend,               // Receive index from previous rank
              1,MPI_COMM_WORLD,&rreq);
    MPI_Wait(&sreq,MPI_STATUS_IGNORE);                          // Wait for send to complete
    MPI_Wait(&rreq,MPI_STATUS_IGNORE);                          // Wait for receive to complete
    Ibody = ibuffer;                                            // Copy index from buffer

    buffer.resize(newSize);                                     // Resize buffer to new number of bodies
    MPI_Isend(&bodies[0].pos[0],oldSize*bytes,MPI_BYTE,irecv,   // Send bodies to next rank
              2,MPI_COMM_WORLD,&sreq);
    MPI_Irecv(&buffer[0].pos[0],newSize*bytes,MPI_BYTE,isend,   // Receive bodies from previous rank
              2,MPI_COMM_WORLD,&rreq);
    MPI_Wait(&sreq,MPI_STATUS_IGNORE);                          // Wait for send to complete
    MPI_Wait(&rreq,MPI_STATUS_IGNORE);                          // Wait for receive to complete
    bodies = buffer;                                            // Copy bodies from buffer
  }

  int splitBodies(Bigints &Ibody, bigint iSplit) {
    int nth = 0;
    while( Ibody[nth] < iSplit && nth < int(Ibody.size() - 1) ) nth++;
    return nth;
  }

  void multisectionGetComm(int l) {
    int  oldcolor = color[0];
    oldnprocs = nprocs[0];
    for( int i=1; i>=0; --i ) {
      int pSplit = (oldnprocs + i) / 2;
      if( RANK - offset[0] < pSplit ) {
        nprocs[i] = pSplit;
        offset[i] = offset[0];
         color[i] =  color[0] * 2;
      } else {
        nprocs[i] = oldnprocs - pSplit;
        offset[i] = offset[0] + pSplit;
         color[i] =  color[0] * 2 + 1;
      }
      key[i] = RANK - offset[i];
    }
    key[2] = color[1] % 2;
    color[2] = key[1] + oldcolor * (1 << (LEVEL - l - 1));
    for( int i=0; i!=3; ++i )
      MPI_Comm_split(MPI_COMM_WORLD,color[i],key[i],&MPI_COMM[i]);
#ifdef DEBUG
    print("level : ",0);
    print(l,0);
    print("\n",0);
    print("key   : \n",0);
    print(key,0,3);
    print("color : \n",0);
    print(color,0,3);
#endif
  }

  void multisectionAlltoall(Bigints &Ibody, Bodies &buffer,
                            int nthLocal, int numLocal, int &newSize) {
    int const MPI_TYPE = getType(Ibody[0]);
    int const bytes = sizeof(bodies[0]);
    int scnt[2] = {nthLocal, numLocal - nthLocal};
    int rcnt[2] = {0, 0};
    MPI_Alltoall(scnt,1,MPI_INT,rcnt,1,MPI_INT,MPI_COMM[2]);
    int sdsp[2] = {0, scnt[0]};
    int rdsp[2] = {0, rcnt[0]};
    newSize = rcnt[0] + rcnt[1];
    if( color[0] != color[1] ) newSize = numLocal;
    ibuffer.resize(newSize);
    buffer.resize(newSize);

    MPI_Alltoallv(&Ibody[0],  scnt,sdsp,MPI_TYPE,
                  &ibuffer[0],rcnt,rdsp,MPI_TYPE,
                  MPI_COMM[2]);
    if( color[0] == color[1] && nthLocal != 0 ) Ibody = ibuffer;

    for(int i=0; i!=2; ++i ) {
      scnt[i] *= bytes;
      sdsp[i] *= bytes;
      rcnt[i] *= bytes;
      rdsp[i] *= bytes;
    }

    MPI_Alltoallv(&bodies[0].pos[0],scnt,sdsp,MPI_BYTE,
                  &buffer[0].pos[0],rcnt,rdsp,MPI_BYTE,
                  MPI_COMM[2]);
    if( color[0] == color[1] && nthLocal != 0 ) bodies = buffer;
    sort(Ibody,bodies,buffer);
  }

  void multisectionScatter(Bigints &Ibody, Bodies &buffer,
                           int nthLocal, int &newSize) {
    int const MPI_TYPE = getType(Ibody[0]);
    int const bytes = sizeof(bodies[0]);
    int numScatter = nprocs[1] - 1;
    int oldSize = newSize;
    int *scnt = new int [nprocs[1]];
    int *sdsp = new int [nprocs[1]];
    int rcnt;
    if( key[1] == numScatter ) {
      sdsp[0] = 0;
      for(int i=0; i!=numScatter; ++i ) {
        scnt[i] = nthLocal / numScatter;
        sdsp[i+1] = sdsp[i] + scnt[i];
      }
      scnt[numScatter] = 0;
      oldSize = 0;
      newSize -= sdsp[numScatter];
      ibuffer.erase(ibuffer.begin(),ibuffer.begin()+sdsp[numScatter]);
       buffer.erase( buffer.begin(), buffer.begin()+sdsp[numScatter]);
    }
    MPI_Scatter(scnt,1,MPI_INT,&rcnt,1,MPI_INT,numScatter,MPI_COMM[1]);
    if( key[1] != numScatter ) {
      newSize += rcnt;
      ibuffer.resize(newSize);
       buffer.resize(newSize);
    }

    MPI_Scatterv(&Ibody[0],        scnt,sdsp,MPI_TYPE,
                 &ibuffer[oldSize],rcnt,     MPI_TYPE,
                 numScatter,MPI_COMM[1]);
    Ibody = ibuffer;

    for(int i=0; i!= nprocs[1]; ++i ) {
      scnt[i] *= bytes;
      sdsp[i] *= bytes;
    }
    rcnt *= bytes;

    MPI_Scatterv(&bodies[0].pos[0],      scnt,sdsp,MPI_BYTE,
                 &buffer[oldSize].pos[0],rcnt,     MPI_BYTE,
                 numScatter,MPI_COMM[1]);
    bodies = buffer;
    if( key[1] != numScatter ) sort(Ibody,bodies,buffer);
    delete[] scnt;
    delete[] sdsp;
  }

  void multisectionGather(Bigints &Ibody, Bodies &buffer,
                          int nthLocal, int &newSize) {
    int const MPI_TYPE = getType(Ibody[0]);
    int const bytes = sizeof(bodies[0]);
    int numGather = nprocs[0] - 1;
    int oldSize = newSize;
    int scnt;
    int *rcnt = new int [nprocs[0]];
    int *rdsp = new int [nprocs[0]];
    if( key[0] != 0 ) {
      scnt = nthLocal / numGather;
      newSize -= scnt;
      ibuffer.erase(ibuffer.begin(),ibuffer.begin()+scnt);
       buffer.erase( buffer.begin(), buffer.begin()+scnt);
    }
    MPI_Gather(&scnt,1,MPI_INT,rcnt,1,MPI_INT,0,MPI_COMM[0]);
    if( key[0] == 0 ) {
      rdsp[0] = 0;
      for(int i=0; i!=numGather; ++i )
        rdsp[i+1] = rdsp[i] + rcnt[i];
      newSize += rdsp[numGather] + rcnt[numGather];
      ibuffer.resize(newSize);
       buffer.resize(newSize);
    }

    MPI_Gatherv(&Ibody[0],        scnt,     MPI_TYPE,
                &ibuffer[oldSize],rcnt,rdsp,MPI_TYPE,
                0,MPI_COMM[0]);
    Ibody = ibuffer;

    scnt *= bytes;
    for(int i=0; i!= nprocs[0]; ++i ) {
      rcnt[i] *= bytes;
      rdsp[i] *= bytes;
    }

    MPI_Gatherv(&bodies[0].pos[0],      scnt,     MPI_BYTE,
                &buffer[oldSize].pos[0],rcnt,rdsp,MPI_BYTE,
                0,MPI_COMM[0]);
    bodies = buffer;
    delete[] rcnt;
    delete[] rdsp;
    if( key[0] == 0 ) sort(Ibody,bodies,buffer);
  }

  void multisection(Bigints &Ibody, Bodies &buffer) {
    int const MPI_TYPE = getType(Ibody[0]);
    nprocs[0] = nprocs[1] = SIZE;                               // Initialize number of processes in groups
    offset[0] = offset[1] = 0;                                  // Initialize offset of body in groups
     color[0] =  color[1] =  color[2] = 0;                      // Initialize color of communicators
       key[0] =    key[1] =    key[2] = 0;                      // Initialize key of communicators
    int newSize;
    bigint numGlobal;
    bigint numLocal = bodies.size();
    MPI_Reduce(&numLocal,&numGlobal,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Bcast(&numGlobal,1,MPI_TYPE,0,MPI_COMM_WORLD);
    bigint nthGlobal = (numGlobal * (nprocs[0] / 2)) / nprocs[0];
    binBodies(Ibody,0);
    sort(Ibody,bodies,buffer);
    bigint iSplit = nth_element(&Ibody[0],numLocal,nthGlobal);
    int nthLocal = splitBodies(Ibody,iSplit);

    for( int l=0; l!=LEVEL; ++l ) {
      multisectionGetComm(l);

      splitDomain(iSplit,l,l%3);

      multisectionAlltoall(Ibody,buffer,nthLocal,numLocal,newSize);

      if( oldnprocs % 2 == 1 && oldnprocs != 1 && nprocs[0] <= nprocs[1] )
        multisectionScatter(Ibody,buffer,nthLocal,newSize);
      if( oldnprocs % 2 == 1 && oldnprocs != 1 && nprocs[0] >= nprocs[1] )
        multisectionGather(Ibody,buffer,nthLocal,newSize);

#ifdef DEBUG
      for(int j=0; j!=SIZE; ++j) {
        MPI_Barrier(MPI_COMM_WORLD);
        usleep(100);
        if( RANK == j ) {
          std::cout << "RANK : " << j << std::endl;
          for(int i=0; i!=9; ++i) std::cout << Ibody[numLocal/10*i] << " ";
          std::cout << std::endl;
        }
      }
#endif
      numLocal = newSize;
      MPI_Reduce(&numLocal,&numGlobal,1,MPI_TYPE,MPI_SUM,0,MPI_COMM[0]);
      MPI_Bcast(&numGlobal,1,MPI_TYPE,0,MPI_COMM[0]);
      nthGlobal = (numGlobal * (nprocs[0] / 2)) / nprocs[0];
      binBodies(Ibody,(l+1) % 3);
      sort(Ibody,bodies,buffer);
      iSplit = nth_element(&Ibody[0],numLocal,nthGlobal,MPI_COMM[0]);
      nthLocal = splitBodies(Ibody,iSplit);
    }
  }

};

#endif
