#ifndef partition_h
#define partition_h
#include "mympi.h"
#include "types.h"

class Partition : public MyMPI {
private:
  Bodies &bodies;                                               // Bodies to be partitioned
  vect XMIN;                                                    // Minimum position vector of bodies
  vect XMAX;                                                    // Maximum position vector of bodies
public:
  Partition(Bodies &b) : bodies(b) {}                           // Constructor
  ~Partition() {}                                               // Destructor

  void setDomain() {                                            // Set bounds of domain to be partitioned
    B_iter B = bodies.begin();                                  // Reset body iterator
    XMIN = XMAX = B->pos;                                       // Initialize xmin,xmax
    int const MPI_TYPE = getType(XMIN[0]);                      // Get MPI data type
    for( B=bodies.begin(); B!=bodies.end(); ++B ) {             // Loop over all bodies
      for( int d=0; d!=3; ++d ) {                               //  Loop over each dimension
        if     (B->pos[d] < XMIN[d]) XMIN[d] = B->pos[d];       //   Determine xmin
        else if(B->pos[d] > XMAX[d]) XMAX[d] = B->pos[d];       //   Determine xmax
      }                                                         //  End loop over each dimension
    }                                                           // End loop over all bodies
    vect buffer;                                                // Define recv buffer
    MPI_Reduce(XMAX,buffer,3,MPI_TYPE,MPI_MAX,0,MPI_COMM_WORLD);// Reduce global maximum
    XMAX = buffer;                                              // Get data from buffer
    MPI_Bcast(XMAX,3,MPI_TYPE,0,MPI_COMM_WORLD);                // Broadcast global maximum
    MPI_Reduce(XMIN,buffer,3,MPI_TYPE,MPI_MIN,0,MPI_COMM_WORLD);// Reduce global minimum
    XMIN = buffer;                                              // Get data from buffer
    MPI_Bcast(XMIN,3,MPI_TYPE,0,MPI_COMM_WORLD);                // Broadcast global minimum
  }

  void binBodies(Bigints &Ibody, int d) {                       // Turn positions into indices of bins
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over all bodies
      Ibody[B-bodies.begin()] = bigint((B->pos[d] - XMIN[d]) / (XMAX[d] - XMIN[d]) * bodies.size());// Bin bodies
    }                                                           // End loop over all bodies
  }

  void shiftBodies(Bodies &buffer) {                            // Send bodies to next rank (wrap around)
    int newSize;                                                // New number of bodies
    int oldSize = bodies.size();                                // Current number of bodies
    int const bytes = sizeof(bodies[0]);                        // Byte size of body structure
    int const isend = (RANK + 1       ) % SIZE;                 // Send to next rank (wrap around)
    int const irecv = (RANK - 1 + SIZE) % SIZE;                 // Receive from previous rank (wrap around)
    MPI_Request sreq,rreq;                                      // Send, receive request handles

    MPI_Isend(&oldSize,1,MPI_INT,irecv,0,MPI_COMM_WORLD,&sreq); // Send current number of bodies
    MPI_Irecv(&newSize,1,MPI_INT,isend,0,MPI_COMM_WORLD,&rreq); // Receive new number of bodies
    MPI_Wait(&sreq,MPI_STATUS_IGNORE);                          // Wait for send to complete
    MPI_Wait(&rreq,MPI_STATUS_IGNORE);                          // Wait for receive to complete

    buffer.resize(newSize);                                     // Resize buffer to new number of bodies
    MPI_Isend(&bodies[0].pos[0],oldSize*bytes,MPI_CHAR,irecv,   // Send bodies to next rank
              1,MPI_COMM_WORLD,&sreq);
    MPI_Irecv(&buffer[0].pos[0],newSize*bytes,MPI_CHAR,isend,   // Receive bodies from previous rank
              1,MPI_COMM_WORLD,&rreq);
    MPI_Wait(&sreq,MPI_STATUS_IGNORE);                          // Wait for send to complete
    MPI_Wait(&rreq,MPI_STATUS_IGNORE);                          // Wait for receive to complete
    bodies = buffer;                                            // Copy bodies from buffer
  }

  int splitBodies(Bigints &Ibody, bigint iSplit) {
    int nth(0);
    while( Ibody[nth] < iSplit ) nth++;
    return nth;
  }

  void multisection(Bigints &Ibody) {
    int level = int(log(SIZE) / M_LN2 - 1e-5) + 1;
    int nprocs[2] = {0, 0};
    int offset[2] = {0, 0};
    int  color[3] = {0, 0, 0};
    int    key[3] = {0, 0, 0};
    int *scnt  = new int [SIZE];
    int *sdsp  = new int [SIZE];
    int *rcnt  = new int [SIZE];
    int *rdsp  = new int [SIZE];
    bigint nthGlobal = bodies.size() * (SIZE / 2);
    bigint iSplit = nth_element(&Ibody[0],bodies.size(),nthGlobal);
    int nthLocal = splitBodies(Ibody,iSplit);
    nprocs[0] = SIZE;
    for( int l=0; l!=level; ++l ) {
      int oldnprocs = nprocs[0];
      int  oldcolor = color[0];
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
      color[2] = key[1] + oldcolor * (1 << (level - l - 1));
      MPI_Comm MPI_COMM[3];
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
    delete[] scnt;
    delete[] sdsp;
    delete[] rcnt;
    delete[] rdsp;
  }

};

#endif
