#ifndef tree_mpi_h
#define tree_mpi_h
#include "kernel.h"
#include "logger.h"
#include "sort.h"
#include <numeric>

#if EXAFMM_COUNT_KERNEL
#define countKernel(N) N++
#else
#define countKernel(N)
#endif

#if EXAFMM_OVERLAP_REMOTE
#include <omp.h>
#define PRAGMA_OMP(x)            _Pragma( #x )
#define get_num_threads()        omp_get_max_threads() 
#define get_thread_num()         omp_get_thread_num()
#define set_num_threads(threads) omp_set_num_threads(threads)
#define wait_omp_tasks           PRAGMA_OMP(omp taskwait)
#define create_omp_task(E)       PRAGMA_OMP(omp task) E()
#else 
#define get_num_threads()        1 
#define get_thread_num()         0
#define set_num_threads(threads) 
#define wait_omp_tasks           
#define create_omp_task(E)       E()  
#endif

namespace exafmm {
  struct Multipole {    
    int      ICHILD;                                            //!< Index of first child cell
    int      NCHILD;                                            //!< Number of child cells
    int      IBODY;                                             //!< Index of first body
    int      NBODY;                                             //!< Number of descendant bodies
    real_t   SCALE;                                             //!< Scale for Helmholtz kernel
    vec3     X;                                                 //!< Cell center
    real_t   R;                                                 //!< Cell radius
    vecP     M;                                                 //!< Multipole coefficients
  public:    
    Multipole() { }
    Multipole(C_iter const& cc):
    ICHILD(cc->ICHILD),NCHILD(cc->NCHILD),IBODY(cc->IBODY),
    NBODY(cc->NBODY),SCALE(cc->SCALE),X(cc->X),R(cc->R),M(cc->M){ }
    Cell cell() const {
      Cell c;
      c.ICHILD = ICHILD;
      c.NCHILD = NCHILD;
      c.IBODY = IBODY;
      c.NBODY = NBODY;
      c.SCALE = SCALE;
      c.X = X;
      c.R = R;
      c.M = M;
      return c;
    }    
    Multipole(Multipole const& multipole) {
      ICHILD = multipole.ICHILD;
      NCHILD = multipole.NCHILD;
      IBODY = multipole.IBODY;
      NBODY = multipole.NBODY;
      SCALE = multipole.SCALE;
      X = multipole.X;
      R = multipole.R;
      M = multipole.M;
    }
  };
  typedef std::vector<Multipole> Multipoles;
  typedef Multipoles::iterator   M_iter;                        //!< Iterator of Multipole vector
}

namespace exafmm {
//! Handles all the communication of local essential trees
class TreeMPI {
protected:  
  typedef std::vector<std::vector<int> > inter_list;          //!< Remote interaction list type  
  typedef std::pair<int,int>  request_pair;                   //!< Cell request int pair type
  typedef std::pair<C_iter,C_iter>  crequest_pair;            //!< Cell request C_iter pair type  
  const int mpirank;                                          //!< Rank of MPI communicator
  const int mpisize;                                          //!< Size of MPI communicator
  const int images;                                           //!< Number of periodic image sublevels
  const int granularity;                                      //!< The granularity of communication
  const uint32_t nulltag;                                     //!< Tag of null message type
  const uint32_t celltag;                                     //!< Tag of cell message type
  const uint32_t childcelltag;                                //!< Tag of child cell message type
  const uint32_t bodytag;                                     //!< Tag of body message type
  const uint32_t flushtag;                                    //!< Tag of flush message type
  const uint32_t maxtag;                                      //!< Max value of tag
  const uint32_t levelshift;                                  //!< Number of level bit shifts
  const uint32_t requestshift;                                //!< Number of request bit shifts
  const uint32_t directionshift;                              //!< Number of directions bit shifts
  const uint32_t levelmask;                                   //!< Mask value of level msg
  const uint32_t requestmask;                                 //!< Mask value of request msg
  const uint32_t directionmask;                               //!< Mask value of direction msg
  const uint32_t grainmask;                                   //!< Mask value of grain msg
  const int    nspawn;                                        //!< Threshold of NBODY for spawning new threads
  const char sendbit;                                         //!< Send bit value
  const char receivebit;                                      //!< Recv bit value
  float (* allBoundsXmin)[3];                                 //!< Array for local Xmin for all ranks
  float (* allBoundsXmax)[3];                                 //!< Array for local Xmax for all ranks
  Bodies sendBodies;                                          //!< Send buffer for bodies
  Bodies recvBodies;                                          //!< Receive buffer for bodies
  Cells sendCells;                                            //!< Send buffer for cells
  Cells recvCells;                                            //!< Receive buffer for cells
  Multipoles* recvRoots;                                      //!< Receive buffer for tree roots
  int * sendBodyCount;                                        //!< Send count
  int * sendBodyDispl;                                        //!< Send displacement
  int * recvBodyCount;                                        //!< Receive count
  int * recvBodyDispl;                                        //!< Receive displacement
  int * sendCellCount;                                        //!< Send count
  int * sendCellDispl;                                        //!< Send displacement
  int * recvCellCount;                                        //!< Receive count
  int * recvCellDispl;                                        //!< Receive displacement  
  C_iter Ci0;                                                 //!< Iterator of first target cell
  C_iter Cj0;                                                 //!< Iterator of first source cell
  int terminated;                                             //!< Number of terminated ranks
  int hitCount;                                               //!< Number of remote requests  
  int cellWordSize;                                           //!< Number of words in cell struct
  int bodyWordSize;                                           //!< Number of words in body struct  
  std::vector<MPI_Request*> pendingCellMPI_Requests;          //!< Buffer for non-blocking requests
  std::vector<Multipoles*>  sendCellBuffers;                  //!< Buffer for non-blocking cell sends
  std::vector<MPI_Request*> pendingBodyMPI_Requests;          //!< Buffer for non-blocking requests
  std::vector<Bodies*> sendBodyBuffers;                       //!< Buffer for non-blocking body sends
  std::vector<int> bodyReferenceCount;                        //!< Reference count for received LET bodies
  std::vector<int> cellReferenceCount;                        //!< Reference count for received LET cells
  int* childCellsReferenceCount;                              //!< Reference count for received LET child cells
  inter_list remoteInteractionList;                           //!< Remote interaction for each leaf cell
  double remoteCommunicationTime;                             //!< Time spent doing remote traversal 
  std::vector<crequest_pair > aggregateP2PPairs;              //!< A buffer to aggregate body request calls
  std::vector<crequest_pair > aggregateCellPairs;             //!< A buffer to aggregate body request calls
  int aggregateP2PCount;                                      //!< A counter to maintain the number of queued leaf cell requests
  int *bodySendRequestBuffer;                                 //!< Send request buffer for aggregated leaf cells
  int *bodyReceiveRequestBuffer;                              //!< Recv request buffer for aggregated leaf cells
#if EXAFMM_OVERLAP_REMOTE  
  std::vector<crequest_pair >* bodyRequests;                  //!< A lock free body requests buffer  
  std::vector<crequest_pair >* cellRequests;                  //!< A lock free cell requests buffer    
  std::vector<bool> fullfilledBodyRequests;                   //!< Reference count for received LET bodies
  std::vector<bool> fullfilledCellRequests;                   //!< Reference count for received LET cells
  pthread_mutex_t* mutex;                                     //!< Pthread mutex
  int threadCount;                                            //!< Number of threads
  volatile bool flushAggregate;                               //!< Trigger to flush aggregate P2P      
  volatile bool finalizeCommunicationThread;                  //!< Trigger end of communication thread   
#endif
#if EXAFMM_COUNT_KERNEL
  real_t numP2P;                                              //!< Number of P2P kernel calls
  real_t numM2L;                                              //!< Number of M2L kernel calls
  real_t numM2LMsgs;                                          //!< Number of M2L cell messages 
  real_t numP2PMsgs;                                          //!< Number of Leaf cells (P2P) messages 
#endif
private:

  //! Exchange send count for bodies
  void alltoall(Bodies & bodies) {
    for (int i = 0; i < mpisize; i++) {                       // Loop over ranks
      sendBodyCount[i] = 0;                                   //  Initialize send counts
    }                                                         // End loop over ranks
    for (B_iter B = bodies.begin(); B != bodies.end(); B++) { // Loop over bodies
      assert(0 <= B->IRANK && B->IRANK < mpisize);            //  Check bounds for process ID
      sendBodyCount[B->IRANK]++;                              //  Fill send count bucket
      //B->IRANK = mpirank;                                     //  Tag for sending back to original rank
    }                                                         // End loop over bodies
    MPI_Alltoall(sendBodyCount, 1, MPI_INT,                   // Communicate send count to get receive count
                 recvBodyCount, 1, MPI_INT, MPI_COMM_WORLD);
    sendBodyDispl[0] = recvBodyDispl[0] = 0;                  // Initialize send/receive displacements
    for (int irank = 0; irank < mpisize - 1; irank++) {       // Loop over ranks
      sendBodyDispl[irank + 1] = sendBodyDispl[irank] + sendBodyCount[irank]; //  Set send displacement
      recvBodyDispl[irank + 1] = recvBodyDispl[irank] + recvBodyCount[irank]; //  Set receive displacement
    }                                                         // End loop over ranks
  }

  //! Exchange bodies
  void alltoallv(Bodies & bodies) {
    assert( (sizeof(bodies[0]) & 3) == 0 );                   // Body structure must be 4 Byte aligned
    int word = sizeof(bodies[0]) / 4;                         // Word size of body structure
    recvBodies.resize(recvBodyDispl[mpisize - 1] + recvBodyCount[mpisize - 1]); // Resize receive buffer
    for (int irank = 0; irank < mpisize; irank++) {           // Loop over ranks
      sendBodyCount[irank] *= word;                           //  Multiply send count by word size of data
      sendBodyDispl[irank] *= word;                           //  Multiply send displacement by word size of data
      recvBodyCount[irank] *= word;                           //  Multiply receive count by word size of data
      recvBodyDispl[irank] *= word;                           //  Multiply receive displacement by word size of data
    }                                                         // End loop over ranks
    MPI_Alltoallv((int*)&bodies[0], sendBodyCount, sendBodyDispl, MPI_INT,// Communicate bodies
                  (int*)&recvBodies[0], recvBodyCount, recvBodyDispl, MPI_INT, MPI_COMM_WORLD);
    for (int irank = 0; irank < mpisize; irank++) {           // Loop over ranks
      sendBodyCount[irank] /= word;                           //  Divide send count by word size of data
      sendBodyDispl[irank] /= word;                           //  Divide send displacement by word size of data
      recvBodyCount[irank] /= word;                           //  Divide receive count by word size of data
      recvBodyDispl[irank] /= word;                           //  Divide receive displacement by word size of data
    }                                                         // End loop over ranks
  }

  //! Exchange bodies in a point-to-point manner
  template<typename T>
  void alltoallv_p2p(T& sendB, T& recvB, int* recvDispl, int* recvCount, int* sendDispl, int* sendCount) {        
    int dataSize = recvDispl[mpisize - 1] + recvCount[mpisize - 1];    
    if (sendB.size() == sendCount[mpirank] && dataSize == sendCount[mpirank]) {
      recvB = sendB;
      return;
    }    
    int* granularRecvDispl = new int [mpisize];
    int* granularSendDispl = new int [mpisize];
    for (int i = 0; i < mpisize; ++i) {
      granularSendDispl[i] = sendDispl[i];
      granularRecvDispl[i] = recvDispl[i];
    }
    recvB.resize(dataSize);
    assert( (sizeof(sendB[0]) & 3) == 0 );
    int word = sizeof(sendB[0]) / 4;
    int* sendBuff = (int*)&sendB[0];
    int* recvBuff = (int*)&recvB[0];
    int grain_size = granularity;
    if(grain_size <= 0) grain_size = 1;    
    for (int i = 0; i < grain_size; ++i) {      
      int sendSize = 0;
      int recvSize = 0;
      MPI_Request* rreq   = new MPI_Request[mpisize - 1];
      MPI_Request* sreq   = new MPI_Request[mpisize - 1];
      MPI_Status* rstatus = new MPI_Status[mpisize - 1];
      MPI_Status* sstatus = new MPI_Status[mpisize - 1];
      for (int irank = 0; irank < mpisize; ++irank) {
        if (irank != mpirank) {
          if (recvCount[irank] > 0) {
            int gRecvCount = recvCount[irank]/grain_size;
            if(i == grain_size - 1) gRecvCount += recvCount[irank]%grain_size;            
            MPI_Irecv(recvBuff + granularRecvDispl[irank]*word,
                      gRecvCount*word, MPI_INT, irank, irank, MPI_COMM_WORLD, &rreq[recvSize]);
            granularRecvDispl[irank]+=gRecvCount;
            recvSize++;
          }
        }
      }
      for (int irank = 0; irank < mpisize; ++irank) {
        if (irank != mpirank) {
          if (sendCount[irank] > 0) {
            int gSendCount = sendCount[irank]/grain_size;
            if(i == grain_size - 1) gSendCount += sendCount[irank]%grain_size;
            MPI_Isend(sendBuff + granularSendDispl[irank]*word,
                      gSendCount*word, MPI_INT, irank, mpirank, MPI_COMM_WORLD, &sreq[sendSize]);
            granularSendDispl[irank]+=gSendCount;
            sendSize++;
          }
        }
      }          
      MPI_Waitall(sendSize, &sreq[0], &sstatus[0]);
      MPI_Waitall(recvSize, &rreq[0], &rstatus[0]);
      delete[] rreq;
      delete[] sreq;
      delete[] rstatus;
      delete[] sstatus;      
    } 
    delete[] granularRecvDispl;
    delete[] granularSendDispl;
    typename T::iterator localBuffer = sendB.begin() + sendDispl[mpirank];
    std::copy(localBuffer, localBuffer + sendCount[mpirank], recvB.begin() + recvBodyDispl[mpirank]);   
  }

  //! Exchange send count for cells
  void alltoall(Cells) {
    MPI_Alltoall(sendCellCount, 1, MPI_INT,                   // Communicate send count to get receive count
                 recvCellCount, 1, MPI_INT, MPI_COMM_WORLD);
    recvCellDispl[0] = 0;                                     // Initialize receive displacements
    for (int irank = 0; irank < mpisize - 1; irank++) {       // Loop over ranks
      recvCellDispl[irank + 1] = recvCellDispl[irank] + recvCellCount[irank]; //  Set receive displacement
    }                                                         // End loop over ranks
  }

  //! Exchange cells
  void alltoallv(Cells & cells) {
    assert( (sizeof(cells[0]) & 3) == 0 );                    // Cell structure must be 4 Byte aligned
    int word = sizeof(cells[0]) / 4;                          // Word size of body structure
    recvCells.resize(recvCellDispl[mpisize - 1] + recvCellCount[mpisize - 1]); // Resize receive buffer
    for (int irank = 0; irank < mpisize; irank++) {           // Loop over ranks
      sendCellCount[irank] *= word;                           //  Multiply send count by word size of data
      sendCellDispl[irank] *= word;                           //  Multiply send displacement by word size of data
      recvCellCount[irank] *= word;                           //  Multiply receive count by word size of data
      recvCellDispl[irank] *= word;                           //  Multiply receive displacement by word size of data
    }                                                         // End loop over ranks
    MPI_Alltoallv((int*)&cells[0], sendCellCount, sendCellDispl, MPI_INT,// Communicate cells
                  (int*)&recvCells[0], recvCellCount, recvCellDispl, MPI_INT, MPI_COMM_WORLD);
    for (int irank = 0; irank < mpisize; irank++) {           // Loop over ranks
      sendCellCount[irank] /= word;                           //  Divide send count by word size of data
      sendCellDispl[irank] /= word;                           //  Divide send displacement by word size of data
      recvCellCount[irank] /= word;                           //  Divide receive count by word size of data
      recvCellDispl[irank] /= word;                           //  Divide receive displacement by word size of data
    }                                                         // End loop over ranks
  }

  Multipoles alltoallv(Multipoles & multipoles) {
    assert( (sizeof(multipoles[0]) & 3) == 0 );               // Cell structure must be 4 Byte aligned
    int word = sizeof(multipoles[0]) / 4;                     // Word size of body structure
    Multipoles recvMults(recvCellDispl[mpisize - 1] + recvCellCount[mpisize - 1]); // Resize receive buffer
    for (int irank = 0; irank < mpisize; irank++) {           // Loop over ranks
      sendCellCount[irank] *= word;                           //  Multiply send count by word size of data
      sendCellDispl[irank] *= word;                           //  Multiply send displacement by word size of data
      recvCellCount[irank] *= word;                           //  Multiply receive count by word size of data
      recvCellDispl[irank] *= word;                           //  Multiply receive displacement by word size of data
    }                                                         // End loop over ranks
    MPI_Alltoallv((int*)&multipoles[0], sendCellCount, sendCellDispl, MPI_INT,// Communicate cells
                  (int*)&recvMults[0], recvCellCount, recvCellDispl, MPI_INT, MPI_COMM_WORLD);
    for (int irank = 0; irank < mpisize; irank++) {           // Loop over ranks
      sendCellCount[irank] /= word;                           //  Divide send count by word size of data
      sendCellDispl[irank] /= word;                           //  Divide send displacement by word size of data
      recvCellCount[irank] /= word;                           //  Divide receive count by word size of data
      recvCellDispl[irank] /= word;                           //  Divide receive displacement by word size of data
    }                                                         // End loop over ranks
    return recvMults;
  }

  inline bool processIncomingMessage(int tag, int source) {
    if ((tag & directionmask) == receivebit)
      return false;    
    uint8_t msgType = getMessageType(tag);
    int nullTag = encryptMessage(nulltag, 0, receivebit);
    MPI_Request* request = new MPI_Request();
    char null;
    bool dataReady = false;
    if (msgType == childcelltag) {
      int childID;
      hitCount++;            
      int level = getMessageLevel(tag);
      MPI_Recv(&childID, 1, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      toggleDirection(tag);      
      Cell& pcell = sendCells[sendCells[childID].IPARENT];
      if (pcell.NCHILD > 0) {
        int sendingSize = 0;
        Multipoles* sendBuffer = new Multipoles();        
        if(granularity > 1) {          
          setLETSubset(sendBuffer, sendCells.begin(), pcell, granularity, sendingSize, source);          
        } else {
          sendingSize = pcell.NCHILD;
          C_iter C0 = sendCells.begin() + pcell.ICHILD;
          for(int i = 0; i < sendingSize; ++ i) sendBuffer->push_back(Multipole(C0+i));                  
        }
        MPI_Isend((int*) & (*sendBuffer)[0],  cellWordSize * sendingSize, MPI_INT, source, tag, MPI_COMM_WORLD, request);
        sendCellBuffers.push_back(sendBuffer);
        pendingCellMPI_Requests.push_back(request);          
        dataReady = true;
      }      
    } else if (msgType == bodytag) {
      hitCount++;            
      MPI_Recv(bodyReceiveRequestBuffer, 2*granularity, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      toggleDirection(tag);
      int ibody,nbody,sendingSize=0;
      Bodies* sendBuffer = new Bodies();
      for (int i = 0; i < granularity; ++i) {
        ibody = bodyReceiveRequestBuffer[i*2];
        if(ibody >=0) {
          nbody = bodyReceiveRequestBuffer[i*2+1];
          sendBuffer->insert(sendBuffer->end(),sendBodies.begin()+ibody,sendBodies.begin()+ibody+nbody);
          sendingSize+=nbody;
        }
      }
      MPI_Isend((int*) & (*sendBuffer)[0], bodyWordSize * sendingSize, MPI_INT, source, tag, MPI_COMM_WORLD, request);      
      sendBodyBuffers.push_back(sendBuffer);
      pendingBodyMPI_Requests.push_back(request);          
      dataReady = true;
    } else if (msgType == celltag) {
      int cellID;
      hitCount++;            
      MPI_Recv(&cellID, 1, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      toggleDirection(tag); 
      Cell& pcell = sendCells[cellID];           
      if(granularity > 1 && pcell.NCHILD > 0) {        
        int sendingSize = 1;
        Multipoles* sendBuffer = new Multipoles();
        sendBuffer->push_back(Multipole(sendCells.begin()+cellID));        
        setLETSubset(sendBuffer, sendCells.begin(), pcell, granularity, sendingSize, source);
        MPI_Isend((int*) & (*sendBuffer)[0],  cellWordSize * sendingSize, MPI_INT, source, tag, MPI_COMM_WORLD, request);
        sendCellBuffers.push_back(sendBuffer);
        pendingCellMPI_Requests.push_back(request);          
        dataReady = true;
      } else {
        Multipoles* sendBuffer = new Multipoles();
        sendBuffer->push_back(Multipole(sendCells.begin()+cellID));
        MPI_Isend((int*) & (*sendBuffer)[0], cellWordSize, MPI_INT, source, tag, MPI_COMM_WORLD, request);
        sendCellBuffers.push_back(sendBuffer);
        pendingCellMPI_Requests.push_back(request);        
        dataReady = true;
      }
    }
    else if (msgType == flushtag) {
      MPI_Recv(&null, 1, MPI_CHAR, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      terminated++;
      dataReady = true;
    }
    else {
      MPI_Isend(&null, 1, MPI_CHAR, source, nullTag, MPI_COMM_WORLD, request);
      dataReady = true;
    }
    deallocateCompletedRequests(pendingCellMPI_Requests,sendCellBuffers);
    deallocateCompletedRequests(pendingBodyMPI_Requests,sendBodyBuffers);    
    return dataReady;
  }

  template<typename MPI_T, typename BUFF_T>
  void deallocateCompletedRequests(MPI_T& requests, BUFF_T& data) {
    int count = requests.size();
    for (int i = 0; i < count; ++i) {
      if (*(requests[i]) != MPI_REQUEST_NULL) {
        int flag;
        MPI_Test(requests[i], &flag, MPI_STATUS_IGNORE);
        if (flag) {
          delete data[i];
          requests.erase(requests.begin() + i);
          data.erase(data.begin() + i);
          count--;
        }
      } else {
        delete data[i];
        requests.erase(requests.begin() + i);
        data.erase(data.begin() + i);
        count--;
      }
    }
  }

protected:
  //! Get distance to other domain
  real_t getDistance(C_iter C, Bounds bounds, vec3 Xperiodic) {
    vec3 dX;                                                  // Distance vector
    for (int d = 0; d < 3; d++) {                             // Loop over dimensions
      dX[d] = (C->X[d] + Xperiodic[d] > bounds.Xmax[d]) *     //  Calculate the distance between cell C and
              (C->X[d] + Xperiodic[d] - bounds.Xmax[d]) +           //  the nearest point in domain [xmin,xmax]^3
              (C->X[d] + Xperiodic[d] < bounds.Xmin[d]) *           //  Take the differnece from xmin or xmax
              (C->X[d] + Xperiodic[d] - bounds.Xmin[d]);            //  or 0 if between xmin and xmax
    }                                                         // End loop over dimensions
    return norm(dX);                                          // Return distance squared
  }

  //! Get distance to other domain
  inline real_t getRankDistance(C_iter C, vec3 Xperiodic, int rank) {
    vec3 dX;                                                  // Distance vector
    for (int d = 0; d < 3; d++) {                             // Loop over dimensions
      dX[d] = (C->X[d] + Xperiodic[d] > allBoundsXmax[rank][d]) * //  Calculate the distance between cell C and
              (C->X[d] + Xperiodic[d] - allBoundsXmax[rank][d]) +   //  the nearest point in domain [xmin,xmax]^3
              (C->X[d] + Xperiodic[d] < allBoundsXmin[rank][d]) *   //  Take the differnece from xmin or xmax
              (C->X[d] + Xperiodic[d] - allBoundsXmin[rank][d]);    //  or 0 if between xmin and xmax
    }                                                         // End loop over dimensions
    return norm(dX);                                          // Return distance squared
  }

  //! Dual tree traversal for a single pair of cells
  void traverseRemote(C_iter Ci, C_iter Cj, bool mutual, real_t remote, int rank) {
    vec3 dX = Ci->X - Cj->X - kernel::Xperiodic;                // Distance vector from source to target
    real_t R2 = norm(dX);                                       // Scalar distance squared
    if (R2 > (Ci->R + Cj->R) * (Ci->R + Cj->R) * (1 - 1e-3)) {  // Distance is far enough
      kernel::M2L(Ci, Cj, false);                               //  M2L kernel            
      countKernel(numM2L);
    } else if (Ci->NCHILD == 0 && Cj->NCHILD == 0) {            // Else if both cells are bodies
      if (Cj->NBODY == 0) {
        kernel::M2L(Ci, Cj, false);                             //   M2L kernel 
        countKernel(numM2L);            
    } else {             
        countKernel(numP2P); 
        P2P(Ci,Cj,rank);            
      }
    } else {                                                    // Else if cells are close but not bodies
      splitCellRemote(Ci, Cj, mutual, remote, rank);            //  Split cell and call function recursively for child
    }                                                           // End if for multipole acceptance
  }

  //! Split cell and call traverse() recursively for child
  void splitCellRemote(C_iter Ci, C_iter Cj, bool mutual, real_t remote, int rank) {    
    bool ready;
    if (Cj->NCHILD == 0) {                                      // If Cj is leaf
      assert(Ci->NCHILD > 0);                                   //  Make sure Ci is not leaf
      for (C_iter ci = Ci0 + Ci->ICHILD; ci != Ci0 + Ci->ICHILD + Ci->NCHILD; ci++) { // Loop over Ci's children
        traverseRemote(ci, Cj, mutual, remote, rank);           //   Traverse a single pair of cells
      }                                                         //  End loop over Ci's children
    } else if (Ci->NCHILD == 0) {                               // Else if Ci is leaf
      assert(Cj->NCHILD > 0);                                   //  Make sure Cj is not leaf
      C_iter cj0 = getChildCells(Ci,Cj,rank,ready);if(!ready) return;
      for (C_iter cj=cj0; cj != cj0 + Cj->NCHILD; ++cj) {       // Loop over Cj's children
        traverseRemote(Ci, cj, mutual, remote, rank);           //   Traverse a single pair of cells
      }
    } else if (Ci->NBODY + Cj->NBODY >= nspawn) {               // Else if cells are still large      
      C_iter cj  = getChildCells(Ci,Cj,rank,ready);if(!ready) return;
      TraverseRemoteRange traverseRange(this, Ci0 + Ci->ICHILD, Ci0 + Ci->ICHILD + Ci->NCHILD, // Instantiate recursive functor
                                        cj, cj + Cj->NCHILD, mutual, remote, rank);
      traverseRange();
    } else if (Ci->R >= Cj->R) {                                // Else if Ci is larger than Cj
      for (C_iter ci = Ci0 + Ci->ICHILD; ci != Ci0 + Ci->ICHILD + Ci->NCHILD; ci++) { // Loop over Ci's children
        traverseRemote(ci, Cj, mutual, remote, rank);           //   Traverse a single pair of cells
      }                                                         //  End loop over Ci's children
    } else {                                                    // Else if Cj is larger than Ci
      C_iter cj0 = getChildCells(Ci,Cj,rank,ready); if(!ready) return;
      for (C_iter cj=cj0; cj != cj0 + Cj->NCHILD; ++cj) {       // Loop over Cj's children
        traverseRemote(Ci, cj, mutual, remote, rank);           //   Traverse a single pair of cells
      }                                                         //  End loop over Cj's children
    }                                                           // End if for leafs and Ci Cj size
  }

  //! Recursive functor for dual tree traversal of a range of Ci and Cj
  struct TraverseRemoteRange {
    TreeMPI * treeMPI;                                          //!< TreeMPI object
    C_iter CiBegin;                                             //!< Begin iterator of target cells
    C_iter CiEnd;                                               //!< End iterator of target cells
    C_iter CjBegin;                                             //!< Begin Iterator of source cells
    C_iter CjEnd;                                               //!< End iterator of source cells
    bool mutual;                                                //!< Flag for mutual interaction
    real_t remote;                                              //!< Weight for remote work load    
    int rank;
    TraverseRemoteRange(TreeMPI * _treeMPI, C_iter _CiBegin, C_iter _CiEnd,// Constructor
                        C_iter _CjBegin, C_iter _CjEnd,
                        bool _mutual, real_t _remote, int _rank) :
      treeMPI(_treeMPI), CiBegin(_CiBegin), CiEnd(_CiEnd),      // Initialize variables
      CjBegin(_CjBegin), CjEnd(_CjEnd),
      mutual(_mutual), remote(_remote), rank(_rank) {}
    void operator() () {                                        // Overload operator()
      Tracer tracer;                                            //  Instantiate tracer
      logger::startTracer(tracer);                              //  Start tracer
      if (CiEnd - CiBegin == 1 || CjEnd - CjBegin == 1) {       //  If only one cell in range
        for (C_iter Ci = CiBegin; Ci != CiEnd; Ci++) {          //    Loop over all Ci cells
          for (C_iter Cj = CjBegin; Cj != CjEnd; Cj++) {        //     Loop over all Cj cells
            treeMPI->traverseRemote(Ci, Cj, mutual, remote, rank); // Call traverse for single pair
          }                                                     //     End loop over all Cj cells
        }                                                       //   End if for Ci == Cj
      } else {                                                  //  If many cells are in the range
        C_iter CiMid = CiBegin + (CiEnd - CiBegin) / 2;         //   Split range of Ci cells in half
        C_iter CjMid = CjBegin + (CjEnd - CjBegin) / 2;         //   Split range of Cj cells in half
        {
          TraverseRemoteRange leftBranch(treeMPI, CiBegin, CiMid,//    Instantiate recursive functor
                                         CjBegin, CjMid, mutual, remote, rank);
          create_omp_task(leftBranch);                            
          TraverseRemoteRange rightBranch(treeMPI, CiMid, CiEnd,//    Instantiate recursive functor
                                          CjMid, CjEnd, mutual, remote, rank);;
          rightBranch();                                        //    Ci:latter Cj:latter
          wait_omp_tasks;
        }
        {
          TraverseRemoteRange leftBranch(treeMPI, CiBegin, CiMid,//    Instantiate recursive functor
                                         CjMid, CjEnd, mutual, remote, rank);;          
          create_omp_task(leftBranch);          
          TraverseRemoteRange rightBranch(treeMPI, CiMid, CiEnd,//    Instantiate recursive functor
                                          CjBegin, CjMid, mutual, remote, rank);;
          rightBranch();                                        //    Ci:latter Cj:former
          wait_omp_tasks;
        }
      }                                                         //  End if for many cells in range
      logger::stopTracer(tracer);                               //  Stop tracer
    }                                                           // End overload operator()
  };

  //! Add cells to send buffer
  void addSendCell(C_iter C, int & irank, int & icell, int & iparent, bool copyData) {
    if (copyData) {                                           // If copying data to send cells
      Cell cell(*C);                                          //  Initialize send cell
      cell.NCHILD = cell.NBODY = 0;                           //  Reset counters
      cell.IPARENT = iparent;                                 //  Index of parent
      sendCells[sendCellDispl[irank] + icell] = cell;         //  Copy cell to send buffer
      C_iter Cparent = sendCells.begin() + sendCellDispl[irank] + iparent;// Get parent iterator
      if (Cparent->NCHILD == 0) Cparent->ICHILD = icell;      //  Index of parent's first child
      Cparent->NCHILD++;                                      //  Increment parent's child counter
    }                                                         // End if for copying data to send cells
    icell++;                                                  // Increment cell counter
  }

  //! Add bodies to send buffer
  void addSendBody(C_iter C, int & irank, int & ibody, int icell, bool copyData) {
    if (copyData) {                                           // If copying data to send bodies
      C_iter Csend = sendCells.begin() + sendCellDispl[irank] + icell; // Send cell iterator
      Csend->NBODY = C->NBODY;                                //  Number of bodies
      Csend->IBODY = ibody;                                   //  Body index per rank
      B_iter Bsend = sendBodies.begin() + sendBodyDispl[irank] + ibody; // Send body iterator
      for (B_iter B = C->BODY; B != C->BODY + C->NBODY; B++, Bsend++) { //  Loop over bodies in cell
        *Bsend = *B;                                          //   Copy body to send buffer
        Bsend->IRANK = irank;                                 //   Assign destination rank
      }                                                       //  End loop over bodies in cell
    }                                                         // End if for copying data to send bodies
    ibody += C->NBODY;                                        // Increment body counter
  }

  inline int encryptMessage(uint8_t requestType, uint8_t level, char direction) {
    int tag = requestType;
    tag <<= levelshift; tag |= level;
    tag <<= directionshift; tag |= direction;
    return tag;
  }

  inline void decryptMessage(int tag, uint8_t& requestType, uint8_t& level, char& direction) {
    direction = tag & directionmask;
    tag >>= directionshift;
    level = tag & levelmask;
    tag >>= levelshift;
    requestType = tag & requestmask;
    tag >>= requestshift;
  }

  inline uint8_t getMessageType(int const& tag) {
    return ((tag >> (directionshift + levelshift)) & requestmask);
  }

  inline uint8_t getMessageLevel(int const& tag) {
    return ((tag >> directionshift) & levelmask);
  }

  inline char getMessageDirection(int const& tag) {
    return (tag & directionmask);
  }

  inline void toggleDirection(int& tag) {
    tag ^= directionmask;
  }

  inline void setLETSubset(Cells* cells, C_iter const& C0, Cell const& C, uint16_t const& grainSize, int& index, int rank) {
    C_iter begin = C0 + C.ICHILD;
    C_iter end   = C0 + C.ICHILD + C.NCHILD;
    int icell = index;
    cells->insert(cells->end(), begin, end);
    index += C.NCHILD;
    real_t R2;
    for (C_iter cc = begin; cc < end; ++cc, ++icell)
      if (index < grainSize) {
        R2 = getRankDistance(cc, kernel::Xperiodic, rank);    //       Get distance to other domain
        if (cc->R * cc->R > R2) {                             //       Divide if cell seems too close
          setLETSubset(cells, C0, *cc, grainSize, index, rank);
        } else {
          (*cells)[icell].NCHILD = 0;
          (*cells)[icell].NBODY = 0;          
        }
      }
  }

  inline void setLETSubset(Multipoles* multipoles, C_iter const& C0, Cell const& C, uint16_t const& grainSize, int& index, int rank) {
    C_iter begin = C0 + C.ICHILD;
    C_iter end   = C0 + C.ICHILD + C.NCHILD;
    int icell = index;        
    for(C_iter cc = begin; cc!=end; ++cc) {
      multipoles->push_back(Multipole(cc));
    }
    index += C.NCHILD;
    real_t R2;
    for (C_iter cc = begin; cc < end; ++cc, ++icell) {              
      if (index < grainSize) {
        R2 = getRankDistance(cc, kernel::Xperiodic, rank);    //       Get distance to other domain
        if (cc->R * cc->R > R2) {                             //       Divide if cell seems too close
          setLETSubset(multipoles, C0, *cc, grainSize, index, rank);
        } else {
          (*multipoles)[icell].NCHILD = 0;
          (*multipoles)[icell].NBODY = 0;          
        }
      }
    }
  }

  //! Determine which cells to send
  void traverseLET(C_iter C, C_iter C0, Bounds bounds, vec3 cycle,
                   int & irank, int & ibody, int & icell, int iparent, bool copyData) {
    int level = int(logf(mpisize - 1) / M_LN2 / 3) + 1;       // Level of local root cell
    if (mpisize == 1) level = 0;                              // Account for serial case
    bool divide[8] = {0, 0, 0, 0, 0, 0, 0, 0};                // Initialize divide flag
    int icells[8] = {0, 0, 0, 0, 0, 0, 0, 0};                 // Initialize icell array
    int cc = 0;                                               // Initialize child index
    for (C_iter CC = C0 + C->ICHILD; CC != C0 + C->ICHILD + C->NCHILD; CC++, cc++) { // Loop over child cells
      icells[cc] = icell;                                     //  Store cell index
      addSendCell(CC, irank, icell, iparent, copyData);       //  Add cells to send
      if (CC->NCHILD == 0) {                                  //  If cell is leaf
        addSendBody(CC, irank, ibody, icell - 1, copyData);   //   Add bodies to send
      } else {                                                //  If cell is not leaf
        vec3 Xperiodic = 0;                                   //   Periodic coordinate offset
        if (images == 0) {                                    //   If free boundary condition
          real_t R2 = getDistance(CC, bounds, Xperiodic);     //    Get distance to other domain
          divide[cc] |= 4 * CC->R * CC->R > R2;               //    Divide if the cell seems too close
        } else {                                              //   If periodic boundary condition
          for (int ix = -1; ix <= 1; ix++) {                  //    Loop over x periodic direction
            for (int iy = -1; iy <= 1; iy++) {                //     Loop over y periodic direction
              for (int iz = -1; iz <= 1; iz++) {              //      Loop over z periodic direction
                Xperiodic[0] = ix * cycle[0];                 //       Coordinate offset for x periodic direction
                Xperiodic[1] = iy * cycle[1];                 //       Coordinate offset for y periodic direction
                Xperiodic[2] = iz * cycle[2];                 //       Coordinate offset for z periodic direction
                real_t R2 = getDistance(CC, bounds, Xperiodic); //       Get distance to other domain
                divide[cc] |= 4 * CC->R * CC->R > R2;         //       Divide if cell seems too close
              }                                               //      End loop over z periodic direction
            }                                                 //     End loop over y periodic direction
          }                                                   //    End loop over x periodic direction
        }                                                     //   Endif for periodic boundary condition
        divide[cc] |= CC->R > (max(cycle) / (1 << (level + 1)));   //   Divide if cell is larger than local root cell
      }                                                       //  Endif for leaf
    }                                                         // End loop over child cells
    cc = 0;                                                   // Initialize child index
    for (C_iter CC = C0 + C->ICHILD; CC != C0 + C->ICHILD + C->NCHILD; CC++, cc++) { // Loop over child cells
      if (divide[cc]) {                                       //  If cell must be divided further
        iparent = icells[cc];                                 //   Parent cell index
        traverseLET(CC, C0, bounds, cycle, irank, ibody, icell, iparent, copyData);// Recursively traverse tree to set LET
      }                                                       //  End if for cell division
    }                                                         // End loop over child cells
  }
  
  void appendLET(Cells const& cells, int& index, int id) {
    Cell parent = cells[id];
    if (index < cells.size()) {
      id = index;
      if (parent.NCHILD > 0) {
        std::copy(cells.begin() + index,cells.begin() + index + parent.NCHILD, recvCells.begin() + parent.ICHILD);       
#if EXAFMM_OVERLAP_REMOTE
        fullfilledCellRequests[parent.ICHILD] = true;
        childCellsReferenceCount[parent.ICHILD] = -1;
#else
        childCellsReferenceCount[parent.ICHILD]++;        
#endif        
        index += parent.NCHILD;
        for (int i = 0; i < parent.NCHILD; ++i) {
          appendLET(cells, index, id++);
        }
      }
    }
  }

  void appendLET(Multipoles const& multipoles, int& index, int id) {
    Multipole parent = multipoles[id];
    if (index < multipoles.size()) {
      id = index;
      if (parent.NCHILD > 0) {
        for (int i = 0; i < parent.NCHILD; ++i) {  
          recvCells[parent.ICHILD + i] = multipoles[index + i].cell();       
        }          
#if EXAFMM_OVERLAP_REMOTE
        fullfilledCellRequests[parent.ICHILD] = true;
        childCellsReferenceCount[parent.ICHILD] = -1;
#else
        childCellsReferenceCount[parent.ICHILD]++;        
#endif        
        index += parent.NCHILD;
        for (int i = 0; i < parent.NCHILD; ++i) {
          appendLET(multipoles, index, id++);
        }
      }
    }
  }


  #if EXAFMM_COUNT_LIST 
  inline void updateInteractionList(int const& cellId, int const& localNumP2P, int const& rank) {
    std::vector<int>& interaction = remoteInteractionList[cellId];    
    interaction[mpirank] = localNumP2P;
    interaction[rank]++;      
  }
#else
  inline void updateInteractionList(int const&, int const&, int const&){ }
#endif

  //! Get bodies underlying a key and a level
  B_iter remoteP2P(C_iter Ci, C_iter Cj, int rank, bool execute = true) {
    assert(rank != mpirank);        
    int key = Cj->IBODY;
    int iparent = Ci->IPARENT;
#if EXAFMM_COUNT_LIST    
      int localNumP2P = Ci->numP2P;
#else 
      int localNumP2P = 0;
#endif   
    updateInteractionList(sendCells[iparent].ICHILD, localNumP2P,rank);
    if (bodyReferenceCount[key] == 0) {                                          
      aggregateP2PPairs.push_back(crequest_pair (Ci,Cj));        
      bodySendRequestBuffer[aggregateP2PCount*2]     = Cj->IBODY;
      bodySendRequestBuffer[aggregateP2PCount*2 + 1] = Cj->NBODY;          
      aggregateP2PCount++;
      if(aggregateP2PCount < granularity) {
        bodyReferenceCount[key] = -1;        
      } else {  
        aggregateP2PComm(rank);
        if(execute)
          flushP2PInteractions();
        else aggregateP2PCount = 0;
      }
    } else if(bodyReferenceCount[key] < 0) {
      aggregateP2PPairs.push_back(crequest_pair (Ci,Cj)); 
    } else {
      Cj->BODY = recvBodies.begin() + Cj->IBODY;
      if(execute)
        kernel::P2P(Ci,Cj,false); 
      else 
        aggregateP2PPairs.push_back(crequest_pair (Ci,Cj));        
    }
  }

  void P2P(C_iter Ci, C_iter Cj, int rank) {
#if EXAFMM_OVERLAP_REMOTE
    requestLeafCell(Ci,Cj);
#else 
    remoteP2P(Ci,Cj,rank);
#endif    
  }

  void flushP2PInteractions() {     
    aggregateP2PCount = 0;
    B_iter B0 = recvBodies.begin();
    int p2pCount = aggregateP2PPairs.size();
    for (int i = 0; i < p2pCount; ++i) {
      C_iter ci = aggregateP2PPairs[i].first;
      C_iter cj = aggregateP2PPairs[i].second;
      cj->BODY = B0 + cj->IBODY;
      kernel::P2P(ci,cj,false);            
    }    
    aggregateP2PPairs.clear();
  }

  inline void recvBodyNCell(MPI_Status* status, int icell, int nchild){
    int receivedTag  = status->MPI_TAG;
    int recvRank     = status->MPI_SOURCE;
    int responseType = getMessageType(receivedTag);
    int recvCount = 0;            
    if (responseType == bodytag) {
      //startCommTimer("Comm LET bodies");
      MPI_Get_count(status, MPI_INT, &recvCount);
      int bodyCount = recvCount / bodyWordSize;
      Bodies recvBuff(bodyCount);                    
      MPI_Recv((int*)&recvBuff[0], recvCount, MPI_INT, recvRank, receivedTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);        
      //stopCommTimer("Comm LET bodies");
      int ibody, nbody, recvIndex = 0;
      for (int i = 0; i < granularity; ++i) {
        ibody = bodySendRequestBuffer[i*2];
        nbody = bodySendRequestBuffer[i*2 + 1];
        for (int j = 0; j < nbody; ++j) {
          assert(ibody+j < recvBodies.size());
          assert(recvIndex < recvBuff.size());
          recvBodies[ibody+j] = recvBuff[recvIndex];            
          recvIndex++;
        }
        bodySendRequestBuffer[i*2]     = -1;
        bodySendRequestBuffer[i*2 + 1] = -1;          
        if(ibody>=0) bodyReferenceCount[ibody] = 1;                  
      }       
    } else if (responseType == childcelltag || responseType == celltag) {
      Multipoles recvData;    
      //startCommTimer("Comm LET cells");
      MPI_Get_count(status, MPI_INT, &recvCount);
      int cellCount = recvCount / cellWordSize;
      recvData.resize(cellCount);
      MPI_Recv(&recvData[0], recvCount, MPI_INT, recvRank, receivedTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      //stopCommTimer("Comm LET cells");
      if (responseType == celltag) {          
        recvCells[icell] = recvData[0].cell();
        cellReferenceCount[icell]++;
        if(granularity > 1) {
          int index = 1;
          appendLET(recvData, index, 0);
        }
      } else if (responseType == childcelltag)  {
        assert(icell + nchild <= recvCells.size());
        for(int i = 0; i < nchild;++i){
          recvCells[icell+i] = recvData[i].cell();
        }
        //std::copy(recvData.begin(),recvData.begin() + nchild, recvCells.begin()+ icell);          
        childCellsReferenceCount[icell]++;
#if EXAFMM_OVERLAP_REMOTE
        fullfilledCellRequests[icell] = true;
        childCellsReferenceCount[icell] = -1;
#endif        
        if (granularity > 1) {
          int index = nchild;
          for (int i = 0; i < nchild; ++i) appendLET(recvData, index, i);            
        }          
      }
    }     
    else {
      std::cout << "bodies or cells not found for requested tag "<< std::endl;
      char null;
      MPI_Recv(&null, 1, MPI_CHAR, recvRank, receivedTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);      
    }      
  }

  void aggregateP2PComm(int rank) {
    if(aggregateP2PCount > 0) {
      countKernel(numP2PMsgs);              
      int requestType = bodytag;        
      int tag = encryptMessage(requestType, 1, sendbit);
      MPI_Request request;
      MPI_Isend(bodySendRequestBuffer, 2*granularity, MPI_INT, rank, tag, MPI_COMM_WORLD, &request);
      MPI_Status status;     
      toggleDirection(tag);      
      int recv = 0;
      while (!recv) {
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);        
        if (tag != status.MPI_TAG || rank != status.MPI_SOURCE) {
          processIncomingMessage(status.MPI_TAG, status.MPI_SOURCE);          
        } else { recv = 1;}
      }      
      recvBodyNCell(&status, 0,0);
    } 
  }

  struct AggregateTraverse {
  private:
      TreeMPI* treeMPI;   
      std::vector<crequest_pair> interactions;      
      const int rank;
      const int interaction_size;
  public:
    AggregateTraverse(TreeMPI* _traversal, std::vector<crequest_pair> _interactions, int _rank):
    treeMPI(_traversal), interactions(_interactions), rank(_rank), interaction_size(_interactions.size()) { }
    void operator() () { 
      for (int i = 0; i < interaction_size; ++i) {  
        crequest_pair pair = interactions[i];         
        treeMPI->splitCellRemote(pair.first,pair.second,false, 1,rank);
      }
    } 
  };

  struct AggregateP2P {
  private:
      TreeMPI * treeMPI;   
      std::vector<crequest_pair> interactions;      
      const int interaction_size;
  public:
    AggregateP2P(TreeMPI* _traversal, std::vector<crequest_pair> _interactions):
    treeMPI(_traversal), interactions(_interactions),interaction_size(_interactions.size()) { }
    void operator() () { 
      B_iter B0 = treeMPI->recvBodies.begin();
      for (int i = 0; i < interaction_size; ++i) {  
        crequest_pair pair = interactions[i]; 
        C_iter ci = pair.first;
        C_iter cj = pair.second;
        cj->BODY = B0 + cj->IBODY;
        kernel::P2P(ci,cj,false);          
      }
    }
  };

#if EXAFMM_OVERLAP_REMOTE  
  inline void read_write_lock(int const& threadId) {    
    pthread_mutex_lock(mutex+threadId); 
  }
  inline void read_write_unlock(int const& threadId) {
    pthread_mutex_unlock(mutex+threadId); 
  }

  void commBodiesCells(int irank) {    
    MPI_Status status;    
    int ready;
    int rsize;
    int csize ;
    char flush ;    
    bool isFlushSent = false;
    std::vector<crequest_pair> cell_pairs;    
    std::vector<crequest_pair> leaf_cell_pairs; 
    bool bufferEmpty = false;           
    while (!finalizeCommunicationThread || (finalizeCommunicationThread && !bufferEmpty)) {            
      MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &ready, &status);
      if(ready) processIncomingMessage(status.MPI_TAG, status.MPI_SOURCE);             
      for (int i = 0; i < threadCount; ++i) {
        read_write_lock(i);
        csize = cellRequests[i].size(); 
        if(csize > 0) {
          cell_pairs = cellRequests[i];          
          cellRequests[i].resize(0);          
        } 
        read_write_unlock(i);
        for(int j = 0; j< csize; ++j) { 
          crequest_pair pair = cell_pairs[j];
          int const& ichild = pair.second->ICHILD;
          int const& nchild = pair.second->NCHILD;
          if(!fullfilledCellRequests[ichild])            
              sendRecvChildCells(ichild,nchild,0,irank,childcelltag);
          else 
              childCellsReferenceCount[ichild] = -1;
        }
        if(csize > 0) {
          aggregateCellPairs.insert(aggregateCellPairs.end(),cell_pairs.begin(),cell_pairs.end());
        }
        read_write_lock(i);
        rsize = bodyRequests[i].size(); 
        if(rsize > 0) {
          leaf_cell_pairs = bodyRequests[i];
          bodyRequests[i].resize(0);          
        } 
        read_write_unlock(i); 
        for(int j = 0; j< rsize; ++j) 
          remoteP2P(leaf_cell_pairs[j].first,leaf_cell_pairs[j].second,irank,false);                                 
      }
      if(finalizeCommunicationThread) {
        bufferEmpty = true;
        for (int i = 0; i < threadCount; ++i) {
          if(bodyRequests[i].size() > 0 || cellRequests[i].size() > 0) {   
            bufferEmpty = false;
            break;
          }
        }        
        if(bufferEmpty) {
          int pendingCellSize = aggregateCellPairs.size();
          int pendingP2PCellSize = aggregateP2PPairs.size();
          if(pendingCellSize > 0) {
            AggregateTraverse aggregate(this,aggregateCellPairs,irank);
            if(pendingCellSize >= nspawn) {
              wait_omp_tasks;
              create_omp_task(aggregate);           
            } else aggregate();
            aggregateCellPairs.clear();            
            bufferEmpty = false;
          } else if(pendingP2PCellSize > 0) {
            aggregateP2PComm(irank);
            AggregateP2P aggregateP2P(this,aggregateP2PPairs);        
            if(pendingP2PCellSize >=nspawn) {
              wait_omp_tasks;              
              create_omp_task(aggregateP2P);              
            } else aggregateP2P();
            aggregateP2PPairs.clear();             
            bufferEmpty = false;
          } else {
            bufferEmpty = true;            
          }                      
        } 
      }    
    }
  }

  bool requestChildCell(C_iter Ci, C_iter Cj) {
    int const& icell = Cj->ICHILD;
    int const& threadId = get_thread_num();
    if (childCellsReferenceCount[icell] < 0){
      return true;
    } else if (childCellsReferenceCount[icell] == 0) {
      childCellsReferenceCount[icell] = Cj->NCHILD;
    }
    read_write_lock(threadId);
    cellRequests[threadId].push_back(crequest_pair(Ci,Cj));
    read_write_unlock(threadId);
    return false;
  }

  void requestLeafCell(C_iter Ci, C_iter Cj) {
    int const& threadId = get_thread_num();
    read_write_lock(threadId);
    bodyRequests[threadId].push_back(crequest_pair(Ci,Cj));
    read_write_unlock(threadId);
  }

#endif

  void sendRecvChildCells(int key, int nchild, int level, int rank, int requestType) {
    MPI_Request request;
    int tag = encryptMessage(requestType, level, sendbit);      
    MPI_Isend(&key, 1, MPI_INT, rank, tag, MPI_COMM_WORLD, &request);
    MPI_Status status;
    int recvRank = rank;
    countKernel(numM2LMsgs);
    int receivedTag;
    toggleDirection(tag);      
    int recv = 0;
    while (!recv) {
      MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);        
      if (tag != status.MPI_TAG || rank != status.MPI_SOURCE) {
        processIncomingMessage(status.MPI_TAG, status.MPI_SOURCE);          
      } else { recv = 1;}
    }
    receivedTag = status.MPI_TAG;
    recvBodyNCell(&status, key,nchild);
  }

  C_iter getRootCell(int const& rank) {
    recvCells[0] = recvRoots[rank][0].cell();
    if(granularity > 1) {
      int index = 1;
      appendLET(recvRoots[rank],index,0);
    }
    return recvCells.begin();
  }

  //! Get cells underlying a key and a level
  C_iter getChildCells(C_iter Ci, C_iter Cj, int rank, bool& ready) {        
    int key = Cj->ICHILD;
#if EXAFMM_OVERLAP_REMOTE      
    ready = requestChildCell(Ci, Cj);        
#else
    if (childCellsReferenceCount[key] == 0) {            
      int nchild = Cj->NCHILD;
      int level = 0;//Cj->LEVEL;
      sendRecvChildCells(key,nchild,level,rank,childcelltag);
      ready = true;
    } else {
      ready = true;
    }
#endif
    return (recvCells.begin() + key);      
  }  

public:
  //! Constructor
  TreeMPI(int _mpirank, int _mpisize, int _images) :
    mpirank(_mpirank), mpisize(_mpisize), images(_images),
    granularity(1), terminated(0), hitCount(0), nulltag(1), celltag(2),
    childcelltag(3), bodytag(8),flushtag(9), maxtag(15), levelshift(5),
    requestshift(4), directionshift(1), levelmask(0x1F),
    requestmask(0xF), directionmask(0x1), grainmask(0xFFFF),sendbit(1), 
    receivebit(0),nspawn(5000)    
#if EXAFMM_COUNT_KERNEL
    ,numP2P(0),numM2L(0),numP2PMsgs(0),numM2LMsgs(0)
#endif
     {                            // Initialize variables
    allBoundsXmin = new float [mpisize][3];                   // Allocate array for minimum of local domains
    allBoundsXmax = new float [mpisize][3];                   // Allocate array for maximum of local domains
    sendBodyCount = new int [mpisize];                        // Allocate send count
    sendBodyDispl = new int [mpisize];                        // Allocate send displacement
    recvBodyCount = new int [mpisize];                        // Allocate receive count
    recvBodyDispl = new int [mpisize];                        // Allocate receive displacement
    sendCellCount = new int [mpisize];                        // Allocate send count
    sendCellDispl = new int [mpisize];                        // Allocate send displacement
    recvCellCount = new int [mpisize];                        // Allocate receive count
    recvCellDispl = new int [mpisize];                        // Allocate receive displacement
    recvRoots = new Multipoles[mpisize];                      // Allocate receive buffer for tree roots
    bodySendRequestBuffer = new int [granularity*2];          // Allocate aggregate bodies send buffer
    bodyReceiveRequestBuffer = new int [granularity*2];       // Allocate aggregate bodies recv buffer
#if EXAFMM_OVERLAP_REMOTE
    threadCount = get_num_threads();
    bodyRequests = new std::vector<crequest_pair >[threadCount];
    cellRequests = new std::vector<crequest_pair >[threadCount];    
    finalizeCommunicationThread = false;  
    mutex = new pthread_mutex_t[threadCount];
    for (int i = 0; i < threadCount; ++i) {
      pthread_mutex_init(mutex+i,NULL);   
    }
#endif    
  }
  //! Temporary constructor added to prevent breaking compatibility with other kernels using this class
  TreeMPI(int _mpirank, int _mpisize, int _images, int _grainSize, int _threads, int _nspawn) :
    mpirank(_mpirank), mpisize(_mpisize), images(_images),
    granularity(_grainSize), terminated(0), hitCount(0), nulltag(1), celltag(2),
    childcelltag(3), bodytag(8),flushtag(9), maxtag(15), levelshift(5),
    requestshift(4), directionshift(1), levelmask(0x1F),
    requestmask(0xF), directionmask(0x1), grainmask(0xFFFF),sendbit(1), 
    receivebit(0),nspawn(_nspawn)
#if EXAFMM_COUNT_KERNEL
    ,numP2P(0),numM2L(0),numP2PMsgs(0),numM2LMsgs(0)
#endif
     {                            // Initialize variables
    allBoundsXmin = new float [mpisize][3];                   // Allocate array for minimum of local domains
    allBoundsXmax = new float [mpisize][3];                   // Allocate array for maximum of local domains
    sendBodyCount = new int [mpisize];                        // Allocate send count
    sendBodyDispl = new int [mpisize];                        // Allocate send displacement
    recvBodyCount = new int [mpisize];                        // Allocate receive count
    recvBodyDispl = new int [mpisize];                        // Allocate receive displacement
    sendCellCount = new int [mpisize];                        // Allocate send count
    sendCellDispl = new int [mpisize];                        // Allocate send displacement
    recvCellCount = new int [mpisize];                        // Allocate receive count
    recvCellDispl = new int [mpisize];                        // Allocate receive displacement
    recvRoots = new Multipoles[mpisize];                      // Allocate receive buffer for tree roots
    bodySendRequestBuffer = new int [granularity*2];          // Allocate aggregate bodies send buffer
    bodyReceiveRequestBuffer = new int [granularity*2];       // Allocate aggregate bodies recv buffer
#if EXAFMM_OVERLAP_REMOTE
    set_num_threads(_threads);
    threadCount = _threads;
    bodyRequests = new std::vector<crequest_pair>[threadCount];
    cellRequests = new std::vector<crequest_pair>[threadCount];    
    finalizeCommunicationThread = false;  
    mutex = new pthread_mutex_t[threadCount];
    for (int i = 0; i < threadCount; ++i) {
      pthread_mutex_init(mutex+i,NULL);   
    }
#endif      
  }
  //! Destructor
  ~TreeMPI() {
    delete[] allBoundsXmin;                                   // Deallocate array for minimum of local domains
    delete[] allBoundsXmax;                                   // Deallocate array for maximum of local domains
    delete[] sendBodyCount;                                   // Deallocate send count
    delete[] sendBodyDispl;                                   // Deallocate send displacement
    delete[] recvBodyCount;                                   // Deallocate receive count
    delete[] recvBodyDispl;                                   // Deallocate receive displacement
    delete[] sendCellCount;                                   // Deallocate send count
    delete[] sendCellDispl;                                   // Deallocate send displacement
    delete[] recvCellCount;                                   // Deallocate receive count
    delete[] recvCellDispl;                                   // Deallocate receive displacement
    delete[] bodySendRequestBuffer;                           // Deallocate aggregate bodies send buffer
    delete[] bodyReceiveRequestBuffer;                        // Deallocate aggregate bodies recv buffer
    delete[] recvRoots;                                       // Deallocate receive buffer for tree roots
#if EXAFMM_OVERLAP_REMOTE
    delete[] bodyRequests;
    delete[] cellRequests;
    delete[] mutex;
#endif        
  }

  //! Allgather bounds from all ranks
  void allgatherBounds(Bounds bounds) {
    float Xmin[3], Xmax[3];
    for (int d = 0; d < 3; d++) {                             // Loop over dimensions
      Xmin[d] = bounds.Xmin[d];                               //  Convert Xmin to float
      Xmax[d] = bounds.Xmax[d];                               //  Convert Xmax to float
    }                                                         // End loop over dimensions
    MPI_Allgather(Xmin, 3, MPI_FLOAT, allBoundsXmin[0], 3, MPI_FLOAT, MPI_COMM_WORLD);// Gather all domain bounds
    MPI_Allgather(Xmax, 3, MPI_FLOAT, allBoundsXmax[0], 3, MPI_FLOAT, MPI_COMM_WORLD);// Gather all domain bounds
  }

  //! Set local essential tree to send to each process
  void setLET(Cells & cells, vec3 cycle) {
    logger::startTimer("Set LET size");                       // Start timer
    C_iter C0 = cells.begin();                                // Set cells begin iterator
    Bounds bounds;                                            // Bounds of local subdomain
    sendBodyDispl[0] = 0;                                     // Initialize body displacement vector
    sendCellDispl[0] = 0;                                     // Initialize cell displacement vector
    for (int irank = 0; irank < mpisize; irank++) {           // Loop over ranks
      if (irank != 0) sendBodyDispl[irank] = sendBodyDispl[irank - 1] + sendBodyCount[irank - 1]; // Update body displacement
      if (irank != 0) sendCellDispl[irank] = sendCellDispl[irank - 1] + sendCellCount[irank - 1]; // Update cell displacement
      sendBodyCount[irank] = 0;                               //  Initialize send body count for current rank
      sendCellCount[irank] = 0;                               //  Initialize send cell count for current rank
      if (irank != mpirank && !cells.empty()) {               //  If not current rank and cell vector is not empty
        int ibody = 0;                                        //   Initialize send body's offset
        int icell = 1;                                        //   Initialize send cell's offset
        for (int d = 0; d < 3; d++) {                         //   Loop over dimensions
          bounds.Xmin[d] = allBoundsXmin[irank][d];           //    Local Xmin for irank
          bounds.Xmax[d] = allBoundsXmax[irank][d];           //    Local Xmax for irank
        }                                                     //   End loop over dimensions
        if (C0->NCHILD == 0) {                                //   If root cell is leaf
          addSendBody(C0, irank, ibody, icell - 1, false);    //    Add bodies to send
        }                                                     //   End if for root cell leaf
        traverseLET(C0, C0, bounds, cycle, irank, ibody, icell, 0, false); // Traverse tree to set LET
        sendBodyCount[irank] = ibody;                         //   Send body count for current rank
        sendCellCount[irank] = icell;                         //   Send cell count for current rank
      }                                                       //  Endif for current rank
    }                                                         // End loop over ranks
    logger::stopTimer("Set LET size");                        // Stop timer
    logger::startTimer("Set LET");                            // Start timer
    int numSendBodies = sendBodyDispl[mpisize - 1] + sendBodyCount[mpisize - 1]; // Total number of send bodies
    int numSendCells = sendCellDispl[mpisize - 1] + sendCellCount[mpisize - 1]; // Total number of send cells
    sendBodies.resize(numSendBodies);                         // Clear send buffer for bodies
    sendCells.resize(numSendCells);                           // Clear send buffer for cells
    MPI_Barrier(MPI_COMM_WORLD);
    for (int irank = 0; irank < mpisize; irank++) {           // Loop over ranks
      if (irank != mpirank && !cells.empty()) {               //  If not current rank and cell vector is not empty
        int ibody = 0;                                        //   Reinitialize send body's offset
        int icell = 0;                                        //   Reinitialize send cell's offset
        for (int d = 0; d < 3; d++) {                         //   Loop over dimensions
          bounds.Xmin[d] = allBoundsXmin[irank][d];           //   Local Xmin for irank
          bounds.Xmax[d] = allBoundsXmax[irank][d];           //   Local Xmax for irank
        }                                                     //   End loop over dimensions
        C_iter Csend = sendCells.begin() + sendCellDispl[irank];//   Send cell iterator
        *Csend = *C0;                                         //   Copy cell to send buffer
        Csend->NCHILD = Csend->NBODY = 0;                     //   Reset link to children and bodies
        icell++;                                              //   Increment send cell counter
        if (C0->NCHILD == 0) {                                //   If root cell is leaf
          addSendBody(C0, irank, ibody, icell - 1, true);     //    Add bodies to send
        }                                                     //   End if for root cell leaf
        traverseLET(C0, C0, bounds, cycle, irank, ibody, icell, 0, true); // Traverse tree to set LET
      }                                                       //  Endif for current rank
    }                                                         // End loop over ranks
    logger::stopTimer("Set LET");                             // Stop timer
  }

  //! Get local essential tree from irank
  void getLET(Cells & cells, int irank) {
    std::stringstream event;                                  // Event name
    event << "Get LET from rank " << irank;                   // Create event name based on irank
    logger::startTimer(event.str());                          // Start timer
    for (int i = 0; i < recvCellCount[irank]; i++) {          // Loop over receive cells
      C_iter C = recvCells.begin() + recvCellDispl[irank] + i;//  Iterator of receive cell
      if (C->NBODY != 0) {                                    //  If cell has bodies
        C->BODY = recvBodies.begin() + recvBodyDispl[irank] + C->IBODY;// Iterator of first body
      }                                                       //  End if for bodies
    }                                                         // End loop over receive cells
    cells.resize(recvCellCount[irank]);                       // Resize cell vector for LET
    cells.assign(recvCells.begin() + recvCellDispl[irank],    // Assign receive cells to vector
                 recvCells.begin() + recvCellDispl[irank] + recvCellCount[irank]);
    logger::stopTimer(event.str());                           // Stop timer
  }

  //! Link LET with received bodies and calcualte NBODY
  void linkLET() {
    logger::startTimer("Link LET");                           // Start timer
    for (int irank = 0; irank < mpisize; irank++) {           // Loop over ranks
      for (int i = 0; i < recvCellCount[irank]; i++) {        //  Loop over receive cells
        C_iter C = recvCells.begin() + recvCellDispl[irank] + i;//   Iterator of receive cell
        if (C->NBODY != 0) {                                  //   If cell has bodies
          C->BODY = recvBodies.begin() + recvBodyDispl[irank] + C->IBODY;// Iterator of first body
        }                                                     //   End if for bodies
      }                                                       //  End loop over receive cells
    }                                                         // End loop over ranks
    logger::stopTimer("Link LET");                            // End timer
  }

  //! Sets local essential tree size
  void setLETSize(Cells & cells, vec3 cycle) {
    C_iter C0 = cells.begin();                                // Set cells begin iterator
    Bounds bounds;                                            // Bounds of local subdomain
    sendBodyDispl[0] = 0;                                     // Initialize body displacement vector
    sendCellDispl[0] = 0;                                     // Initialize cell displacement vector
    for (int irank = 0; irank < mpisize; irank++) {           // Loop over ranks
      if (irank != 0) sendBodyDispl[irank] = sendBodyDispl[irank - 1] + sendBodyCount[irank - 1]; // Update body displacement
      if (irank != 0) sendCellDispl[irank] = sendCellDispl[irank - 1] + sendCellCount[irank - 1]; // Update cell displacement
      sendBodyCount[irank] = 0;                               //  Initialize send body count for current rank
      sendCellCount[irank] = 0;                               //  Initialize send cell count for current rank
      if (irank != mpirank && !cells.empty()) {               //  If not current rank and cell vector is not empty
        int ibody = 0;                                        //   Initialize send body's offset
        int icell = 1;                                        //   Initialize send cell's offset
        for (int d = 0; d < 3; d++) {                         //   Loop over dimensions
          bounds.Xmin[d] = allBoundsXmin[irank][d];           //    Local Xmin for irank
          bounds.Xmax[d] = allBoundsXmax[irank][d];           //    Local Xmax for irank
        }                                                     //   End loop over dimensions
        if (C0->NCHILD == 0) {                                //   If root cell is leaf
          addSendBody(C0, irank, ibody, icell - 1, false);    //    Add bodies to send
        }                                                     //   End if for root cell leaf
        traverseLET(C0, C0, bounds, cycle, irank, ibody, icell, 0, false); // Traverse tree to set LET
        sendBodyCount[irank] = ibody;                         //   Send body count for current rank
        sendCellCount[irank] = icell;                         //   Send cell count for current rank
      }                                                       //  Endif for current rank
    }                                                         // End loop over ranks
  }

  //! Copy remote root cells to body structs (for building global tree)
  Bodies root2body() {
    logger::startTimer("Root to body");                       // Start timer
    Bodies bodies;                                            // Bodies to contain remote root coordinates
    bodies.reserve(mpisize - 1);                              // Reserve size of body vector
    for (int irank = 0; irank < mpisize; irank++) {           // Loop over ranks
      if (irank != mpirank) {                                 //  If not current rank
        C_iter C0 = recvCells.begin() + recvCellDispl[irank]; //   Root cell iterator for irank
        Body body;                                            //   Body to contain remote root coordinates
        body.X = C0->X;                                       //   Copy remote root coordinates
        body.IBODY = recvCellDispl[irank];                    //   Copy remote root displacement in vector
        bodies.push_back(body);                               //   Push this root cell to body vector
      }                                                       //  End if for not current rank
    }                                                         // End loop over ranks
    logger::stopTimer("Root to body");                        // Stop timer
    return bodies;                                            // Return body vector
  }

  //! Graft remote trees to global tree
  void attachRoot(Cells & cells) {
    logger::startTimer("Attach root");                        // Start timer
    int globalCells = cells.size();                           // Number of global cells
    cells.insert(cells.end(), recvCells.begin(), recvCells.end()); // Join LET cell vectors
    for (C_iter C = cells.begin(); C != cells.begin() + globalCells; C++) { // Loop over global cells
      if (C->NCHILD == 0) {                                   // If leaf cell
        int offset = globalCells + C->BODY->IBODY;            //  Offset of received root cell index
        C_iter C0 = cells.begin() + offset;                   //  Root cell iterator
        C0->IPARENT = C->IPARENT;                             //  Link remote root to global leaf
        *C = *C0;                                             //  Copy remote root to global leaf
        C->ICHILD += offset;                                  //  Add offset to child index
      } else {                                                // If not leaf cell
        C->BODY = recvBodies.end();                           //  Use BODY as flag to indicate non-leaf global cell
      }                                                       // End if for leaf cell
    }                                                         // End loop over global cells
    for (int irank = 0; irank < mpisize; irank++) {           // Loop over ranks
      if (irank != mpirank && recvCellCount[irank] > 0) {     //  If not current rank
        C_iter C0 = cells.begin() + globalCells + recvCellDispl[irank];// Root cell iterator for irank
        for (C_iter C = C0 + 1; C != C0 + recvCellCount[irank]; C++) { //    Loop over cells received from irank
          int offset = globalCells + recvCellDispl[irank];    //     Offset of received root cell index
          C->IPARENT += offset;                               //     Add offset to parent index
          C->ICHILD += offset;                                //     Add offset to child index
        }                                                     //    End loop over cells received from irank
      }                                                       //  End if for not current rank
    }                                                         // End loop over ranks
    C_iter C0 = cells.begin();                                // Root cell iterator of entire LET
    for (int i = globalCells - 1; i >= 0; i--) {              // Loop over global cells bottom up
      C_iter C = cells.begin() + i;                           //  Iterator of current cell
      if (C->BODY == recvBodies.end()) {                      //  If non-leaf global cell
        vec3 Xmin = C->X, Xmax = C->X;                        //   Initialize Xmin, Xmax
        for (C_iter CC = C0 + C->ICHILD; CC != C0 + C->ICHILD + C->NCHILD; CC++) { // Loop over child cells
          Xmin = min(CC->X - CC->R, Xmin);                    //    Update Xmin
          Xmax = max(CC->X + CC->R, Xmax);                    //    Update Xmax
        }                                                     //   End loop over child cells
        C->X = (Xmax + Xmin) / 2;                             //   Calculate center of domain
        for (int d = 0; d < 3; d++) {                         //   Loop over dimensions
          C->R = std::max(C->X[d] - Xmin[d], C->R);           //    Calculate min distance from center
          C->R = std::max(Xmax[d] - C->X[d], C->R);           //    Calculate max distance from center
        }                                                     //   End loop over dimensions
        C->M = 0;                                             //   Reset multipoles
        kernel::M2M(C, C0);                                   //   M2M kernel
      }                                                       //  End if for non-leaf global cell
    }                                                         // End loop over global cells bottom up
    logger::stopTimer("Attach root");                         // Stop timer
  }

  //! Send bodies
  Bodies commBodies(Bodies bodies) {
    logger::startTimer("Comm partition");                     // Start timer
    alltoall(bodies);                                         // Send body count    
#if EXAFMM_USE_ALLTOALL    
    alltoallv(bodies);                                        // Send bodies        
#else    
    alltoallv_p2p(bodies,recvBodies,recvBodyDispl,recvBodyCount,sendBodyDispl,sendBodyCount);
#endif    
    logger::stopTimer("Comm partition");                      // Stop timer
    return recvBodies;                                        // Return received bodies
  }

  //! Send bodies
  Bodies commBodies() {
    logger::startTimer("Comm LET bodies");                    // Start timer
    alltoall(sendBodies);                                     // Send body count
#if EXAFMM_USE_ALLTOALL        
    alltoallv(sendBodies);                                    // Send bodies
#else        
    alltoallv_p2p(sendBodies,recvBodies,recvBodyDispl,recvBodyCount,sendBodyDispl,sendBodyCount);
#endif        
    logger::stopTimer("Comm LET bodies");                     // Stop timer
    return recvBodies;                                        // Return received bodies
  }

  //! Send cells
  void commCells() {
    logger::startTimer("Comm LET cells");                     // Start timer
    alltoall(sendCells);                                      // Send cell count
#if EXAFMM_USE_ALLTOALL            
    alltoallv(sendCells);                                     // Senc cells
#else
    alltoallv_p2p(sendCells,recvCells,recvCellDispl,recvCellCount,sendCellDispl,sendCellCount);
#endif    
    logger::stopTimer("Comm LET cells");                      // Stop timer
  }

  //! Copy recvBodies to bodies
  Bodies getRecvBodies() {
    return recvBodies;                                        // Return recvBodies
  }

#if EXAFMM_COUNT_LIST
  inter_list& getRemoteInteractionList() {
    return remoteInteractionList;
  }
#endif

#if EXAFMM_COUNT_KERNEL
  real_t& getRemoteP2PCount() {
    return numP2P;
  }

  real_t& getRemoteM2LCount() {
    return numM2L;
  }
#endif


  //! Send bodies to next rank (round robin)
  void shiftBodies(Bodies & bodies) {
    int newSize;                                              // New number of bodies
    int oldSize = bodies.size();                              // Current number of bodies
    const int word = sizeof(bodies[0]) / 4;                   // Word size of body structure
    const int isend = (mpirank + 1          ) % mpisize;      // Send to next rank (wrap around)
    const int irecv = (mpirank - 1 + mpisize) % mpisize;      // Receive from previous rank (wrap around)
    MPI_Request sreq, rreq;                                   // Send, receive request handles

    MPI_Isend(&oldSize, 1, MPI_INT, irecv, 0, MPI_COMM_WORLD, &sreq);// Send current number of bodies
    MPI_Irecv(&newSize, 1, MPI_INT, isend, 0, MPI_COMM_WORLD, &rreq);// Receive new number of bodies
    MPI_Wait(&sreq, MPI_STATUS_IGNORE);                       // Wait for send to complete
    MPI_Wait(&rreq, MPI_STATUS_IGNORE);                       // Wait for receive to complete

    recvBodies.resize(newSize);                               // Resize buffer to new number of bodies
    MPI_Isend((int*)&bodies[0], oldSize * word, MPI_INT, irecv, // Send bodies to next rank
              1, MPI_COMM_WORLD, &sreq);
    MPI_Irecv((int*)&recvBodies[0], newSize * word, MPI_INT, isend, // Receive bodies from previous rank
              1, MPI_COMM_WORLD, &rreq);
    MPI_Wait(&sreq, MPI_STATUS_IGNORE);                       // Wait for send to complete
    MPI_Wait(&rreq, MPI_STATUS_IGNORE);                       // Wait for receive to complete
    bodies = recvBodies;                                      // Copy bodies from buffer
  }

  //! Allgather bodies
  Bodies allgatherBodies(Bodies & bodies) {
    const int word = sizeof(bodies[0]) / 4;                   // Word size of body structure
    sendBodyCount[0] = bodies.size();                         // Determine send count
    MPI_Allgather(sendBodyCount, 1, MPI_INT,                  // Allgather number of bodies
                  recvBodyCount, 1, MPI_INT, MPI_COMM_WORLD);
    recvBodyDispl[0] = 0;                                     // Initialize receive displacement
    for (int irank = 0; irank < mpisize - 1; irank++) {       // Loop over ranks
      recvBodyDispl[irank + 1] = recvBodyDispl[irank] + recvBodyCount[irank]; // Set receive displacement
    }                                                         // End loop over ranks
    recvBodies.resize(recvBodyDispl[mpisize - 1] + recvBodyCount[mpisize - 1]); // Resize receive buffer
    for (int irank = 0; irank < mpisize; irank++) {           // Loop over ranks
      recvBodyCount[irank] *= word;                           //  Multiply receive count by word size of data
      recvBodyDispl[irank] *= word;                           //  Multiply receive displacement by word size of data
    }                                                         // End loop over ranks
    MPI_Allgatherv((int*)&bodies[0], sendBodyCount[0]*word, MPI_INT,// Allgather bodies
                   (int*)&recvBodies[0], recvBodyCount, recvBodyDispl, MPI_INT, MPI_COMM_WORLD);
    return recvBodies;                                        // Return bodies
  }

  void setSendBodyNCells(Cells const& cells, Bodies const& bodies) {
    sendCells = cells;
    sendBodies = bodies;
  }

  void exchangeTreeRoots() {
    logger::startTimer("Set Tree Roots");
    Multipoles* roots = new Multipoles();
    int previous = 0;
    int index = 0;
    sendCellCount[0] = 0;
    sendCellDispl[0] = 0;
    C_iter C0 = sendCells.begin();
    for(int i = 0; i< mpisize; ++i) {
      if(i!=mpirank) {        
        roots->push_back(Multipole(C0));
        index++;
        if(granularity > 1)
          setLETSubset(roots, C0, sendCells[0], granularity, index, i);
      }      
      sendCellCount[i] = index - previous;
      if(i != 0) sendCellDispl[i] = sendCellDispl[i-1] + sendCellCount[i-1];
      previous = index;
    }
    logger::stopTimer("Set Tree Roots");
    logger::startTimer("Exchange Tree Roots");
    alltoall(sendCells);
    Multipoles recvMults = alltoallv(*roots); 
    M_iter M0 =recvMults.begin();   
    for (int i = 0; i < mpisize; ++i) 
      if(i!=mpirank && recvCellCount[i] > 0)
        recvRoots[i] = Multipoles(M0 + recvCellDispl[i], M0 + recvCellDispl[i] + recvCellCount[i]);          
    logger::stopTimer("Exchange Tree Roots"); 
  }  

  void setSendLET() {
    exchangeTreeRoots();
    logger::startTimer("Set LET");    
    Multipole _multipole;   
    cellWordSize = sizeof(Multipole) / 4;
    bodyWordSize = sizeof(sendBodies[0]) / 4;
    int counts[2] = {sendBodies.size(),sendCells.size()};
    int maxCounts[2];
    MPI_Allreduce(counts, maxCounts, 2, MPI_INT, MPI_MAX, MPI_COMM_WORLD);// Reduce domain Xmin
    int const& maxBodyCount = maxCounts[0];
    recvBodies.resize(maxBodyCount);
    bodyReferenceCount.resize(maxBodyCount);
#if EXAFMM_OVERLAP_REMOTE
    fullfilledBodyRequests.resize(maxBodyCount);
#endif    
    int const& maxCellCount = maxCounts[1];
    aggregateP2PCount = 0;
    terminated = 0;
    remoteCommunicationTime = 0;                                
    recvCells.resize(maxCellCount);
    cellReferenceCount.resize(maxCellCount);              
    childCellsReferenceCount = new int[maxCellCount];
#if EXAFMM_OVERLAP_REMOTE
    fullfilledCellRequests.resize(maxCellCount);
#endif
    remoteInteractionList.resize(maxCellCount);
    remoteInteractionList.assign(maxCellCount,std::vector<int>(mpisize));
    for (int i = 0; i < 2*granularity; ++i) {
      bodySendRequestBuffer[i]=-1;
    }
    logger::stopTimer("Set LET");    
  }

  inline void startCommTimer(const char*  event) {
#if EXAFMM_TIME_COMM
      logger::stopTimer("Traverse Remote",0);                       
      logger::startTimer(event);
#endif
  }

  inline void stopCommTimer(const char* event) {
#if EXAFMM_TIME_COMM
      remoteCommunicationTime += logger::stopTimer(event, 0);
      logger::startTimer("Traverse Remote");      
#endif
  }

  void sendFlushRequest(char* null_request) {
    MPI_Request request;
    int tag = encryptMessage(flushtag, 0, sendbit);    
    for (int i = 0; i < mpisize; ++i)
      if (i != mpirank) 
        MPI_Isend(&null_request, 1, MPI_CHAR, i, tag, MPI_COMM_WORLD, &request);         
  }

  void flushAllRequests() { 
    char flush;
    sendFlushRequest(&flush);
    MPI_Status status;    
    while (terminated < (mpisize - 1)) {      
      MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
      processIncomingMessage(status.MPI_TAG, status.MPI_SOURCE);      
    }          
    deallocateCompletedRequests(pendingCellMPI_Requests,sendCellBuffers);    
    deallocateCompletedRequests(pendingBodyMPI_Requests,sendBodyBuffers);
    MPI_Barrier(MPI_COMM_WORLD);                                //<! needed to finalize flush calls for remote traversal prior to going to the direct loop
  }
   //! Evaluate P2P and M2L using dual tree traversal  
  void dualTreeTraversalRemote(Cells & icells, Bodies& ibodies, int mpirank, int mpisize, real_t remote = 1) {
#if EXAFMM_COUNT_KERNEL
   numP2P = 0;                                                  //!< Reset number of P2P kernel calls
   numM2L = 0;                                                  //!< Reset number of M2L kernel calls
   numP2PMsgs = 0;                                              //!< Reset number of P2P Message requests
   numM2LMsgs = 0;                                              //!< Reset number of M2L Message requests
#endif
    if (icells.empty()) return;                                 // Quit if either of the cell vectors are empty  
    logger::startTimer("Traverse Remote");                      // Start timer
    kernel::Xperiodic = 0;    
    logger::initTracer();                                       // Initialize tracer        
    Ci0 = icells.begin();                                       // Set iterator of target root cell  
#if EXAFMM_OVERLAP_REMOTE
#pragma omp parallel sections
{
  #pragma omp section
  {
    for (int irank = 0; irank < mpisize; ++irank) {
      if (irank != mpirank) {         
        bodyReferenceCount.assign(bodyReferenceCount.size(),0);
        cellReferenceCount.assign(cellReferenceCount.size(),0);
        const int N = cellReferenceCount.size();
        for (int j = 0; j < N; ++j)  childCellsReferenceCount[j] = 0;        
        fullfilledCellRequests.assign(fullfilledCellRequests.size(),false);
        fullfilledBodyRequests.assign(fullfilledBodyRequests.size(),false);        
        for(int i = 0; i < threadCount; ++i) {          
          cellRequests[i].clear();
          bodyRequests[i].clear();          
        }
        finalizeCommunicationThread = false;
        #pragma omp task
        {          
          commBodiesCells(irank);
        }
        #pragma omp task
        {
          Cj0 = getRootCell(irank);
          traverseRemote(Ci0, Cj0, false, remote, irank);
          finalizeCommunicationThread = true;
          #pragma omp flush(finalizeCommunicationThread)
        }
        #pragma omp taskwait
      }
    }
    
  }
}
#else 
     for (int irank = 0; irank < mpisize; ++irank) {
      if (irank != mpirank) {         
        bodyReferenceCount.assign(bodyReferenceCount.size(),0);
        cellReferenceCount.assign(cellReferenceCount.size(),0);
        const int N = cellReferenceCount.size();
        for (int j = 0; j < N; ++j)  childCellsReferenceCount[j] = 0;     
          Cj0 = getRootCell(irank);
          traverseRemote(Ci0, Cj0, false, remote, irank);        
          aggregateP2PComm(irank);
          flushP2PInteractions();
        }
      }
#endif
#if EXAFMM_TIME_COMM
    logger::printTime("Comm LET bodies");      
    logger::printTime("Comm LET cells");      
#endif    
    logger::stopTimer("Traverse Remote");                       // Stop timer
    logger::writeTracer();                                      // Write tracer to file
  }


  void writeRemoteTraversalData(int mpirank){
#if EXAFMM_COUNT_KERNEL    
    std::stringstream name;                                   // File name
    name << "num" << std::setfill('0') << std::setw(6)        // Set format
         << mpirank << ".dat";                                // Create file name for list
    std::ofstream listFile(name.str().c_str(), std::ios::app);// Open list log file
    listFile << std::setw(logger::stringLength) << std::left  //  Set format
      << "Remote P2P calls" << " " << numP2P << std::endl;    //  Print event and timer
    listFile << std::setw(logger::stringLength) << std::left  //  Set format
      << "Remote M2L calls" << " " << numM2L << std::endl;    //  Print event and timer
    listFile << std::setw(logger::stringLength) << std::left  //  Set format
      << "Remote P2P Msgs" << " " << numP2PMsgs << std::endl; //  Print event and timer
    listFile << std::setw(logger::stringLength) << std::left  //  Set format
      << "Remote M2L Msgs" << " " << numM2LMsgs << std::endl; //  Print event and timer
#endif
  }
};
}
#endif
