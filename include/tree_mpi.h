#ifndef tree_mpi_h
#define tree_mpi_h
#include "kernel.h"
#include "logger.h"
#include "sort.h"
#include <set>
#include <queue>

namespace exafmm {
//! Handles all the communication of local essential trees
class TreeMPI {
protected:  
  const int mpirank;                                          //!< Rank of MPI communicator
  const int mpisize;                                          //!< Size of MPI communicator
  const int images;                                           //!< Number of periodic image sublevels
  int granularity;                                            //!< The granularity of communication
  float (* allBoundsXmin)[3];                                 //!< Array for local Xmin for all ranks
  float (* allBoundsXmax)[3];                                 //!< Array for local Xmax for all ranks
  Bodies sendBodies;                                          //!< Send buffer for bodies
  Bodies recvBodies;                                          //!< Receive buffer for bodies
  Cells sendCells;                                            //!< Send buffer for cells
  Cells recvCells;                                            //!< Receive buffer for cells
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
  typedef std::vector<std::vector<int> > NeighborList;        //!< Type of Neighbor list  
  struct GraphNode
  {
  public:
    GraphNode(int p, int number, int src)
    :parent(p),id(number),source(src){ }
    int parent;
    int id;
    int source;
  };
  typedef std::vector<GraphNode> TreeLevel;                  //! The FMM communication tree level 
  typedef std::vector<TreeLevel> CommTree;                   //! The FMM communication tree
  typedef std::vector<CommTree> CommGraph;                   //! The Hierarchical communication tree
  
private:
  //! Exchange send count for bodies
  void alltoall(Bodies & bodies) {
    for (int i = 0; i < mpisize; i++) {                       // Loop over ranks
      sendBodyCount[i] = 0;                                   //  Initialize send counts
    }                                                         // End loop over ranks
    for (B_iter B = bodies.begin(); B != bodies.end(); B++) { // Loop over bodies
      assert(0 <= B->IRANK && B->IRANK < mpisize);            //  Check bounds for process ID
      sendBodyCount[B->IRANK]++;                              //  Fill send count bucket
      B->IRANK = mpirank;                                     //  Tag for sending back to original rank
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
    std::copy(localBuffer, localBuffer + sendCount[mpirank], recvB.begin() + recvDispl[mpirank]);   
  }

  //! Exchange bodies in a point-to-point manner
  template<typename T>
  void alltoallv_p2p_onesided(T& sendB, T& recvB, int* recvDispl, int* recvCount, int* sendDispl, int* sendCount) { 
    //std::vector<int> v = getNeighborRanks();       
    int dataSize = recvDispl[mpisize - 1] + recvCount[mpisize - 1];    
    if (sendB.size() == sendCount[mpirank] && dataSize == sendCount[mpirank]) {
      recvB = sendB;
      return;
    }    
    MPI_Win win;          
    int* remoteDispl = new int[mpisize];    
    MPI_Alltoall(sendDispl, 1, MPI_INT, remoteDispl, 1, MPI_INT, MPI_COMM_WORLD);                
    recvB.resize(dataSize);
    assert( (sizeof(sendB[0]) & 3) == 0 );
    int word = sizeof(sendB[0]) / 4;    
    MPI_Win_create((int*)&sendB[0], sendB.size()*word*sizeof(int), sizeof(int), MPI_INFO_NULL, MPI_COMM_WORLD, &win);          
    for (int irank = 0; irank < mpisize; ++irank) {
      int sendRank = (mpirank + 1 + irank) % mpisize;
      if (sendRank != mpirank) {
        MPI_Win_lock(MPI_LOCK_SHARED, sendRank, 0, win);
        MPI_Get((int*)&recvB[0] + recvDispl[sendRank]*word, recvCount[sendRank]*word,MPI_INT, sendRank, remoteDispl[sendRank]*word,recvCount[sendRank]*word,MPI_INT,win);
        MPI_Win_unlock(sendRank, win);        
      }                     
    }             
    MPI_Win_free(&win);
    delete[] remoteDispl;        
    typename T::iterator localBuffer = sendB.begin() + sendDispl[mpirank];
    std::copy(localBuffer, localBuffer + sendCount[mpirank], recvB.begin() + recvDispl[mpirank]);       
  }

  template<typename PathT, typename Msg>
  void getMessageInfo(PathT const& path, Msg& content, int rank , int stage, int const& logP) {
    if(stage < logP) {
      for (int l = stage; l < logP; ++l) {
        int newDest = path[rank][l];
        content.push_back(newDest);
        getMessageInfo(path, content, newDest ,l+1, logP);
      }
    }
  }

  template<typename PathT>
  bool searchMessagePath(PathT const& path, int rank , int stage, int const& logP, int const& dest) {    
    int newDest = path[rank][stage];                
    if(stage >= logP)
      return false;
    else if(newDest == dest)
      return true;
    else if(!searchMessagePath(path, rank , stage + 1, logP, dest))
      return searchMessagePath(path, newDest , stage + 1, logP, dest);
    else 
      return true;
  }

  std::vector<std::vector<int> > getHypercubeMatrix(int P) {
    int logP = log2(P);
    std::vector<std::vector<int> >path(P,std::vector<int> (logP));
    for (int irank = 0; irank < P; ++irank) {
      for (int i = 0; i < logP; ++i) 
        path[irank][i] = irank ^ (1<<i); 
    }
    return path;
  }

  int getTag(int source, int destination, int const& shift) {
    int tag = source;
    tag<<=shift;
    tag|=destination;
    return tag;
  }

  void getOriginNDestination(int tag, int const& mask, int const& shift, int& origin, int& dest) {
    dest = tag&mask;
    origin = tag>>shift;
  }

  //! Exchange bodies in a point-to-point manner following hypercube pattern    
  template<typename T>  
  void alltoallv_p2p_hypercube(T& sendB, T& recvB, int* recvDispl, int* recvCount, int* sendDispl, int* sendCount) {        
    int dataSize = recvDispl[mpisize - 1] + recvCount[mpisize - 1];       
    int* maxCount = new int[mpisize];
    MPI_Allreduce(sendCount, maxCount, mpisize, MPI_INT, MPI_MAX, MPI_COMM_WORLD);// Reduce domain Xmin
    int maxSize = maxCount[0];
    for (int i = 1; i < mpisize; ++i) 
      if(maxCount[i] > maxSize)
        maxSize = maxCount[i]; 
    delete[] maxCount;
    std::vector<MPI_Request*> pendingRequests;                 //!< Buffer for non-blocking requests
    std::vector<T*>  pendingBuffers;                           //!< Buffer for non-blocking cell sends
    std::vector<std::vector<int> > path = getHypercubeMatrix(mpisize);
    recvB.clear();
    assert( (sizeof(sendB[0]) & 3) == 0 );
    int word = sizeof(sendB[0]) / 4;
    int* sendBuff = (int*)&sendB[0];    
    int logP = log2(mpisize);
    int previousSendBufferEnd = 0;
    int previousSendBufferSize = 0;
    std::vector<std::pair<int,T*> > sendRecvBuffer;
    std::vector<int> dataToSend;
    int ownDataIndex = 0; 
    int ownDataCursor = 0;   
    const int logPReceiveSize = mpisize >> 1;
    MPI_Request* rreq = new MPI_Request[logPReceiveSize];            
    T** irecvBuff = new T*[logPReceiveSize];    
    for (int i = 0; i < logPReceiveSize; ++i) 
      irecvBuff[i] = new T(maxSize);    

    for (int i = 0; i < logP; ++i) {
      int sendRecvSize = mpisize>>(i+1);      
      int sendRecvRank = path[mpirank][i];
      std::vector<int> route;
      route.push_back(sendRecvRank);
      getMessageInfo(path, route,sendRecvRank, i+1,  logP);
      MPI_Request request;
      for (int s = 0; s < logPReceiveSize; ++s) {
        MPI_Irecv((int*)&(*irecvBuff[s])[0], maxSize*word , MPI_INT, sendRecvRank, MPI_ANY_TAG, MPI_COMM_WORLD, &rreq[s]);
      }
      // Sending rank's own data
      MPI_Isend(sendBuff + sendDispl[sendRecvRank]*word, sendCount[sendRecvRank]*word,MPI_INT,sendRecvRank,sendRecvRank,MPI_COMM_WORLD,&request);      
      for (int s = 1; s < sendRecvSize; ++s) {                     
          MPI_Isend(sendBuff + sendDispl[route[s]]*word, sendCount[route[s]]*word,MPI_INT,sendRecvRank,route[s],MPI_COMM_WORLD,&request);        
      }      

      // end of sending my rank's own data
      int dataToSendSize = dataToSend.size();
      int rerouteSize = 0;
      // sending other rank's data
      for(int s = 0; s < dataToSendSize; ++s) {
        int pp = dataToSend[s];
        std::pair<int, T*> const& dataPair = sendRecvBuffer[pp];
        T* data = dataPair.second;
        int dest = dataPair.first;
        bool reroute = true;
        if(dest != sendRecvRank)
          reroute = searchMessagePath(path, sendRecvRank,i+1,logP,dest);
        if(reroute) {            
          int size = data->size();
          MPI_Request* pendingRequest = new MPI_Request();
          pendingRequests.push_back(pendingRequest);          
          pendingBuffers.push_back(data);
          MPI_Isend((int*)&(*data)[0],size*word,MPI_INT,sendRecvRank,dest,MPI_COMM_WORLD,pendingRequest);
          dataToSend.erase(dataToSend.begin() + s); 
          dataToSendSize--;
          s--;   
          rerouteSize++;       
        } 
      }
      for (int pp = previousSendBufferEnd; pp < previousSendBufferSize; ++pp) {
        std::pair<int, T*> const& dataPair = sendRecvBuffer[pp];
        T* data = dataPair.second;
        int dest = dataPair.first;
        bool reroute = true;
        if(dest != sendRecvRank)
          reroute = searchMessagePath(path, sendRecvRank,i+1,logP,dest);
        if(reroute) {  
          int size = data->size();
          MPI_Request* pendingRequest = new MPI_Request();
          pendingRequests.push_back(pendingRequest);          
          pendingBuffers.push_back(data);
          MPI_Isend((int*)&(*data)[0],size*word,MPI_INT,sendRecvRank,dest,MPI_COMM_WORLD,pendingRequest);
          rerouteSize++;
        } else {
          dataToSend.push_back(pp);
        }
      }
      for (int s = 0; s < logPReceiveSize; ++s) {
        int index; 
        MPI_Status status;
        MPI_Waitany(logPReceiveSize, rreq, &index, &status);
        int intCount; 
        MPI_Get_count(&status, MPI_INT, &intCount);     
        int dest = status.MPI_TAG;
        int count = intCount/word;
        if(dest == mpirank) {
          if(count>0) {
            recvB.insert(recvB.end(),irecvBuff[index]->begin(),irecvBuff[index]->begin() + count);
            recvDispl[irecvBuff[index]->begin()->IRANK] = ownDataCursor;
            recvCount[irecvBuff[index]->begin()->IRANK] = count;
            ownDataCursor+=count;
          }  
          ownDataIndex++;
        } else {         
          sendRecvBuffer.push_back(std::pair<int,T*>(dest,new T(irecvBuff[index]->begin(),irecvBuff[index]->begin() + count)));
        }        
      }      
      previousSendBufferEnd = previousSendBufferSize;
      previousSendBufferSize = sendRecvBuffer.size();    
      deallocateCompletedRequests(pendingRequests, pendingBuffers);
    }
    deallocateCompletedRequests(pendingRequests, pendingBuffers, true,true);
    assert(ownDataIndex == mpisize - 1);
    assert(pendingBuffers.size() == 0);
    typename T::iterator localBuffer = sendB.begin() + sendDispl[mpirank];
    recvB.insert(recvB.end(),localBuffer, localBuffer + sendCount[mpirank]);
    recvDispl[mpirank] = ownDataCursor;
    recvCount[mpirank] = sendCount[mpirank];
    sendB.clear();
    for (int i = 0; i < logPReceiveSize; ++i) {
      delete irecvBuff[i];  
    }
    delete[] irecvBuff;    
  }

  template<typename T>
  void alltoallv_heirarchical(T& sendB, T& recvB, int* recvDispl, int* recvCount, int* sendDispl, int* sendCount) { 
    int logP = log2(mpisize);
    int depth = 0;
    int dataSize = recvDispl[mpisize - 1] + recvCount[mpisize - 1];    
    if (sendB.size() == sendCount[mpirank] && dataSize == sendCount[mpirank]) {
      recvB = sendB;
      return;
    }    
    recvB.resize(dataSize);
    int word = sizeof(sendB[0]) / 4;    
    int* sendBuff = (int*)&sendB[0];
    int* recvBuff = (int*)&recvB[0];
    int* groupSendCount = new int[mpisize];
    int* groupRecvCount = new int[mpisize];
    int* groupSendDispl = new int[mpisize];
    int* groupRecvDispl = new int[mpisize];
    for (int i = 0; i < logP; ++i) {
      int partitions = mpisize<<i;      
      if(depth == 0) {
        int dest = mpirank ^ 1; 
        MPI_Request rreq;
        MPI_Request sreq;
        MPI_Irecv(recvBuff + recvDispl[dest]*word,
                    recvCount[dest]*word, MPI_INT, dest, dest, MPI_COMM_WORLD, &rreq);
        MPI_Isend(sendBuff + sendDispl[dest]*word,
                    sendCount[dest]*word, MPI_INT, dest, mpirank, MPI_COMM_WORLD, &sreq);        
        MPI_Wait(&sreq,MPI_STATUS_IGNORE);
        MPI_Wait(&rreq,MPI_STATUS_IGNORE);
      }
      else { 
        MPI_Comm commSplit;
        MPI_Comm interComm;
        int split = mpirank >> depth;
        MPI_Comm_split(MPI_COMM_WORLD, split, mpirank, &commSplit);                                
        int intraCount = mpirank >> (depth+1);
        int dest = mpirank ^ (1<<i); 
        int lead_dest = (dest >> depth) << depth ;
        MPI_Intercomm_create(commSplit, 0, MPI_COMM_WORLD, lead_dest, 
                            intraCount, &interComm); 
        int sendRecvSize = 1 << depth;      
        for (int i = 0; i < sendRecvSize; ++i) {
          groupSendCount[i] = sendCount[lead_dest+i] * word;          
          groupRecvCount[i] = recvCount[lead_dest+i] * word;
          groupSendDispl[i] = sendDispl[lead_dest+i] * word;
          groupRecvDispl[i] = recvDispl[lead_dest+i] * word;
        }
        MPI_Alltoallv(sendBuff,groupSendCount,groupSendDispl,MPI_INT, recvBuff, groupRecvCount, groupRecvDispl,MPI_INT, interComm);
        MPI_Comm_free(&commSplit); 
        MPI_Comm_free(&interComm); 
      }      
      depth++;
    }
    delete[] groupRecvCount;        
    delete[] groupSendCount;
    delete[] groupRecvDispl;
    delete[] groupSendDispl;
    typename T::iterator localBuffer = sendB.begin() + sendDispl[mpirank];
    std::copy(localBuffer, localBuffer + sendCount[mpirank], recvB.begin() + recvDispl[mpirank]);  
  }       
  

  //! Exchange send count for cells
  void alltoall(Cells& cells) {
    MPI_Alltoall(sendCellCount, 1, MPI_INT,                   // Communicate send count to get receive count
                 recvCellCount, 1, MPI_INT, MPI_COMM_WORLD);
    recvCellDispl[0] = 0;                                     // Initialize receive displacements
    for (int irank = 0; irank < mpisize - 1; irank++) {       // Loop over ranks
      recvCellDispl[irank + 1] = recvCellDispl[irank] + recvCellCount[irank]; //  Set receive displacement
    }                                                         // End loop over ranks    
    C_iter c0 = cells.begin();
    for (C_iter cc=c0; cc!= cells.end(); ++cc) {
      cc->IRANK = mpirank;
    }
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

  template<typename MPI_T, typename BUFF_T>
  void deallocateCompletedRequests(MPI_T& requests, BUFF_T& data, bool eraseRequests=true, bool finalize = false) {
    int count = requests.size();
    for (int i = 0; i < count; ++i) {
      if (*(requests[i]) != MPI_REQUEST_NULL) {
        int flag = 0;
        if(!finalize) {          
          MPI_Test(requests[i], &flag, MPI_STATUS_IGNORE);          
        } else {
          MPI_Wait(requests[i], MPI_STATUS_IGNORE);
          flag = 1;          
        }
        if(flag) {          
          delete data[i];
          if(eraseRequests) {
            requests.erase(requests.begin() + i);
            data.erase(data.begin() + i);
            count--;
            i--;
          }
        }
      } else {
        delete data[i];
        if(eraseRequests) {
          requests.erase(requests.begin() + i);
          data.erase(data.begin() + i);
          count--;
          i--;
        } 
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

  template<typename BoundsArr>
  NeighborList setNeighbors(BoundsArr minBounds, BoundsArr maxBounds, Bounds globalBounds) {
    NeighborList neighbors(mpisize);
    real_t const leeway[3] = {
      0.001*(globalBounds.Xmax[0] - globalBounds.Xmin[0]), 
      0.001*(globalBounds.Xmax[1] - globalBounds.Xmin[1]),
      0.001*(globalBounds.Xmax[2] - globalBounds.Xmin[2])
    };
    for (int i = 0; i < mpisize; ++i) {
      real_t* localXmin = minBounds[i];
      real_t* localXmax = maxBounds[i];    
      for (int irank = 0; irank < mpisize; ++irank) {
         if(i!=irank) {
          int dim;
          for (dim = 0; dim < 3; ++dim) {            
            if((localXmin[dim] - leeway[dim] > maxBounds[irank][dim]) || (minBounds[irank][dim] - leeway[dim] > localXmax[dim])) break;
          }
          if(dim == 3) neighbors[i].push_back(irank);
        }
      }      
    }  
    return neighbors;
  }

  template <typename TraversedSet>
  void updateNeighbors(NeighborList const& neighbors, int parent, int depth, TraversedSet& traversedSet, NeighborList& output) {
    if(traversedSet.find(parent) == traversedSet.end()) {
      traversedSet.insert(parent);
      if(output.size() > depth) {
        output[depth].push_back(parent);
      }
      else {
        output.push_back(std::vector<int>());
        output[depth].push_back(parent);
      }
      std::vector<int> r_neighbors = neighbors[parent];
      for (int i = 0; i < r_neighbors.size(); ++i) {      
        updateNeighbors(neighbors, r_neighbors[i], depth+1 , traversedSet, output);          
      }
    }
  }


  void setSendRecvTrees(CommGraph const& graph, int const& level, std::map<int, TreeLevel>& sendList, std::map<int, TreeLevel>& recvList) {
    assert(level >= 2);    
    if(level < graph[mpirank].size()) {
      TreeLevel treeLevel = graph[mpirank][level];    
      for (int i = 0; i < treeLevel.size(); ++i) {  
        GraphNode node = treeLevel[i];
        if(recvList.find(node.source) == recvList.end()) recvList[node.source] = TreeLevel();          
        recvList[node.source].push_back(node);
      }
    }    
    for (int irank = 0; irank < mpisize; ++irank) {
      if(irank!=mpirank) {
        if(level < graph[irank].size()) {
          TreeLevel treeLevel = graph[irank][level];
          for (int i = 0; i < treeLevel.size(); ++i) {  
            GraphNode node = treeLevel[i];
            if(node.source == mpirank) {
              if(sendList.find(irank) == sendList.end()) 
                sendList[irank] = TreeLevel();          
              sendList[irank].push_back(node);
            }          
          }
        }   
      }       
    }
  }

  CommGraph setFMMTreeCommunicationGraph(NeighborList const& neighbors) {    
    CommGraph graph(mpisize);
    for (int source = 0; source < mpisize; ++source) {
      std::queue<GraphNode> traversalQ;
      CommTree output;    
      traversalQ.push(GraphNode(-1,source,0));
      std::set<int> traversedSet;
      int currentDepth = 0, elemToDepthInc = 1, nextElemToDepthInc = 0;
      bool visit = false;
      while(!traversalQ.empty() && traversedSet.size() < mpisize){
        GraphNode current = traversalQ.front();   
        traversalQ.pop();        
        visit = false;
        if(traversedSet.find(current.id) == traversedSet.end()) {
          visit = true;
          traversedSet.insert(current.id);
          if(output.size() > currentDepth) {
            output[currentDepth].push_back(current);
          }
          else {
            output.push_back(std::vector<GraphNode>());
            output[currentDepth].push_back(current);
          }      
        }
        if(visit) {
          std::vector<int> r_neighbors = neighbors[current.id];
          for (int i = 0; i < r_neighbors.size(); ++i) {    
              traversalQ.push(GraphNode(current.id, r_neighbors[i],(currentDepth>1)?current.source:current.id));
          }        
          nextElemToDepthInc += neighbors[current.id].size();   
        }
        elemToDepthInc--;
        if(elemToDepthInc == 0) {
          currentDepth++;
          elemToDepthInc = nextElemToDepthInc;
          nextElemToDepthInc = 0;
        }
      }
      graph[source] = output;
    }        
    return graph;
  }
  
public:
  //! Constructor
  TreeMPI(int _mpirank, int _mpisize, int _images) :
    mpirank(_mpirank), mpisize(_mpisize), images(_images),
    granularity(1)   
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
    }

  //! Temporary constructor added to prevent breaking compatibility with other kernels using this class
  TreeMPI(int _mpirank, int _mpisize, int _images, int _grainSize) :
    mpirank(_mpirank), mpisize(_mpisize), images(_images),
    granularity(_grainSize)
  {                                                           // Initialize variables
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


  void setSourcesNDestinations(NeighborList const& neighbors, int rank, int*& sources, int*& destinations,int& indegree,int& outdegree) {
    std::vector<int> inwardNeighbors = neighbors[rank];
    indegree = inwardNeighbors.size();
    sources = new int[indegree];
    for (int i = 0; i < indegree; ++i) {      
      sources[i] = inwardNeighbors[i];
    }
    std::vector<int> outwardNeighbors;
    for (int i = 0; i < neighbors.size(); ++i) {
      if(i != rank) {
        for (int j = 0; j < neighbors[i].size(); ++j)
        {
          if(neighbors[i][j] == rank){
            outwardNeighbors.push_back(i);
          }
        }
      }
    }
    outdegree = outwardNeighbors.size();
    destinations = new int[outdegree];
    for (int i = 0; i < outdegree; ++i) {
      destinations[i] = outwardNeighbors[i];
    }
  }

  NeighborList getAllNeighbors(Bounds globalBounds){
    return setNeighbors(allBoundsXmin, allBoundsXmax, globalBounds);       
  }

  template<typename Cycle>
  void commDistGraph(Cycle cycle, Bounds globalBounds) {
    logger::startTimer("Comm LET");
    real_t* localXmin = allBoundsXmin[mpirank];
    real_t* localXmax = allBoundsXmax[mpirank];        
    alltoall(sendBodies);
    alltoall(sendCells);
    int bodyDataSize = recvBodyDispl[mpisize - 1] + recvBodyCount[mpisize - 1];    
    int cellDataSize = recvCellDispl[mpisize - 1] + recvCellCount[mpisize - 1];    
    recvBodies.resize(bodyDataSize);
    recvCells.resize(cellDataSize);    
    int bodyWord = sizeof(sendBodies[0]) / 4;    
    int cellWord = sizeof(sendCells[0]) / 4;    
    Bodies tempBodies;
    Cells tempCells;
    NeighborList neighbors = setNeighbors(allBoundsXmin, allBoundsXmax, globalBounds);         
    CommGraph res = setFMMTreeCommunicationGraph(neighbors);  
    int* sources;
    int* destinations;
    int indegree;
    int outdegree;
    MPI_Comm graphComm;
    setSourcesNDestinations(neighbors,mpirank,sources,destinations,indegree, outdegree);
    int* inWeights = new int[indegree];
    int* outWeights = new int[outdegree];
    for (int i = 0; i < indegree; ++i)  inWeights[i] = 0;
    for (int i = 0; i < outdegree; ++i) outWeights[i] = 0;    
    MPI_Dist_graph_create_adjacent(MPI_COMM_WORLD, indegree, sources,inWeights,outdegree, destinations, outWeights, MPI_INFO_NULL, 0, &graphComm);
    int topo_type;
    MPI_Topo_test(graphComm, &topo_type ); // Get the type of Toplogy we are using.
    assert(topo_type==MPI_DIST_GRAPH);    
    MPI_Dist_graph_neighbors(graphComm, indegree, sources, inWeights,outdegree,destinations,outWeights);
    std::map<int,int> ranktoindex; 
    for (int i = 0; i < outdegree; ++i)
      ranktoindex[destinations[i]] = i;        
    
    int* sendBodyCountNeighbor = new int[outdegree];
    int* sendBodyDisplNeighbor = new int[outdegree];
    int* sendCellCountNeighbor = new int[outdegree];
    int* sendCellDisplNeighbor = new int[outdegree];
    for (int i = 0; i < outdegree; ++i) {
      sendBodyCountNeighbor[i] = sendBodyCount[destinations[i]] * bodyWord;   
      sendBodyDisplNeighbor[i] = sendBodyDispl[destinations[i]] * bodyWord;      
      sendCellCountNeighbor[i] = sendCellCount[destinations[i]] * cellWord; 
      sendCellDisplNeighbor[i] = sendCellDispl[destinations[i]] * cellWord; 
    }
    int* recvBodyCountNeighbor = new int[indegree];
    int* recvBodyDisplNeighbor = new int[indegree];
    int* recvCellCountNeighbor = new int[indegree];
    int* recvCellDisplNeighbor = new int[indegree];
    std::vector<int> receivedLevelTrees;
    int receivedTreesCount = indegree;
    for (int i = 0; i < indegree; ++i) {
      int irank = sources[i];                
      recvBodyCountNeighbor[i] = recvBodyCount[sources[i]] * bodyWord; 
      recvBodyDisplNeighbor[i] = recvBodyDispl[sources[i]] * bodyWord; 
      recvCellCountNeighbor[i] = recvCellCount[sources[i]] * cellWord; 
      recvCellDisplNeighbor[i] = recvCellDispl[sources[i]] * cellWord;       
      receivedLevelTrees.push_back(irank);
    }
    MPI_Neighbor_alltoallv((int*)&sendBodies[0], sendBodyCountNeighbor, sendBodyDisplNeighbor, MPI_INT, (int*)&recvBodies[0],recvBodyCountNeighbor,recvBodyDisplNeighbor,MPI_INT,graphComm);
    MPI_Neighbor_alltoallv((int*)&sendCells[0] , sendCellCountNeighbor, sendCellDisplNeighbor, MPI_INT, (int*)&recvCells[0],recvCellCountNeighbor,recvCellDisplNeighbor,MPI_INT,graphComm);                  
    typedef std::map<int, TreeLevel> nodemap;
    typedef nodemap::iterator map_iterator;
    int* sendCountB = new int[outdegree];
    int* sendDisplB = new int[outdegree];
    int* sendCountC = new int[outdegree];
    int* sendDisplC = new int[outdegree];
    
    int* recvCountB = new int[indegree];
    int* recvDisplB = new int[indegree];
    int* recvCountC = new int[indegree];
    int* recvDisplC = new int[indegree];
    int isend;
    int maxLevels = 0; 
    for (int i = 0; i < res.size(); ++i) {
      if(res[i].size() > maxLevels)
        maxLevels = res[i].size();
    }
    for (int level = 2; level < maxLevels; ++level) {      
      for (int i = 0; i < outdegree; ++i) {
        sendCountB[i] = 0; sendDisplB[i] = 0;    
        sendCountC[i] = 0; sendDisplC[i] = 0;    
      }
      
      std::map<int, TreeLevel> sendList, recvList;       
      setSendRecvTrees(res, level, sendList, recvList);  
      std::vector<int> sendRanks;    
      std::vector<std::vector<int> > sendLETs;
      std::vector<std::vector<int> > recvLETs;
      std::vector<std::vector<int> > cellCounts(outdegree);                  
      std::vector<std::vector<int> > bodyCounts(outdegree);                  

      for (int i = 0; i < outdegree; ++i) {
        int irank = destinations[i];
        sendRanks.push_back(irank);
        sendLETs.push_back(std::vector<int>());
        if(sendList.find(irank) != sendList.end()) {
          TreeLevel tLevel = sendList[irank];
          for (int j = 0; j < tLevel.size(); ++j) {
            GraphNode node = tLevel[j];
            sendLETs[i].push_back(node.id);
          }
        }
      }
      for (int i = 0; i < indegree; ++i) {
        int irank = sources[i];                
        recvLETs.push_back(std::vector<int>());
        if(recvList.find(irank) != recvList.end()) {
          TreeLevel tLevel = recvList[irank];
          for (int j = 0; j < tLevel.size(); ++j) {
            GraphNode node = tLevel[j];
            recvLETs[i].push_back(node.id);
          }
        }
      }
      if(sendList.size() > 0) {
        linkLET(receivedLevelTrees);
        setLET(sendLETs, cycle, sendRanks, cellCounts, bodyCounts);
      }
      std::vector<int> flatSendCellCounts;
      std::vector<int> flatSendBodyCounts;
      std::vector<int> flatRecvCellCounts;
      std::vector<int> flatRecvBodyCounts;      
      for (int i = 0; i < outdegree; ++i) {
        if(i!=0) sendDisplC[i] = sendDisplC[i-1] + sendCountC[i-1];
        if(i!=0) sendDisplB[i] = sendDisplB[i-1] + sendCountB[i-1];
        sendCountC[i] = cellCounts[i].size();
        for (int j = 0; j < sendCountC[i]; ++j) {
          int prev = (j!=0)? cellCounts[i][j-1]:0;
          flatSendCellCounts.push_back(cellCounts[i][j] - prev);
        }
        sendCountB[i] = bodyCounts[i].size();
        for (int j = 0; j < sendCountB[i]; ++j) {
          int prev = (j!=0)? bodyCounts[i][j-1]:0;
          flatSendBodyCounts.push_back(bodyCounts[i][j] - prev);
        }
      }
      MPI_Neighbor_alltoall(sendCountC, 1, MPI_INT, recvCountC, 1, MPI_INT, graphComm);
      MPI_Neighbor_alltoall(sendCountB, 1, MPI_INT, recvCountB, 1, MPI_INT, graphComm);
      recvDisplC[0] = 0;
      recvDisplB[0] = 0;
      for (int i = 1; i < indegree; ++i) {
        recvDisplC[i] = recvDisplC[i-1] + recvCountC[i-1];
        recvDisplB[i] = recvDisplB[i-1] + recvCountB[i-1];
      }
      flatRecvCellCounts.resize(recvDisplC[indegree-1] + recvCountC[indegree-1]);
      flatRecvBodyCounts.resize(recvDisplB[indegree-1] + recvCountB[indegree-1]);
      MPI_Neighbor_alltoallv((int*)&flatSendCellCounts[0], sendCountC, sendDisplC, MPI_INT, (int*)&flatRecvCellCounts[0],recvCountC,recvDisplC,MPI_INT,graphComm);
      MPI_Neighbor_alltoallv((int*)&flatSendBodyCounts[0], sendCountB, sendDisplB, MPI_INT, (int*)&flatRecvBodyCounts[0],recvCountB,recvDisplB,MPI_INT,graphComm);
      if(sendList.size() > 0) {
        for (int i = 0; i < outdegree; ++i) {
          sendBodyCountNeighbor[i] = sendBodyCount[destinations[i]] * bodyWord;   
          sendBodyDisplNeighbor[i] = sendBodyDispl[destinations[i]] * bodyWord;      
          sendCellCountNeighbor[i] = sendCellCount[destinations[i]] * cellWord; 
          sendCellDisplNeighbor[i] = sendCellDispl[destinations[i]] * cellWord; 
        }
      } else {
        for (int i = 0; i < outdegree; ++i) {
          sendBodyCountNeighbor[i] = 0;   
          sendBodyDisplNeighbor[i] = 0;      
          sendCellCountNeighbor[i] = 0; 
          sendCellDisplNeighbor[i] = 0; 
        }
      }
      recvCellDisplNeighbor[0] = 0;
      recvBodyDisplNeighbor[0] = 0;
      for (int i = 0; i < indegree; ++i) {
        recvBodyCountNeighbor[i] = 0;
        recvCellCountNeighbor[i] = 0;
        for(int j = 0; j < recvCountB[i]; ++j){
          recvBodyCountNeighbor[i] += flatRecvBodyCounts[recvDisplB[i] + j];
        }
        for(int j = 0; j < recvCountC[i]; ++j){
          recvCellCountNeighbor[i] += flatRecvCellCounts[recvDisplC[i] + j];
        }
        recvBodyCountNeighbor[i] *= bodyWord;
        recvCellCountNeighbor[i] *= cellWord;
        if(i!=0) {
          recvCellDisplNeighbor[i] = recvCellDisplNeighbor[i-1] + recvCellCountNeighbor[i-1];
          recvBodyDisplNeighbor[i] = recvBodyDisplNeighbor[i-1] + recvBodyCountNeighbor[i-1];
        }        
      }
      bodyDataSize = recvBodyDisplNeighbor[indegree - 1] + recvBodyCountNeighbor[indegree - 1];    
      cellDataSize = recvCellDisplNeighbor[indegree - 1] + recvCellCountNeighbor[indegree - 1];              
      tempBodies.resize(bodyDataSize/bodyWord);
      tempCells.resize(cellDataSize/cellWord);
      MPI_Neighbor_alltoallv((int*)&sendBodies[0], sendBodyCountNeighbor, sendBodyDisplNeighbor, MPI_INT, (int*)&tempBodies[0],recvBodyCountNeighbor,recvBodyDisplNeighbor,MPI_INT,graphComm);
      MPI_Neighbor_alltoallv((int*)&sendCells[0], sendCellCountNeighbor, sendCellDisplNeighbor, MPI_INT, (int*)&tempCells[0],recvCellCountNeighbor,recvCellDisplNeighbor,MPI_INT,graphComm);                    
      receivedLevelTrees.clear();
      for (int i = 0; i < indegree; ++i) {
        int cDispl = 0;
        int bDispl = 0;
        recvBodyDisplNeighbor[i]/=bodyWord;
        recvCellDisplNeighbor[i]/=cellWord;
        for (int j = 0; j < recvLETs[i].size(); ++j) {
          int irank = recvLETs[i][j];                    
          int cCount = flatRecvCellCounts[recvDisplC[i] + j];
          int bCount = flatRecvBodyCounts[recvDisplB[i] + j];
          receivedLevelTrees.push_back(irank);
          receivedTreesCount++;
          std::copy(tempBodies.begin() + recvBodyDisplNeighbor[i] + bDispl, tempBodies.begin() + recvBodyDisplNeighbor[i] + bDispl + bCount, recvBodies.begin() + recvBodyDispl[irank]);
          recvBodyCount[irank] = bCount;          
          C_iter const& CC  = tempCells.begin() + recvCellDisplNeighbor[i] + cDispl;
          C_iter const& CR0 = recvCells.begin() + recvCellDispl[irank];
          for (int k = 0; k < cCount; ++k) {
            C_iter recvIter = CR0 + k;
            *recvIter = *(CC+k);
            if(recvIter->NBODY > 0)
              recvIter->IBODY  = recvIter->IBODY - bDispl;
            if(recvIter->NCHILD > 0)
              recvIter->ICHILD = recvIter->ICHILD - cDispl;
            if(recvIter->IPARENT > 0)
              recvIter->IPARENT= recvIter->IPARENT - cDispl;
          }            
          assert(cCount<=recvCellCount[irank]);
          assert(bCount<=recvBodyCount[irank]);
          recvCellCount[irank] = cCount;
          recvBodyCount[irank] = bCount;
          cDispl+=cCount;
          bDispl+=bCount;
        }
      }    
    }
    assert(receivedTreesCount == mpisize - 1);     
    delete[] sources;
    delete[] destinations;
    delete[] inWeights;
    delete[] outWeights;
    delete[] sendBodyCountNeighbor;
    delete[] sendBodyDisplNeighbor;
    delete[] sendCellCountNeighbor;
    delete[] sendCellDisplNeighbor;
    delete[] recvBodyCountNeighbor;
    delete[] recvBodyDisplNeighbor;
    delete[] recvCellCountNeighbor;
    delete[] recvCellDisplNeighbor;
    delete[] sendCountB;
    delete[] sendDisplB;
    delete[] sendCountC;
    delete[] sendDisplC;
    delete[] recvCountB;
    delete[] recvDisplB;
    delete[] recvCountC;
    delete[] recvDisplC;
    logger::stopTimer("Comm LET");
    logger::printTime("Link LET");                           
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
        traverseLET(C0, C0, bounds, cycle, irank, ibody, icell, icell - 1, false); // Traverse tree to set LET
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
        traverseLET(C0, C0, bounds, cycle, irank, ibody, icell, icell - 1, true); // Traverse tree to set LET
      }                                                       //  Endif for current rank
    }                                                         // End loop over ranks
    logger::stopTimer("Set LET");                             // Stop timer
  }

  //! Set local essential tree to send to each process
  void setLET(std::vector<std::vector<int> > cells, vec3 cycle, std::vector<int> ranks, std::vector<std::vector<int> >& cellCounts, std::vector<std::vector<int> >& bodyCounts) {
    logger::startTimer("Set LET size");                       // Start timer
    Bounds bounds;                                            // Bounds of local subdomain
    for (int irank = 0; irank < mpisize; ++irank) {
      sendBodyDispl[0] = 0;                                     // Initialize body displacement vector
      sendCellDispl[0] = 0;                                     // Initialize cell displacement vector    
      sendBodyCount[0] = 0;                                     // Initialize body displacement vector
      sendCellCount[0] = 0;                                     // Initialize cell displacement vector    
    }  
    int size = ranks.size();
    cellCounts.resize(size);
    bodyCounts.resize(size);
    for (int i = 0; i < size; ++i) {
      int irank = ranks[i];
      if (i != 0) {
        int prevRank = ranks[i-1];
        sendBodyDispl[irank] = sendBodyDispl[prevRank] + sendBodyCount[prevRank]; // Update body displacement
        sendCellDispl[irank] = sendCellDispl[prevRank] + sendCellCount[prevRank]; // Update cell displacement
      }
      for (int d = 0; d < 3; d++) {                         //   Loop over dimensions
        bounds.Xmin[d] = allBoundsXmin[irank][d];           //    Local Xmin for irank
        bounds.Xmax[d] = allBoundsXmax[irank][d];           //    Local Xmax for irank
      }                                                     //   End loop over dimensions
      int ibody = 0;                                        //   Initialize send body's offset
      int icell = 1;                                        //   Initialize send cell's offset
      for (int j = 0; j < cells[i].size(); ++j) {
        int r = cells[i][j];
        C_iter C0 = recvCells.begin() + recvCellDispl[r];
        if (irank != mpirank) {                                 //  If not current rank and cell vector is not empty
          if (C0->NCHILD == 0) {                                //   If root cell is leaf
            addSendBody(C0, irank, ibody, icell - 1, false);    //    Add bodies to send
          }                                                     //   End if for root cell leaf
          traverseLET(C0, C0, bounds, cycle, irank, ibody, icell, 0, false); // Traverse tree to set LET
        }
        cellCounts[i].push_back(icell);
        bodyCounts[i].push_back(ibody);
        icell++;            
      }              
      sendBodyCount[irank] = ibody;                         //   Send body count for current rank
      sendCellCount[irank] = icell - 1;                     //   Send cell count for current rank
    }
    logger::stopTimer("Set LET size");                        // Stop timer
    logger::startTimer("Set LET");                            // Start timer
    int lastRank = ranks[size-1];
    int numSendBodies = sendBodyDispl[lastRank] + sendBodyCount[lastRank]; // Total number of send bodies
    int numSendCells = sendCellDispl[lastRank] + sendCellCount[lastRank]; // Total number of send cells
    sendBodies.resize(numSendBodies);                         // Clear send buffer for bodies
    sendCells.resize(numSendCells);                           // Clear send buffer for cells
    for (int i = 0; i < size; ++i) {
      int irank = ranks[i];
      for (int d = 0; d < 3; d++) {                         //   Loop over dimensions
        bounds.Xmin[d] = allBoundsXmin[irank][d];           //    Local Xmin for irank
        bounds.Xmax[d] = allBoundsXmax[irank][d];           //    Local Xmax for irank
      }                                                     //   End loop over dimensions
      int ibody = 0;                                        //   Initialize send body's offset
      int icell = 0;                                        //   Initialize send cell's offset
      for (int j = 0; j < cells[i].size(); ++j) {
        int r = cells[i][j];
        C_iter C0 = recvCells.begin() + recvCellDispl[r];
        if (irank != mpirank) {                                //  If not current rank and cell vector is not empty
          C_iter Csend = sendCells.begin() + sendCellDispl[irank] + icell;//   Send cell iterator
         *Csend = *C0;                                         //   Copy cell to send buffer
          Csend->NCHILD = Csend->NBODY = 0;                     //   Reset link to children and bodies
          icell++;                                              //   Increment send cell counter
          if (C0->NCHILD == 0) {                                //   If root cell is leaf
            addSendBody(C0, irank, ibody, icell - 1, true);    //    Add bodies to send
          }                                                     //   End if for root cell leaf
          traverseLET(C0, C0, bounds, cycle, irank, ibody, icell, icell - 1, true); // Traverse tree to set LET
        }
      }              
    }
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

  //! Link LET with received bodies and calcualte NBODY
  void linkLET(std::vector<int> ranks) {
    logger::startTimer("Link LET");                           // Start timer
    for (int c = 0; c < ranks.size(); ++c) {                 // Loop over ranks
      int irank = ranks[c];
      for (int i = 0; i < recvCellCount[irank]; i++) {        //  Loop over receive cells
        C_iter C = recvCells.begin() + recvCellDispl[irank] + i;//   Iterator of receive cell
        if (C->NBODY != 0) {                                  //   If cell has bodies
          C->BODY = recvBodies.begin() + recvBodyDispl[irank] + C->IBODY;// Iterator of first body
        }                                                     //   End if for bodies
      }                                                       //  End loop over receive cells
    }                                                         // End loop over ranks
    logger::stopTimer("Link LET",0);                          // End timer
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
#if EXAFMM_USE_H_ALLTOALL    
    alltoall(bodies);                                         // Send body count    
    alltoallv_heirarchical(bodies,recvBodies,recvBodyDispl,recvBodyCount,sendBodyDispl,sendBodyCount);    
#elif EXAFMM_USE_BUTTERFLY   
    for (int i = 0; i < mpisize; i++) {                       // Loop over ranks
      sendBodyCount[i] = 0;                                   //  Initialize send counts
      recvBodyCount[i] = 0;
      recvBodyDispl[i] = 0;
    }                                                         // End loop over ranks
    for (B_iter B = bodies.begin(); B != bodies.end(); B++) { // Loop over bodies
      assert(0 <= B->IRANK && B->IRANK < mpisize);            //  Check bounds for process ID
      sendBodyCount[B->IRANK]++;                              //  Fill send count bucket
      B->IRANK = mpirank;                                     //  Tag for sending back to original rank
    }                                                         // End loop over bodies
    sendBodyDispl[0]  = 0;                                    // Initialize send displacements
    for (int irank = 0; irank < mpisize - 1; irank++) {       // Loop over ranks
      sendBodyDispl[irank + 1] = sendBodyDispl[irank] + sendBodyCount[irank]; //  Set send displacement
    }
    alltoallv_p2p_hypercube(bodies,recvBodies,recvBodyDispl,recvBodyCount,sendBodyDispl,sendBodyCount);
// #elif EXAFMM_USE_ONESIDED
//     alltoall(bodies);                                     // Send body count       
//     alltoallv_p2p_onesided(bodies,recvBodies,recvBodyDispl,recvBodyCount,sendBodyDispl,sendBodyCount);
#else 
    alltoall(bodies);                                         // Send body count    
    alltoallv(bodies);                                        // Send bodies        
#endif    
    logger::stopTimer("Comm partition");                      // Stop timer
    return recvBodies;                                        // Return received bodies
  }

  //! Send bodies
  Bodies commBodies() {
    logger::startTimer("Comm LET bodies");                    // Start timer
#if EXAFMM_USE_H_ALLTOALL 
    alltoall(sendBodies);                                     // Send body count 
    alltoallv_heirarchical(sendBodies,recvBodies,recvBodyDispl,recvBodyCount,sendBodyDispl,sendBodyCount);
#elif EXAFMM_USE_BUTTERFLY    
    for (int i = 0; i < mpisize; i++) {                       // Loop over ranks
      sendBodyCount[i] = 0;                                   //  Initialize send counts
      recvBodyCount[i] = 0;
      recvBodyDispl[i] = 0;
    }                                                       // End loop over ranks
    for (B_iter B = sendBodies.begin(); B != sendBodies.end(); B++) { // Loop over bodies
      assert(0 <= B->IRANK && B->IRANK < mpisize);            //  Check bounds for process ID
      sendBodyCount[B->IRANK]++;                              //  Fill send count bucket
      B->IRANK = mpirank;                                     //  Tag for sending back to original rank
    }  
    sendBodyDispl[0] = 0;                                     // Initialize send displacements
    for (int irank = 0; irank < mpisize - 1; irank++) {       // Loop over ranks
      sendBodyDispl[irank + 1] = sendBodyDispl[irank] + sendBodyCount[irank]; //  Set send displacement
    }      
    alltoallv_p2p_hypercube(sendBodies,recvBodies,recvBodyDispl,recvBodyCount,sendBodyDispl,sendBodyCount);
   #elif EXAFMM_USE_ONESIDED
     alltoall(sendBodies);                                     // Send body count       
     alltoallv_p2p_onesided(sendBodies,recvBodies,recvBodyDispl,recvBodyCount,sendBodyDispl,sendBodyCount);
#else    
    alltoall(sendBodies);                                     // Send body count       
    alltoallv(sendBodies);                                    // Send bodies
#endif        
    logger::stopTimer("Comm LET bodies");                     // Stop timer
    return recvBodies;                                        // Return received bodies
  }

  //! Send cells
  void commCells() {
    logger::startTimer("Comm LET cells");                     // Start timer
#if EXAFMM_USE_H_ALLTOALL
    alltoall(sendCells);                                      // Send cell count
    alltoallv_heirarchical(sendCells,recvCells,recvCellDispl,recvCellCount,sendCellDispl,sendCellCount);      
#elif EXAFMM_USE_BUTTERFLY
    C_iter c0 = sendCells.begin();
    for (C_iter cc=c0; cc!= sendCells.end(); ++cc) {
      cc->IRANK = mpirank;
    }
    for (int i = 0; i < mpisize; i++) {                       // Loop over ranks
      recvCellCount[i] = 0;
      recvCellDispl[i] = 0;
    }  
    alltoallv_p2p_hypercube(sendCells,recvCells,recvCellDispl,recvCellCount,sendCellDispl,sendCellCount);
 #elif EXAFMM_USE_ONESIDED
    alltoall(sendCells);                                     // Send body count       
    alltoallv_p2p_onesided(sendCells,recvCells,recvCellDispl,recvCellCount,sendCellDispl,sendCellCount);      
#else 
    alltoall(sendCells);                                      // Send cell count            
    alltoallv(sendCells);                                     // Senc cells
#endif    
    logger::stopTimer("Comm LET cells");                      // Stop timer
  }

  //! Copy recvBodies to bodies
  Bodies getRecvBodies() {
    return recvBodies;                                        // Return recvBodies
  }

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
};
}
#endif
