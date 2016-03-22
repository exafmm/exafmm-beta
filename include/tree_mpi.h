#ifndef tree_mpi_h
#define tree_mpi_h
#include "kernel.h"
#include "logger.h"

namespace exafmm {
  //! Handles all the communication of local essential trees
  class TreeMPI {
  protected:
    typedef std::map<int, Cell>    CellMap;                     //!< Type of cell hash map
    typedef std::map<int, Cells>   ChildCellsMap;               //!< Type of child cells hash map
    typedef std::map<int, Bodies>  BodiesMap;                   //!< Type of bodies hash map
    const int mpirank;                                          //!< Rank of MPI communicator
    const int mpisize;                                          //!< Size of MPI communicator
    const int images;                                           //!< Number of periodic image sublevels
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
    const char sendbit;                                         //!< Send bit value
    const char receivebit;                                      //!< Recv bit value
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
    CellMap cellsMap;                                           //!< mapper to keep track of received individual cells
    ChildCellsMap childrenMap;                                  //!< mapper to keep track of received child cells
    BodiesMap bodyMap;                                          //!< mapper to keep track of received bodies
    Cells LETCells;                                             //!< LET cells
    Bodies LETBodies;                                           //!< LET bodies
    C_iter Ci0;                                                 //!< Iterator of first target cell
    C_iter Cj0;                                                 //!< Iterator of first source cell
    int nspawn;                                                 //!< Threshold of NBODY for spawning new threads
    int terminated;                                             //!< Number of terminated requests
    size_t hitCount;                                            //!< Number of remote requests
    int sendIndex;                                              //!< LET Send cursor
    int cellWordSize;                                           //!< Number of words in cell struct
    int bodyWordSize;                                           //!< Number of words in body struct
    int granularity;                                            //!< The granularity of communication
    std::vector<MPI_Request*> pendingRequests;
    std::vector<Cells*> sendBuffers;

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
    void alltoallv_p2p(Bodies& bodies) {
      int dataSize = recvBodyDispl[mpisize - 1] + recvBodyCount[mpisize - 1];
      if (bodies.size() == sendBodyCount[mpirank] && dataSize == sendBodyCount[mpirank])
        return;
      recvBodies.resize(dataSize);
      assert( (sizeof(bodies[0]) & 3) == 0 );
      int word = sizeof(bodies[0]) / 4;
      for (int irank = 0; irank < mpisize; ++irank) {
        sendBodyCount[irank] *= word;
        sendBodyDispl[irank] *= word;
        recvBodyCount[irank] *= word;
        recvBodyDispl[irank] *= word;
      }
      int* sendBuff = (int*)&bodies[0];
      int* recvBuff = (int*)&recvBodies[0];
      size_t sendSize = 0;
      size_t recvSize = 0;
      MPI_Request* rreq   = new MPI_Request[mpisize - 1];
      MPI_Request* sreq   = new MPI_Request[mpisize - 1];
      MPI_Status* rstatus = new MPI_Status[mpisize - 1];
      MPI_Status* sstatus = new MPI_Status[mpisize - 1];

      for (int irank = 0; irank < mpisize; ++irank) {
        if (irank != mpirank) {
          if (recvBodyCount[irank] > 0) {
            MPI_Irecv(recvBuff + recvBodyDispl[irank],
                      recvBodyCount[irank], MPI_INT, irank, irank, MPI_COMM_WORLD, &rreq[recvSize]);
            recvSize++;
          }
        }
      }
      for (int irank = 0; irank < mpisize; ++irank) {
        if (irank != mpirank) {
          if (sendBodyCount[irank] > 0) {
            MPI_Isend(sendBuff + sendBodyDispl[irank],
                      sendBodyCount[irank], MPI_INT, irank, mpirank, MPI_COMM_WORLD, &sreq[sendSize]);
            sendSize++;
          }
        }
      }
      for (int irank = 0; irank < mpisize; ++irank) {
        sendBodyCount[irank] /= word;
        sendBodyDispl[irank] /= word;
        recvBodyCount[irank] /= word;
        recvBodyDispl[irank] /= word;
      }
      B_iter localBuffer = bodies.begin() + sendBodyDispl[mpirank];
      std::copy(localBuffer, localBuffer + sendBodyCount[mpirank], recvBodies.begin() + recvBodyDispl[mpirank]);
      MPI_Waitall(sendSize, &sreq[0], &sstatus[0]);
      MPI_Waitall(recvSize, &rreq[0], &rstatus[0]);
      delete[] rreq;
      delete[] sreq;
      delete[] rstatus;
      delete[] sstatus;
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
        
    inline bool processIncomingMessage(int tag, int source) {
      if ((tag & directionmask) == receivebit)
        return false;
      hitCount++;
      uint8_t msgType = getMessageType(tag);
      int nullTag = encryptMessage(1, nulltag, 0, receivebit);
      MPI_Request* request = new MPI_Request();
      char null;
      bool dataReady = false;
      if (msgType == childcelltag) {
        int childID;
        int level = getMessageLevel(tag);
        MPI_Recv(&childID, 1, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        toggleDirection(tag);
        Cell& cell = LETCells[childID];
        Cell& pcell = LETCells[cell.IPARENT];
        if (pcell.NCHILD > 0) {
          uint16_t grainSize = getGrainSize(tag);
          size_t sendingSize = 0;
          Cells* sendBuffer = new Cells();
          packSubtree(sendBuffer, LETCells.begin(), pcell, grainSize, sendingSize, source);
          MPI_Isend((int*)&(*sendBuffer)[0],  cellWordSize * sendingSize, MPI_INT, source, tag, MPI_COMM_WORLD, request);
          sendBuffers.push_back(sendBuffer);
          pendingRequests.push_back(request);        
          sendIndex += sendingSize;
          dataReady = true;
        }
        dataReady = false;
      }
      else if (msgType == bodytag) {
        int recvBuff[2];
        MPI_Recv(recvBuff, 2, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        toggleDirection(tag);
        MPI_Isend((int*)&LETBodies[recvBuff[0]], bodyWordSize * recvBuff[1], MPI_INT, source, tag, MPI_COMM_WORLD, request);
        dataReady = true;
      }
      else if (msgType == celltag) {
        int cellID;
        MPI_Recv(&cellID, 1, MPI_INT, source, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        toggleDirection(tag);
        MPI_Isend((int*)&LETCells[cellID], cellWordSize, MPI_INT, source, tag, MPI_COMM_WORLD, request);
        dataReady = true;
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
      deallocateCompletedRequests();
      return dataReady;    
    }

    void deallocateCompletedRequests() {
      size_t count = pendingRequests.size();
      for(int i=0;i<count;++i) {
        if(*(pendingRequests[i])!=MPI_REQUEST_NULL) {
          int flag;
          MPI_Test(pendingRequests[i],&flag,MPI_STATUS_IGNORE);
          if(flag) { 
            delete sendBuffers[i];        
            pendingRequests.erase(pendingRequests.begin() + i);
            sendBuffers.erase(sendBuffers.begin() + i);
            count--;
          }
        } else {
          pendingRequests.erase(pendingRequests.begin() + i);
          sendBuffers.erase(sendBuffers.begin() + i);
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
    real_t getRankDistance(C_iter C, vec3 Xperiodic, int rank) {
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
    void traverseRemote(C_iter Ci, C_iter Cj, bool mutual, real_t remote, size_t rank) {
      vec3 dX = Ci->X - Cj->X - kernel::Xperiodic;                // Distance vector from source to target
      real_t R2 = norm(dX);                                       // Scalar distance squared
      if (R2 > (Ci->R + Cj->R) * (Ci->R + Cj->R) * (1 - 1e-3)) {  // Distance is far enough
        kernel::M2L(Ci, Cj, false);                               //  M2L kernel
      } else if (Ci->NCHILD == 0 && Cj->NCHILD == 0) {            // Else if both cells are bodies
        double commtime;
        Bodies bodies = getBodies(Cj->IBODY, Cj->NBODY, Cj->LEVEL, rank, bodytag, commtime);
#if WEIGH_COM
        Ci->WEIGHT += commtime;
#endif
        if (bodies.size() > 0) {
          Cj->BODY = bodies.begin();
          kernel::P2P(Ci, Cj, false);                             //    P2P kernel for pair of cells
#if 0
          countWeight(Ci, remoteWeight);                          //   Increment P2P weight
#endif
        } else {
          kernel::M2L(Ci, Cj, false);                             //   M2L kernel
          //countWeight(Ci, remote,remoteWeight*0.25);            //   Increment M2L weight
        }
      } else {                                                    // Else if cells are close but not bodies
        splitCellRemote(Ci, Cj, mutual, remote, rank);            //  Split cell and call function recursively for child
      }                                                           // End if for multipole acceptance
    }

    //! Split cell and call traverse() recursively for child
    void splitCellRemote(C_iter Ci, C_iter Cj, bool mutual, real_t remote, size_t rank) {
      double commtime;
      if (Cj->NCHILD == 0) {                                      // If Cj is leaf
        assert(Ci->NCHILD > 0);                                   //  Make sure Ci is not leaf
        for (C_iter ci = Ci0 + Ci->ICHILD; ci != Ci0 + Ci->ICHILD + Ci->NCHILD; ci++) { // Loop over Ci's children
          traverseRemote(ci, Cj, mutual, remote, rank);           //   Traverse a single pair of cells
        }                                                         //  End loop over Ci's children
      } else if (Ci->NCHILD == 0) {                               // Else if Ci is leaf
        assert(Cj->NCHILD > 0);                                   //  Make sure Cj is not leaf
        Cells cells = getCell(Cj->ICHILD, Cj->NCHILD, Cj->LEVEL, rank, childcelltag, commtime);
#if WEIGH_COM
        Ci->WEIGHT += commtime;
#endif
        for (C_iter cj = cells.begin(); cj != cells.end(); ++cj) {// Loop over Cj's children
          traverseRemote(Ci, cj, mutual, remote, rank);           //   Traverse a single pair of cells
        }
      } else if (Ci->NBODY + Cj->NBODY >= nspawn) {               // Else if cells are still large
        Cells cells = getCell(Cj->ICHILD, Cj->NCHILD, Cj->LEVEL, rank, childcelltag, commtime);
#if WEIGH_COM
        Ci->WEIGHT += commtime;
#endif
        TraverseRemoteRange traverseRange(this, Ci0 + Ci->ICHILD, Ci0 + Ci->ICHILD + Ci->NCHILD, // Instantiate recursive functor
                                          cells.begin(), cells.end(), mutual, remote, rank);
        traverseRange();
      } else if (Ci->R >= Cj->R) {                                // Else if Ci is larger than Cj
        for (C_iter ci = Ci0 + Ci->ICHILD; ci != Ci0 + Ci->ICHILD + Ci->NCHILD; ci++) { // Loop over Ci's children
          traverseRemote(ci, Cj, mutual, remote, rank);           //   Traverse a single pair of cells
        }                                                         //  End loop over Ci's children
      } else {                                                    // Else if Cj is larger than Ci
        Cells cells = getCell(Cj->ICHILD, Cj->NCHILD, Cj->LEVEL, rank, childcelltag, commtime);
#if WEIGH_COM
        Ci->WEIGHT += commtime;
#endif
        for (C_iter cj = cells.begin(); cj != cells.end(); ++cj) {// Loop over Cj's children
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
      size_t rank;
      TraverseRemoteRange(TreeMPI * _treeMPI, C_iter _CiBegin, C_iter _CiEnd,// Constructor
                          C_iter _CjBegin, C_iter _CjEnd,
                          bool _mutual, real_t _remote, size_t _rank) :
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
            leftBranch();
            TraverseRemoteRange rightBranch(treeMPI, CiMid, CiEnd,//    Instantiate recursive functor
                                            CjMid, CjEnd, mutual, remote, rank);;
            rightBranch();                                        //    Ci:latter Cj:latter
          }
          {
            TraverseRemoteRange leftBranch(treeMPI, CiBegin, CiMid,//    Instantiate recursive functor
                                           CjMid, CjEnd, mutual, remote, rank);;
            leftBranch();
            TraverseRemoteRange rightBranch(treeMPI, CiMid, CiEnd,//    Instantiate recursive functor
                                            CjBegin, CjMid, mutual, remote, rank);;
            rightBranch();                                        //    Ci:latter Cj:former
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

    inline int encryptMessage(uint16_t grainSize, uint8_t requestType, uint8_t level, char direction) {
      int tag = int(grainSize);
      tag <<= requestshift; tag |= requestType;
      tag <<= levelshift; tag |= level;
      tag <<= directionshift; tag |= direction;
      return tag;
    }

    inline void decryptMessage(int tag, uint16_t& grainSize, uint8_t& requestType, uint8_t& level, char& direction) {
      direction = tag & directionmask;
      tag >>= directionshift;
      level = tag & levelmask;
      tag >>= levelshift;
      requestType = tag & requestmask;
      tag >>= requestshift;
      grainSize = tag & grainmask;
    }

    inline uint16_t getGrainSize(int const& tag) {
      return ((tag >> (requestshift + directionshift + levelshift)) & grainmask);
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

    inline void packSubtree(Cells* cells, C_iter const& C0, Cell const& C, uint16_t const& grainSize, size_t& index, int rank) {
      C_iter begin = C0 + C.ICHILD;
      C_iter end   = C0 + C.ICHILD + C.NCHILD;
      cells->insert(cells->end(), begin, end);      
      index += C.NCHILD;    
      for (C_iter cc = begin; cc < end; ++cc)
        if(index < grainSize)             
          packSubtree(cells, C0, *cc, grainSize, index, rank);    
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

    template <typename MapperType>
    void appendDataToCellMapper(MapperType& map, Cells const& cells, int& index, int id) {
      Cell parent = cells[id];
      if (index < cells.size()) {
        id = index;
        if (parent.NCHILD > 0) {        
          assert((index + parent.NCHILD) <= cells.size());
          assert(map.find(parent.ICHILD) == map.end());
          map[parent.ICHILD] = Cells(cells.begin() + index, cells.begin() + index + parent.NCHILD);
          index += parent.NCHILD;
          for (size_t i = 0; i < parent.NCHILD; ++i) {
            appendDataToCellMapper(map, cells, index, id++);
          }
        }
      }
    }

    void sendFlushRequest() {
      MPI_Request request;
      int tag = encryptMessage(1, flushtag, 0, sendbit);
      char null;
      for (int i = 0; i < mpisize; ++i)
        if (i != mpirank)
          MPI_Isend(&null, 1, MPI_CHAR, i, tag, MPI_COMM_WORLD, &request);
    }

  public:
    //! Constructor
    TreeMPI(int _mpirank, int _mpisize, int _images) :
      mpirank(_mpirank), mpisize(_mpisize), images(_images),
      terminated(0), hitCount(0), sendIndex(0),
      nulltag(1), celltag(2),
      childcelltag(3), bodytag(8),
      flushtag(9), maxtag(15), levelshift(5),
      requestshift(4), directionshift(1), levelmask(0x1F),
      requestmask(0xF), directionmask(0x1), grainmask(0xFFFF)
      , sendbit(1), receivebit(0) {                             // Initialize variables
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
      alltoallv(bodies);                                        // Send bodies
      logger::stopTimer("Comm partition");                      // Stop timer
      return recvBodies;                                        // Return received bodies
    }

    //! Send bodies
    Bodies commBodies() {
      logger::startTimer("Comm LET bodies");                    // Start timer
      alltoall(sendBodies);                                     // Send body count
      alltoallv(sendBodies);                                    // Send bodies
      logger::stopTimer("Comm LET bodies");                     // Stop timer
      return recvBodies;                                        // Return received bodies
    }

    //! Send cells
    void commCells() {
      logger::startTimer("Comm LET cells");                     // Start timer
      alltoall(sendCells);                                      // Send cell count
      alltoallv(sendCells);                                     // Senc cells
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

    void setSendLET(Cells& cells, Bodies& bodies) {
      LETCells = cells;
      LETBodies = bodies;    
      cellWordSize = sizeof(cells[0]) / 4;
      bodyWordSize = sizeof(bodies[0]) / 4;
    }

    void flushAllRequests() {
      sendFlushRequest();
      int ready;
      MPI_Status status;
      while (terminated < (mpisize - 1)) {
        ready = 0;
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &ready, &status);
        if (ready) {
          processIncomingMessage(status.MPI_TAG, status.MPI_SOURCE);
        }
      }        
    }

    //! Evaluate P2P and M2L using dual tree traversal
    void dualTreeTraversalRemote(Cells & icells, Bodies& ibodies, int mpirank, int mpisize, int spwan, int commgrain, real_t remote = 1) {
      if (icells.empty()) return;                                 // Quit if either of the cell vectors are empty
      nspawn = spwan;
      granularity = commgrain;
      logger::startTimer("Traverse Remote");                      // Start timer
      kernel::Xperiodic = 0;
      logger::initTracer();                                       // Initialize tracer    
      ///setSendLET(icells, ibodies);
      Ci0 = icells.begin();                                       // Set iterator of target root cell
      for (int i = 0; i < mpisize; ++i) {
        if (i != mpirank) {
          double commtime;
          Cells cells = getCell(0, 1, 0, i, celltag, commtime);
          assert(cells.size() > 0);
          traverseRemote(Ci0, cells.begin(), false, remote, i);
          logger::startTimer("Clear cache");                      // Start timer
          clearCellCache(i);
          logger::stopTimer("Clear cache", 0);                    // Start timer
        }
      }
      logger::printTime("Clear cache");
      logger::stopTimer("Traverse Remote");                       // Stop timer
      logger::writeTracer();                                      // Write tracer to file
    }

  //! Get bodies underlying a key and a level
  Bodies getBodies(int key, int nchild, size_t level, int rank, int requestType, double& commtime) {
    assert(rank != mpirank);
    commtime = 0.0;
    Bodies recvData;
    if (requestType == bodytag && bodyMap.find(key) ==  bodyMap.end()) {
      MPI_Request request;
      int tag = encryptMessage(1, requestType, level, sendbit);
#if CALC_COM_COMP
      logger::startTimer("Communication");
#endif
      int sendBuffer[2] = {key, nchild};
      MPI_Isend(sendBuffer, 2, MPI_INT, rank, tag, MPI_COMM_WORLD, &request);
      MPI_Status status;
      int recvRank = rank;
      int receivedTag;
      toggleDirection(tag);
      int ready = 0;
      while (!ready) {
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &ready, &status);
        if (ready) {
          if (tag != status.MPI_TAG || rank != status.MPI_SOURCE) {
            processIncomingMessage(status.MPI_TAG, status.MPI_SOURCE);
            ready = 0;
          }
        }
      }
      receivedTag = status.MPI_TAG;
      int responseType = getMessageType(receivedTag);
      int recvCount = 0;
      assert(responseType == bodytag || responseType == nulltag);
      if (responseType == bodytag) {
        MPI_Get_count(&status, MPI_INT, &recvCount);
        int bodyCount = recvCount / bodyWordSize;
        recvData.resize(bodyCount);
        MPI_Recv(&recvData[0], recvCount, MPI_INT, recvRank, receivedTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#if CALC_COM_COMP
        commtime = logger::stopTimer("Communication", 0);
#endif
        bodyMap[key] = recvData;
      } else {
        std::cout << "bodies not found for cell " << key << std::endl;
        char null;
        MPI_Recv(&null, 1, MPI_CHAR, recvRank, receivedTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#if CALC_COM_COMP
        commtime = logger::stopTimer("Communication", 0);
#endif
      }
    }
    else if (requestType == bodytag) {
      recvData = bodyMap[key];
    }
    return recvData;
  }

  //! Get cells underlying a key and a level
  Cells getCell(int key, int nchild, size_t level, int rank, int requestType, double& commtime) {
    assert(rank != mpirank);
    commtime = 0.0;
    Cells recvData;
    if ((requestType == childcelltag && childrenMap.find(key) == childrenMap.end()) ||
        (requestType == celltag      &&  cellsMap.find(key) ==  cellsMap.end())) {
      MPI_Request request;
      assert(requestType <= maxtag);
      int tag = encryptMessage(granularity, requestType, level, sendbit);
#if CALC_COM_COMP
      logger::startTimer("Communication");
#endif
      MPI_Isend(&key, 1, MPI_INT, rank, tag, MPI_COMM_WORLD, &request);
      MPI_Status status;
      int recvRank = rank;
      int receivedTag;
      toggleDirection(tag);
      int ready = 0;
      while (!ready) {
        MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &ready, &status);
        if (ready) {
          if (tag != status.MPI_TAG || rank != status.MPI_SOURCE) {
            processIncomingMessage(status.MPI_TAG, status.MPI_SOURCE);
            ready = 0;
          }
        }
      }
      receivedTag = status.MPI_TAG;
      int responseType = getMessageType(receivedTag);
      int recvCount = 0;
      assert(responseType != bodytag);
      if (responseType == childcelltag || responseType == celltag) {
        MPI_Get_count(&status, MPI_INT, &recvCount);
        int cellCount = recvCount / cellWordSize;
        recvData.resize(cellCount);
        MPI_Recv(&recvData[0], recvCount, MPI_INT, recvRank, receivedTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#if CALC_COM_COMP
        commtime = logger::stopTimer("Communication", 0);
#endif
        if (responseType == celltag) cellsMap[key] = recvData[0];
        else if (responseType == childcelltag)  {
          if (granularity > 1) {
            childrenMap[key] = Cells(recvData.begin(), recvData.begin() + nchild);
            int index = nchild;
            for (int i = 0; i < nchild; ++i) appendDataToCellMapper(childrenMap, recvData, index, i);
            recvData = childrenMap[key];
          }
          else childrenMap[key] = recvData;
        }
      } else if (responseType == nulltag) {
        char null;
        MPI_Recv(&null, 1, MPI_CHAR, recvRank, receivedTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
#if CALC_COM_COMP
        commtime = logger::stopTimer("Communication", 0);
#endif
      }
    }
    else if (requestType == celltag) {
      recvData.push_back(cellsMap[key]);
    }
    else if (requestType == childcelltag) {
      recvData = childrenMap[key];
    }
    return recvData;
  }

  void clearCellCache(int rank) {
    cellsMap.clear();
    bodyMap.clear();
    childrenMap.clear();
  }
};
}
#endif
