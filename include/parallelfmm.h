/*
Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
#ifndef parallelfmm_h
#define parallelfmm_h
#include "partition.h"

//! Handles all the communication of local essential trees
template<Equation equation>
class ParallelFMM : public Partition<equation> {
private:
  std::vector<int>    sendBodyCnt;                              //!< Vector of body send counts
  std::vector<int>    sendBodyDsp;                              //!< Vector of body send displacements
  std::vector<int>    recvBodyCnt;                              //!< Vector of body recv counts
  std::vector<int>    recvBodyDsp;                              //!< Vector of body recv displacements
  std::vector<int>    sendBodyRanks;                            //!< Vector of ranks to send bodies to
  std::vector<int>    sendBodyCellCnt;                          //!< Vector of send counts for cells of bodies
  std::vector<C_iter> sendBodyCells;                            //!< Vector of cell iterators for cells of bodies to send
  std::vector<int>    sendCellCnt;                              //!< Vector of cell send counts
  std::vector<int>    sendCellDsp;                              //!< Vector of cell send displacements
  std::vector<int>    recvCellCnt;                              //!< Vector of cell recv counts
  std::vector<int>    recvCellDsp;                              //!< Vector of cell recv displacements
  std::vector<vect>   xminAll;                                  //!< Buffer for gathering XMIN
  std::vector<vect>   xmaxAll;                                  //!< Buffer for gathering XMAX

  JBodies sendBodies;                                           //!< Send buffer for bodies
  JBodies recvBodies;                                           //!< Recv buffer for bodies
  JCells  sendCells;                                            //!< Send buffer for cells
  JCells  recvCells;                                            //!< Recv buffer for cells

public:
  using Kernel<equation>::printNow;                             //!< Switch to print timings
  using Kernel<equation>::startTimer;                           //!< Start timer for given event
  using Kernel<equation>::stopTimer;                            //!< Stop timer for given event
  using Kernel<equation>::sortBodies;                           //!< Sort bodies according to cell index
  using Kernel<equation>::sortCells;                            //!< Sort cells according to cell index
  using Kernel<equation>::R0;                                   //!< Radius of root cell
  using TreeStructure<equation>::buffer;                        //!< Buffer for MPI communication & sorting
  using TreeStructure<equation>::getLevel;                      //!< Get level from cell index
  using TreeStructure<equation>::getCenter;                     //!< Get cell center and radius from cell index
  using TreeStructure<equation>::bodies2twigs;                  //!< Group bodies into twig cells
  using TreeStructure<equation>::twigs2cells;                   //!< Link twigs bottomup to create all cells in tree
  using Partition<equation>::isPowerOfTwo;                      //!< If n is power of two return true
  using Partition<equation>::splitRange;                        //!< Split range and return partial range
  using Partition<equation>::print;                             //!< Print in MPI
  using Partition<equation>::LEVEL;                             //!< Level of the MPI process binary tree
  using Partition<equation>::XMIN;                              //!< Minimum position vector of bodies
  using Partition<equation>::XMAX;                              //!< Maximum position vector of bodies
  using Partition<equation>::nprocs;                            //!< Number of processes in the two split groups
  using Partition<equation>::color;                             //!< Color for hypercube communicators
  using Partition<equation>::key;                               //!< Key for hypercube communicators
  using Partition<equation>::MPI_COMM;                          //!< Hypercube communicators

private:
//! Gather bounds of other domain
  void gatherBounds() {
    xminAll.resize(MPISIZE);                                    // Resize buffer for gathering xmin
    xmaxAll.resize(MPISIZE);                                    // Resize buffer for gathering xmax
    sendBodyCnt.resize(MPISIZE);                                // Resize vector of body send counts
    sendBodyDsp.resize(MPISIZE);                                // Resize vector of body send displacements
    recvBodyCnt.resize(MPISIZE);                                // Resize vector of body recv counts
    recvBodyDsp.resize(MPISIZE);                                // Resize vector of body recv displacements
    sendCellCnt.resize(MPISIZE);                                // Resize vector of cell send counts
    sendCellDsp.resize(MPISIZE);                                // Resize vector of cell send displacements
    recvCellCnt.resize(MPISIZE);                                // Resize vector of cell recv counts
    recvCellDsp.resize(MPISIZE);                                // Resize vector of cell recv displacements
    MPI_Datatype MPI_TYPE = getType(XMIN[LEVEL][0]);            // Get MPI data type
    MPI_Allgather(&XMIN[LEVEL][0],3,MPI_TYPE,                   // Gather XMIN
                  &xminAll[0][0],3,MPI_TYPE,MPI_COMM_WORLD);
    MPI_Allgather(&XMAX[LEVEL][0],3,MPI_TYPE,                   // Gather XMAX
                  &xmaxAll[0][0],3,MPI_TYPE,MPI_COMM_WORLD);
  }

//! Get neighbor ranks to send to
  void getSendRank(Cells &cells) {
    sendBodyRanks.clear();                                      // Clear send ranks
    sendBodyCellCnt.clear();                                    // Clear send counts
    sendBodyCells.clear();                                      // Clear send body cells
    int oldsize = 0;                                            // Per rank offset of the number of cells to send
    for( int irank=0; irank!=MPISIZE; ++irank ) {               // Loop over ranks
      int ic = 0;                                               //  Initialize neighbor dimension counter
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimensions
        if(xminAll[irank][d] < 2 * XMAX[LEVEL][d] - XMIN[LEVEL][d] &&// If the two domains are touching or overlapping
           xmaxAll[irank][d] > 2 * XMIN[LEVEL][d] - XMAX[LEVEL][d]) {// in all dimensions, they are neighboring domains
          ic++;                                                 //    Increment neighbor dimension counter
        }                                                       //   Endif for overlapping domains
      }                                                         //  End loop over dimensions
      ic = 3;
      if( ic == 3 && irank != MPIRANK ) {                       //  If ranks are neighbors in all dimensions
        for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {    //   Loop over cells
          if( C->NCHILD == 0 ) {                                //    If cell is a twig
            bool send = false;                                  //     Initialize logical for sending
            if( IMAGES == 0 ) {                                 //     If free boundary condition
              Xperiodic = 0;                                    //      Set periodic coordinate offset
              real R = getDistance(C,xminAll[irank],xmaxAll[irank]);//  Get distance to other domain
              send |= CLET * C->R > THETA * R - EPS2;           //      If the cell seems close enough for P2P
            } else {                                            //     If periodic boundary condition
              for( int ix=-1; ix<=1; ++ix ) {                   //      Loop over x periodic direction
                for( int iy=-1; iy<=1; ++iy ) {                 //       Loop over y periodic direction
                  for( int iz=-1; iz<=1; ++iz ) {               //        Loop over z periodic direction
                    Xperiodic[0] = ix * 2 * R0;                 //         Coordinate offset for x periodic direction
                    Xperiodic[1] = iy * 2 * R0;                 //         Coordinate offset for y periodic direction
                    Xperiodic[2] = iz * 2 * R0;                 //         Coordinate offset for z periodic direction
                    real R = getDistance(C,xminAll[irank],xmaxAll[irank]);// Get distance to other domain
                    send |= CLET * C->R > THETA * R - EPS2;     //         If the cell seems close enough for P2P
                  }                                             //        End loop over z periodic direction
                }                                               //       End loop over y periodic direction
              }                                                 //      End loop over x periodic direction
            }                                                   //     Endif for periodic boundary condition
            if( send ) {                                        //     If the cell seems close enough for P2P
              sendBodyCells.push_back(C);                       //      Add cell iterator to scells
            }                                                   //     Endif for cell distance
          }                                                     //    Endif for twigs
        }                                                       //   End loop over cells
        sendBodyRanks.push_back(irank);                         //   Add current rank to sendBodyRanks
        sendBodyCellCnt.push_back(sendBodyCells.size()-oldsize);//   Add current cell count to sendBodyCellCnt
        oldsize = sendBodyCells.size();                         //   Set new offset for cell count
      }                                                         //  Endif for neighbor ranks
    }                                                           // End loop over ranks
  }

//! Get size of data to send
  void getSendCount(bool comm=true) {
    int ic = 0, ssize = 0;                                      // Initialize counter and offset for scells
    sendBodyCnt.assign(MPISIZE,0);                              // Initialize send count
    sendBodyDsp.assign(MPISIZE,0);                              // Initialize send displacement
    for( int i=0; i!=int(sendBodyRanks.size()); ++i ) {         // Loop over ranks to send to & recv from
      int irank = sendBodyRanks[i];                             //  Rank to send to & recv from
      for( int c=0; c!=sendBodyCellCnt[i]; ++c,++ic ) {         //  Loop over cells to send to that rank
        C_iter C = sendBodyCells[ic];                           //   Set cell iterator
        for( B_iter B=C->LEAF; B!=C->LEAF+C->NDLEAF; ++B ) {    //   Loop over bodies in that cell
          JBody body;                                           //    Set compact body type for sending
          body.ICELL = B->ICELL;                                //    Set cell index of compact body type
          body.X     = B->X;                                    //    Set position of compact body type
          body.SRC   = B->SRC;                                  //    Set source values of compact body type
          sendBodies.push_back(body);                           //    Push it into the send buffer
        }                                                       //   End loop over bodies
      }                                                         //  End loop over cells
      sendBodyCnt[irank] = sendBodies.size()-ssize;             //  Set send count of current rank
      sendBodyDsp[irank] = ssize;                               //  Set send displacement of current rank
      ssize += sendBodyCnt[irank];                              //  Increment offset for vector scells
    }                                                           // End loop over ranks
    if( comm ) {                                                // If communication is necessary
      MPI_Alltoall(&sendBodyCnt[0],1,MPI_INT,&recvBodyCnt[0],1,MPI_INT,MPI_COMM_WORLD);// Communicate the send counts
      int rsize = 0;                                            // Initialize total recv count
      for( int i=0; i!=MPISIZE; ++i ) {                         // Loop over ranks to recv from
        recvBodyDsp[i] = rsize;                                 //  Set recv displacements
        rsize += recvBodyCnt[i];                                //  Accumulate recv counts
      }                                                         // End loop over ranks to recv from
      recvBodies.resize(rsize);                                 // Resize recv buffer
    }
  }

//! Communicate cells by one-to-one MPI_Alltoallv
  void commBodiesAlltoall() {
    assert(isPowerOfTwo(MPISIZE));                              // Make sure the number of processes is a power of two
    int bytes = sizeof(sendBodies[0]);                          // Byte size of jbody structure
    int *scntd = new int [MPISIZE];                             // Permuted send count
    int *rcntd = new int [MPISIZE];                             // Permuted recv count
    int *rdspd = new int [MPISIZE];                             // Permuted recv displacement
    int *irev  = new int [MPISIZE];                             // Map original to compressed index
    JBodies sendBuffer = sendBodies;                            // Send buffer
    JBodies recvBuffer;                                         // Recv buffer
    for( int l=0; l!=LEVEL; ++l ) {                             // Loop over levels of N-D hypercube communication
      int npart = 1 << (LEVEL - l - 1);                         // Size of partition block
      int scnt2[2], sdsp2[2], rcnt2[2], rdsp2[2];               // Send/recv counts/displacements per level
      int ic = 0;                                               // Initialize counter
      for( int i=0; i!=2; ++i ) {                               // Loop over the two blocks
        scnt2[i] = 0;                                           //  Initialize send count per level
        for( int irank=0; irank!=MPISIZE/2; ++irank ) {         //  Loop over ranks in each block
          int idata = (irank / npart) * 2 * npart + irank % npart + i * npart;// Original index
          int isend = i * MPISIZE / 2 + irank;                  //   Compressed index
          irev[idata] = isend;                                  //   Map original to compressed index
          scntd[isend] = sendBodyCnt[idata];                    //   Permuted send count
          scnt2[i] += sendBodyCnt[idata] * bytes;               //   Send count per block
          for( int id=sendBodyDsp[idata]; id!=sendBodyDsp[idata]+sendBodyCnt[idata]; ++id,++ic ) {// Loop over bodies
            sendBuffer[ic] = sendBodies[id];                    //    Fill send buffer
          }                                                     //   End loop over bodies
        }                                                       //  End loop over ranks in each block
      }                                                         // End loop over blocks
      MPI_Alltoall(scntd,MPISIZE/2,MPI_INT,rcntd,MPISIZE/2,MPI_INT,MPI_COMM[l+1][2]);// Comm permuted count
      MPI_Alltoall(scnt2,1,MPI_INT,rcnt2,1,MPI_INT,MPI_COMM[l+1][2]);// Comm send count per block
      sdsp2[0] = 0; sdsp2[1] = scnt2[0];                        // Set send displacement
      rdsp2[0] = 0; rdsp2[1] = rcnt2[0];                        // Set recv displacement
      int rsize = (rdsp2[1] + rcnt2[1]) / bytes;                // Size of recieved bodies
      sendBodies.resize(rsize);                                 // Resize send bodies
      sendBuffer.resize(rsize);                                 // Resize send buffer
      recvBuffer.resize(rsize);                                 // Resize recv buffer
      MPI_Alltoallv(&sendBuffer[0],scnt2,sdsp2,MPI_BYTE,        // Communicate cells
                    &recvBuffer[0],rcnt2,rdsp2,MPI_BYTE,MPI_COMM[l+1][2]);// MPI_COMM[2] is for the one-to-one pair
      rdspd[0] = 0;                                             // Initialize permuted recv displacement
      for( int irank=0; irank!=MPISIZE-1; ++irank ) {           // Loop over ranks
        rdspd[irank+1] = rdspd[irank] + rcntd[irank];           //  Set permuted recv displacement
      }                                                         // End loop over ranks
      ic = 0;                                                   // Initiaize counter
      for( int i=0; i!=2; ++i ) {                               // Loop over the two blocks
        for( int irank=0; irank!=MPISIZE/2; ++irank ) {         //  Loop over ranks in each block
          int idata = (irank / npart) * 2 * npart + irank % npart + i * npart;// Original index
          int irecv = i * MPISIZE / 2 + irank;                  //   Compressed index
          recvBodyCnt[idata] = rcntd[irecv];                    //   Set recv cound
          idata = irev[irecv];                                  //   Set premuted index
          for( int id=rdspd[idata]; id!=rdspd[idata]+rcntd[idata]; ++id,++ic ) {// Loop over bodies
            sendBodies[ic] = recvBuffer[id];                    //    Get data from recv buffer
          }                                                     //   End loop over bodies
        }                                                       //  End loop over ranks in each block
      }                                                         // End loop over blocks
      recvBodyDsp[0] = 0;                                       // Initialize recv displacement
      for( int irank=0; irank!=MPISIZE-1; ++irank ) {           // Loop over ranks
        recvBodyDsp[irank+1] = recvBodyDsp[irank] + recvBodyCnt[irank];//  Set recv displacement
      }                                                         // End loop over ranks
      for( int irank=0; irank!=MPISIZE; ++irank ) {             // Loop over ranks
        sendBodyCnt[irank] = recvBodyCnt[irank];                //  Get next send count
        sendBodyDsp[irank] = recvBodyDsp[irank];                //  Get next send displacement
      }                                                         // End loop over ranks
    }                                                           // End loop over levels of N-D hypercube communication
    recvBodies = sendBodies;                                    // Copy send bodies to recv bodies
    delete[] scntd;                                             // Delete permuted send count
    delete[] rcntd;                                             // Delete permuted recv count
    delete[] rdspd;                                             // Delete permuted recv displacement
    delete[] irev;                                              // Delete map from original to compressed index
  }

//! Get boundries of domains on other processes
  void getOtherDomain(vect &xmin, vect &xmax, int l) {
    startTimer("Get domain   ");                                // Start timer
    MPI_Datatype MPI_TYPE = getType(XMIN[l][0]);                // Get MPI data type
    vect send[2],recv[2];                                       // Send and recv buffer
    MPI_Request req;                                            // MPI requests
    send[0] = send[1] = XMIN[l];                                // Set XMIN into send buffer
    recv[0] = recv[1] = 0;                                      // Initialize recv buffer
    MPI_Alltoall(send,3,MPI_TYPE,recv,3,MPI_TYPE,MPI_COMM[l][2]);// Communicate XMIN
    xmin = recv[1-key[l][2]];                                   // Copy xmin from recv buffer
    send[0] = send[1] = XMAX[l];                                // Set XMAX into send buffer
    recv[0] = recv[1] = 0;                                      // Initialize recv buffer
    MPI_Alltoall(send,3,MPI_TYPE,recv,3,MPI_TYPE,MPI_COMM[l][2]);// Communicate XMAX
    xmax = recv[1-key[l][2]];                                   // Copy xmax from recv buffer
    if( nprocs[l-1][0] % 2 == 1 && nprocs[l][0] >= nprocs[l][1] ) {// If right half of odd group
      int isend = (key[l][0] + 1            ) % nprocs[l][0];   //  Send to next rank (wrapped)
      int irecv = (key[l][0] - 1 + nprocs[l][0]) % nprocs[l][0];//  Recv from previous rank (wrapped)
      send[0] = xmin;                                           //  Set xmin in send buffer
      MPI_Isend(send,3,MPI_TYPE,isend,0,MPI_COMM[l][0],&req);   //  Send to next rank
      MPI_Irecv(recv,3,MPI_TYPE,irecv,0,MPI_COMM[l][0],&req);   //  Recv from previous rank
      MPI_Wait(&req,MPI_STATUS_IGNORE);                         //  Wait for recv to finish
      if( color[l][0] != color[l][1] ) xmin = recv[0];          //  Update only if leftover process of odd group
      send[0] = xmax;                                           //  Set xmax in send buffer
      MPI_Isend(send,3,MPI_TYPE,isend,0,MPI_COMM[l][0],&req);   //  Send to next rank
      MPI_Irecv(recv,3,MPI_TYPE,irecv,0,MPI_COMM[l][0],&req);   //  Recv from previous rank
      MPI_Wait(&req,MPI_STATUS_IGNORE);                         //  Wait for recv to finish
      if( color[l][0] != color[l][1] ) xmax = recv[0];          //  Update only if leftover process of odd group
    }                                                           // Endif for right half of odd group
    if( nprocs[l-1][0] == 1 ) {                                 // If previous group had one process
      xmin = XMIN[l];                                           //  Use own XMIN value for xmin
      xmax = XMAX[l];                                           //  Use own XMAX value for xmax
    }                                                           // Endif for isolated process
    stopTimer("Get domain   ",printNow);                        // Stop timer 
  }

//! Get disatnce to other domain
  real getDistance(C_iter C, vect xmin, vect xmax) {
    vect dist;                                                  // Distance vector
    for( int d=0; d!=3; ++d ) {                                 // Loop over dimensions
      dist[d] = (C->X[d] + Xperiodic[d] > xmax[d])*             //  Calculate the distance between cell C and
                (C->X[d] + Xperiodic[d] - xmax[d])+             //  the nearest point in domain [xmin,xmax]^3
                (C->X[d] + Xperiodic[d] < xmin[d])*             //  Take the differnece from xmin or xmax
                (C->X[d] + Xperiodic[d] - xmin[d]);             //  or 0 if between xmin and xmax
    }                                                           // End loop over dimensions
    real R = std::sqrt(norm(dist));                             // Scalar distance
    return R;
  }

//! Determine which cells to send
  void getLET(C_iter C0, C_iter C, vect xmin, vect xmax) {
    int level = int(log(MPISIZE-1) / M_LN2 / 3) + 1;            // Level of local root cell
    if( MPISIZE == 1 ) level = 0;                               // Account for serial case
    for( int i=0; i!=C->NCHILD; i++ ) {                         // Loop over child cells
      C_iter CC = C0+C->CHILD+i;                                //  Iterator for child cell
      bool divide = false;                                      //  Initialize logical for dividing
      if( IMAGES == 0 ) {                                       //  If free boundary condition
        Xperiodic = 0;                                          //   Set periodic coordinate offset
        real R = getDistance(CC,xmin,xmax);                     //   Get distance to other domain
        divide |= CLET * CC->R > THETA * R - EPS2;              //   If the cell seems too close and not twig
      } else {                                                  //  If periodic boundary condition
        for( int ix=-1; ix<=1; ++ix ) {                         //   Loop over x periodic direction
          for( int iy=-1; iy<=1; ++iy ) {                       //    Loop over y periodic direction
            for( int iz=-1; iz<=1; ++iz ) {                     //     Loop over z periodic direction
              Xperiodic[0] = ix * 2 * R0;                       //      Coordinate offset for x periodic direction
              Xperiodic[1] = iy * 2 * R0;                       //      Coordinate offset for y periodic direction
              Xperiodic[2] = iz * 2 * R0;                       //      Coordinate offset for z periodic direction
              real R = getDistance(CC,xmin,xmax);               //      Get distance to other domain
              divide |= CLET * CC->R > THETA * R - EPS2;        //      If the cell seems too close and not twig
            }                                                   //     End loop over z periodic direction
          }                                                     //    End loop over y periodic direction
        }                                                       //   End loop over x periodic direction
      }                                                         //  Endif for periodic boundary condition
      divide |= R0 / (1 << level) + 1e-5 < CC->R;               //  If the cell is larger than the local root cell
      if( divide && CC->NCHILD != 0 ) {                         //  If the cell seems too close and not twig
        getLET(C0,CC,xmin,xmax);                                //   Traverse the tree further
      } else {                                                  //  If the cell is far or a twig
        assert( R0 / (1 << level) + 1e-5 > CC->R );             //   Can't send cells that are larger than local root
        JCell cell;                                             //   Set compact cell type for sending
        cell.ICELL = CC->ICELL;                                 //   Set index of compact cell type
        cell.M     = CC->M;                                     //   Set Multipoles of compact cell type
        sendCells.push_back(cell);                              //    Push cell into send buffer vector
      }                                                         //  Endif for interaction
    }                                                           // End loop over child cells
    if( C->ICELL == 0 && C->NCHILD == 0 ) {                     // If the root cell has no children
      JCell cell;                                               //  Set compact cell type for sending
      cell.ICELL = C->ICELL;                                    //  Set index of compact cell type
      cell.M     = C->M;                                        //  Set Multipoles of compact cell type
      sendCells.push_back(cell);                                //  Push cell into send buffer vector
    }                                                           // Endif for root cells children
  }

//! Communicate cells by one-to-one MPI_Alltoallv
  void commCellsAlltoall(int l) {
    const int bytes = sizeof(sendCells[0]);                     // Byte size of JCell structure
    int rcnt[2], scnt[2] = {0, 0};                              // Recv count, send count
    scnt[1-key[l+1][2]] = sendCells.size()*bytes;               // Set send count to size of send buffer * bytes
    MPI_Alltoall(scnt,1,MPI_INT,rcnt,1,MPI_INT,MPI_COMM[l+1][2]);// Communicate the send count to get recv count
    int sdsp[2] = {0, scnt[0]};                                 // Send displacement
    int rdsp[2] = {0, rcnt[0]};                                 // Recv displacement
    if( color[l+1][0] != color[l+1][1] ) {                      // If leftover process of odd group
      rcnt[1-key[l+1][2]] = 0;                                  //  It won't have a pair for this communication
    }                                                           // Endif for leftover process of odd group
    recvCells.resize(rcnt[1-key[l+1][2]]/bytes);                // Resize recv buffer
    MPI_Alltoallv(&sendCells[0],scnt,sdsp,MPI_BYTE,             // Communicate cells
                  &recvCells[0],rcnt,rdsp,MPI_BYTE,MPI_COMM[l+1][2]);// MPI_COMM[2] is for the one-to-one pair
  }

//! Communicate cells by scattering from leftover processes
  void commCellsScatter(int l) {
    const int bytes = sizeof(sendCells[0]);                     // Byte size of JCell structure
    int numScatter = nprocs[l+1][1] - 1;                        // Number of processes to scatter to
    int oldSize = recvCells.size();                             // Size of recv buffer before communication
    int *scnt = new int [nprocs[l+1][1]];                       // Send count
    int *sdsp = new int [nprocs[l+1][1]];                       // Send displacement
    int rcnt;                                                   // Recv count
    if( key[l+1][1] == numScatter ) {                           // If this is the leftover proc to scatter from
      sdsp[0] = 0;                                              //  Initialize send displacement
      for(int i=0; i!=numScatter; ++i ) {                       //  Loop over processes to send to
        int begin = 0, end = sendCells.size();                  //   Set begin and end of range to send
        splitRange(begin,end,i,numScatter);                     //   Split range into numScatter and get i-th range
        scnt[i] = end - begin;                                  //   Set send count based on range
        sdsp[i+1] = sdsp[i] + scnt[i];                          //   Set send displacement based on send count
      }                                                         //  End loop over processes to send to
      scnt[numScatter] = 0;                                     //  Send count to self should be 0
    }                                                           // Endif for leftover proc
    MPI_Scatter(scnt,1,MPI_INT,&rcnt,1,MPI_INT,numScatter,MPI_COMM[l+1][1]);// Scatter the send count to get recv count

    recvCells.resize(oldSize+rcnt);                             // Resize recv buffer based on recv count
    for(int i=0; i!= nprocs[l+1][1]; ++i ) {                    // Loop over group of processes
      scnt[i] *= bytes;                                         //  Multiply send count by byte size of data
      sdsp[i] *= bytes;                                         //  Multiply send displacement by byte size of data
    }                                                           // End loop over group of processes
    rcnt *= bytes;                                              // Multiply recv count by byte size of data
    MPI_Scatterv(&sendCells[0],      scnt,sdsp,MPI_BYTE,        // Communicate cells via MPI_Scatterv
                 &recvCells[oldSize],rcnt,     MPI_BYTE,        // Offset recv buffer by oldSize
                 numScatter,MPI_COMM[l+1][1]);
    delete[] scnt;                                              // Delete send count
    delete[] sdsp;                                              // Delete send displacement
  }

//! Turn recv bodies to twigs
  void rbodies2twigs(Bodies &bodies, Cells &twigs) {
    startTimer("Recv bodies  ");                                //  Start timer
    for( JB_iter JB=recvBodies.begin(); JB!=recvBodies.end(); ++JB ) {// Loop over recv bodies
      Body body;                                                //  Body structure
      body.IBODY = 0;                                           //  Initialize body index
      body.IPROC = 0;                                           //  Initialize proc index
      body.TRG   = 0;                                           //  Initialize target values
      body.ICELL = JB->ICELL;                                   //  Set index of cell
      body.X     = JB->X;                                       //  Set position of body
      body.SRC   = JB->SRC;                                     //  Set source values of body
      bodies.push_back(body);                                   //  Push body into bodies vector
    }                                                           // End loop over recv bodies
    buffer.resize(bodies.size());                               // Resize sort buffer
    stopTimer("Recv bodies  ",printNow);                        //  Stop timer 
    sortBodies(bodies,buffer,false);                            // Sort bodies in descending order
    bodies2twigs(bodies,twigs);                                 // Turn bodies to twigs
  }

//! Turn cells to twigs
  void cells2twigs(Cells &cells, Cells &twigs, bool last) {
    while( !cells.empty() ) {                                   // While cell vector is not empty
      if( cells.back().NCHILD == 0 ) {                          //  If cell has no child
        if( cells.back().NDLEAF == 0 || !last ) {               //   If cell has no leaf or is not last iteration
          cells.back().NDLEAF = 0;                              //    Set number of leafs to 0
          twigs.push_back(cells.back());                        //    Push cell into twig vector
        }                                                       //   Endif for no leaf
      }                                                         //  Endif for no child
      cells.pop_back();                                         //  Pop last element from cell vector
    }                                                           // End while for cell vector
  }

//! Turn send buffer to twigs
  void send2twigs(Bodies &bodies, Cells &twigs, int offTwigs) {
    for( JC_iter JC=sendCells.begin(); JC!=sendCells.begin()+offTwigs; ++JC ) {// Loop over send buffer
      Cell cell;                                                //  Cell structure
      cell.ICELL = JC->ICELL;                                   //  Set index of cell
      cell.M     = JC->M;                                       //  Set multipole of cell
      cell.NDLEAF = cell.NCHILD = 0;                            //  Set number of leafs and children
      cell.LEAF  = bodies.end();                                //  Set pointer to first leaf
      getCenter(cell);                                          //  Set center and radius
      twigs.push_back(cell);                                    //  Push cell into twig vector
    }                                                           // End loop over send buffer
    sendCells.clear();                                          // Clear send buffer
  }

//! Turn recv buffer to twigs
  void recv2twigs(Bodies &bodies, Cells &twigs) {
    for( JC_iter JC=recvCells.begin(); JC!=recvCells.end(); ++JC ) {// Loop over recv buffer
      Cell cell;                                                //  Cell structure
      cell.ICELL = JC->ICELL;                                   //  Set index of cell
      cell.M     = JC->M;                                       //  Set multipole of cell
      cell.NDLEAF = cell.NCHILD = 0;                            //  Set number of leafs and children
      cell.LEAF  = bodies.end();                                //  Set pointer to first leaf
      getCenter(cell);                                          //  Set center and radius
      twigs.push_back(cell);                                    //  Push cell into twig vector
    }                                                           // End loop over recv buffer
  }

//! Zip two groups of twigs that overlap
  void zipTwigs(Cells &twigs, Cells &cells, Cells &sticks, bool last) {
    startTimer("Sort resize  ");                                // Start timer
    Cells cbuffer = twigs;                                      // Sort buffer for cells
    stopTimer("Sort resize  ",printNow);                        // Stop timer 
    sortCells(twigs,cbuffer);                                   // Sort twigs in ascending order
    startTimer("Ziptwigs     ");                                // Start timer
    bigint index = -1;                                          // Initialize index counter
    while( !twigs.empty() ) {                                   // While twig vector is not empty
      if( twigs.back().ICELL != index ) {                       //  If twig's index is different from previous
        cells.push_back(twigs.back());                          //   Push twig into cell vector
        index = twigs.back().ICELL;                             //   Update index counter
      } else if ( twigs.back().NDLEAF == 0 || !last ) {         //  Elseif twig-twig collision
        cells.back().M += twigs.back().M;                       //   Accumulate the multipole
      } else if ( cells.back().NDLEAF == 0 ) {                  //  Elseif twig-body collision
        Mset M;                                                 //   Multipole for temporary storage
        M = cells.back().M;                                     //   Save multipoles from cells
        cells.back() = twigs.back();                            //   Copy twigs to cells
        cells.back().M = M;                                     //   Copy back multipoles to cells
        twigs.back().M = M - twigs.back().M;                    //   Take the difference of the two
        if( std::abs(twigs.back().M[0]/M[0]) > EPS ) {          //   If the difference is non-zero
          sticks.push_back(twigs.back());                       //    Save this difference in the sticks vector
        }                                                       //   Endif for non-zero difference
      } else {                                                  //  Else body-body collision (don't do anything)
      }                                                         //  Endif for collision type
      twigs.pop_back();                                         //  Pop last element from twig vector
    }                                                           // End while for twig vector
    stopTimer("Ziptwigs     ",printNow);                        // Stop timer 
    sortCells(cells,cbuffer);                                   // Sort cells in ascending order
    startTimer("Ziptwigs     ");                                // Start timer
    twigs = cells;                                              // Copy cells to twigs
    cells.clear();                                              // Clear cells
    stopTimer("Ziptwigs     ",printNow);                        // Stop timer 
  }

//! Re-index bodies
  void reindexBodies(Bodies &bodies, Cells &twigs, Cells &cells ,Cells &sticks) {
    startTimer("Reindex      ");                                // Start timer
    while( !twigs.empty() ) {                                   // While twig vector is not empty
      if( twigs.back().NDLEAF == 0 ) {                          //  If twig has no leafs
        cells.push_back(twigs.back());                          //   Push twig into cell vector
      }                                                         //  Endif for no leafs
      twigs.pop_back();                                         //  Pop last element from twig vector
    }                                                           // End while for twig vector
//    BottomUp::setIndex(bodies,-1,0,0,true);                     // Set index of bodies
    buffer.resize(bodies.size());                               // Resize sort buffer
    stopTimer("Reindex      ",printNow);                        // Stop timer 
//    sortBodies(bodies,buffer,false);                            // Sort bodies in descending order
//    BottomUp::grow(bodies);                                     // Grow tree structure
    sortBodies(bodies,buffer,false);                              // Sort bodies in descending order
    bodies2twigs(bodies,twigs);                                 // Turn bodies to twigs
    startTimer("Reindex      ");                                // Start timer
    for( C_iter C=twigs.begin(); C!=twigs.end(); ++C ) {        // Loop over cells
      if( sticks.size() > 0 ) {                                 //  If stick vector is not empty
        if( C->ICELL == sticks.back().ICELL ) {                 //   If twig's index is equal to stick's index
          C->M += sticks.back().M;                              //    Accumulate multipole
          sticks.pop_back();                                    //    Pop last element from stick vector
        }                                                       //   Endif for twig's index
      }                                                         //  Endif for stick vector
    }                                                           // End loop over cells
    cells.insert(cells.begin(),twigs.begin(),twigs.end());      // Add twigs to the end of cell vector
    cells.insert(cells.begin(),sticks.begin(),sticks.end());    // Add remaining sticks to the end of cell vector
    sticks.clear();                                             // Clear sticks
    Cells cbuffer = cells;                                      // Sort buffer for cells
    stopTimer("Reindex      ",printNow);                        // Stop timer 
    sortCells(cells,cbuffer);                                   // Sort cells in ascending order
    startTimer("Reindex      ");                                // Start timer
    twigs = cells;                                              // Copy cells to twigs
    cells.clear();                                              // Clear cells
    stopTimer("Reindex      ",printNow);                        // Stop timer 
  }

//! Turn sticks to send buffer
  void sticks2send(Cells &sticks, int &offTwigs) {
    while( !sticks.empty() ) {                                  // While stick vector is not empty
      JCell cell;                                               //  Cell structure
      cell.ICELL = sticks.back().ICELL;                         //  Set index of cell
      cell.M     = sticks.back().M;                             //  Set multipole of cell
      sendCells.push_back(cell);                                //  Push cell into send buffer
      sticks.pop_back();                                        //  Pop last element of stick vector
    }                                                           // End while for stick vector
    offTwigs = sendCells.size();                                // Keep track of current send buffer size
  }

//! Validate number of send cells
  void checkNumCells(int l) {                                   // Only works with octsection
    int maxLevel = int(log(MPISIZE-1) / M_LN2 / 3) + 1;
    if( MPISIZE == 1 ) maxLevel = 0;
    int octant0 = -1;
    int numCells = 0;
    for( JC_iter JC=sendCells.begin(); JC!=sendCells.end(); ++JC ) {
      int level = getLevel(JC->ICELL);
      int index = JC->ICELL - ((1 << 3*level) - 1) / 7;
      int octant = index / (1 << 3 * (level - maxLevel));
      if( octant != octant0 ) {
        octant0 = octant;
        numCells++;
      }
    }
    int numCellsExpect = (1 << (3 * maxLevel - 1)) / (1 << l);  // Isn't true for far domains
    if( numCellsExpect != numCells && MPIRANK == 0) std::cout << numCells << " " << numCellsExpect << std::endl;
  }

//! Check total charge
  void checkSumMass(Cells &cells) {
    real localMass = 0;
    for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {
      if( C->NCHILD == 0 ) {
        localMass += std::abs(C->M[0]);
      }
    }
    real globalMass;
    MPI_Allreduce(&localMass,&globalMass,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    print("localMass : ",0);
    print(localMass);
    print("globalMass : ",0);
    print(globalMass,0);
    print("\n",0);
  }

public:
//! Constructor
  ParallelFMM() : Partition<equation>() {}
//! Destructor
  ~ParallelFMM() {}

//! Set bodies to communicate
  void setCommBodies(Cells &cells) {
    startTimer("Gather bounds");                                // Start timer
    gatherBounds();                                             // Gather bounds of other domain
    stopTimer("Gather bounds",printNow);                        // Stop timer 
    startTimer("Get send rank");                                // Start timer
    getSendRank(cells);                                         // Get neighbor ranks to send to
    stopTimer("Get send rank",printNow);                        // Stop timer 
  }

//! Update bodies using the previous send count
  void updateBodies(bool comm=true) {
    startTimer("Get send cnt ");                                // Start timer
    getSendCount(comm);                                         // Get size of data to send
    stopTimer("Get send cnt ",printNow);                        // Stop timer 
    startTimer("Alltoall B   ");                                // Start timer
#if 1
    int bytes = sizeof(sendBodies[0]);                          // Byte size of jbody structure
    for( int i=0; i!=MPISIZE; ++i ) {                           // Loop over ranks
      sendBodyCnt[i] *= bytes;                                  //  Multiply by bytes
      sendBodyDsp[i] *= bytes;                                  //  Multiply by bytes
      recvBodyCnt[i] *= bytes;                                  //  Multiply by bytes
      recvBodyDsp[i] *= bytes;                                  //  Multiply by bytes
    }                                                           // End loop over ranks
    MPI_Alltoallv(&sendBodies[0],&sendBodyCnt[0],&sendBodyDsp[0],MPI_BYTE,
                  &recvBodies[0],&recvBodyCnt[0],&recvBodyDsp[0],MPI_BYTE,MPI_COMM_WORLD);
    for( int i=0; i!=MPISIZE; ++i ) {                           // Loop over ranks
      sendBodyCnt[i] /= bytes;                                  //  Divide by bytes
      sendBodyDsp[i] /= bytes;                                  //  Divide by bytes
      recvBodyCnt[i] /= bytes;                                  //  Divide by bytes
      recvBodyDsp[i] /= bytes;                                  //  Divide by bytes
    }                                                           // End loop over ranks
#else
    commBodiesAlltoall();
#endif
    sendBodies.clear();                                         // Clear send buffer for bodies
    stopTimer("Alltoall B   ",printNow);                        // Stop timer 
  }

//! Communicate bodies in the local essential tree
  void commBodies(Cells &cells) {
    setCommBodies(cells);                                       // Set bodies to communicate
    updateBodies();                                             // Update bodies with alltoall
  }

//! Convert recvBodies to cells
  void bodies2cells(Bodies &bodies, Cells &cells) {
    Cells twigs,sticks;                                         // Twigs and sticks are special types of cells
    rbodies2twigs(bodies,twigs);                                // Put recv bodies into twig vector
    twigs2cells(twigs,cells,sticks);                            // Turn twigs to cells
  }

//! Communicate cells in the local essential tree
  void commCells(Bodies &bodies, Cells &cells) {
    vect xmin = 0, xmax = 0;                                    // Initialize domain boundaries
    Cells twigs,sticks;                                         // Twigs and sticks are special types of cells

#if 1
    startTimer("Get LET      ");                                // Start timer
    int ssize = 0;                                              // Initialize offset for send cells
    sendCellCnt.assign(MPISIZE,0);                              // Initialize cell send count
    sendCellDsp.assign(MPISIZE,0);                              // Initialize cell send displacement
    for( int irank=0; irank!=MPISIZE; ++irank ) {               // Loop over ranks to send to
      getLET(cells.begin(),cells.end()-1,xminAll[irank],xmaxAll[irank]);//  Determine which cells to send
      sendCellCnt[irank] = sendCells.size()-ssize;              //  Set cell send count of current rank
      sendCellDsp[irank] = ssize;                               //  Set cell send displacement of current rank
      ssize += sendCellCnt[irank];                              //  Increment offset for vector send cells
    }                                                           // End loop over ranks
    stopTimer("Get LET      ",printNow);                        // Stop timer 
    startTimer("Alltoall C   ");                                // Start timer
    MPI_Alltoall(&sendCellCnt[0],1,MPI_INT,&recvCellCnt[0],1,MPI_INT,MPI_COMM_WORLD);// Communicate the send counts
    int rsize = 0;                                              // Initialize total recv count
    for( int irank=0; irank!=MPISIZE; ++irank ) {               // Loop over ranks to recv from
      recvCellDsp[irank] = rsize;                               //  Set recv displacements
      rsize += recvCellCnt[irank];                              //  Accumulate recv counts
    }                                                           // End loop over ranks to recv from
    recvCells.resize(rsize);                                    // Resize recv buffer
    int bytes = sizeof(sendCells[0]);                           // Byte size of jbody structure
    for( int i=0; i!=MPISIZE; ++i ) {                           // Loop over ranks
      sendCellCnt[i] *= bytes;                                  //  Multiply by bytes
      sendCellDsp[i] *= bytes;                                  //  Multiply by bytes
      recvCellCnt[i] *= bytes;                                  //  Multiply by bytes
      recvCellDsp[i] *= bytes;                                  //  Multiply by bytes
    }                                                           // End loop over ranks
    MPI_Alltoallv(&sendCells[0],&sendCellCnt[0],&sendCellDsp[0],MPI_BYTE,
                  &recvCells[0],&recvCellCnt[0],&recvCellDsp[0],MPI_BYTE,MPI_COMM_WORLD);
    for( int i=0; i!=MPISIZE; ++i ) {                           // Loop over ranks
      sendCellCnt[i] /= bytes;                                  //  Divide by bytes
      sendCellDsp[i] /= bytes;                                  //  Divide by bytes
      recvCellCnt[i] /= bytes;                                  //  Divide by bytes
      recvCellDsp[i] /= bytes;                                  //  Divide by bytes
    }                                                           // End loop over ranks
    stopTimer("Alltoall C   ",printNow);                        // Stop timer 
    rbodies2twigs(bodies,twigs);                                // Put recv bodies into twig vector
    startTimer("Cells2twigs  ");                                // Start timer
    cells2twigs(cells,twigs,true);                              // Put cells into twig vector
    stopTimer("Cells2twigs  ",printNow);                        // Stop timer 
    startTimer("Recv2twigs   ");                                // Start timer
    recv2twigs(bodies,twigs);                                   // Put recv buffer into twig vector
    stopTimer("Recv2twigs   ",printNow);                        // Stop timer 
    zipTwigs(twigs,cells,sticks,true);                          // Zip two groups of twigs that overlap
    reindexBodies(bodies,twigs,cells,sticks);                   // Re-index bodies
    twigs2cells(twigs,cells,sticks);                            // Turn twigs to cells
    sendCells.clear();                                          // Clear send buffer
    recvCells.clear();                                          // Clear recv buffer
#else
    int offTwigs = 0;                                           // Initialize offset of twigs
    for( int l=0; l!=LEVEL; ++l ) {                             // Loop over levels of N-D hypercube communication
      getOtherDomain(xmin,xmax,l+1);                            //  Get boundries of domains on other processes
      startTimer("Get LET      ");                              //  Start timer
      getLET(cells.begin(),cells.end()-1,xmin,xmax);            //  Determine which cells to send
#ifdef DEBUG
      checkNumCells(LEVEL-l-1);
      checkSumMass(cells);
#endif
      stopTimer("Get LET      ",printNow);                      //  Stop timer 
      startTimer("Alltoall C   ");                              //  Start timer
      commCellsAlltoall(l);                                     //  Communicate cells by one-to-one MPI_Alltoallv
      if( nprocs[l][0] % 2 == 1 && nprocs[l][0] != 1 && nprocs[l+1][0] <= nprocs[l+1][1] ) {// If scatter is necessary
        commCellsScatter(l);                                    //   Communicate cells by scattering from leftover proc
      }                                                         //  Endif for odd number of procs
      stopTimer("Alltoall C   ",printNow);                      //  Stop timer 
      if( l == LEVEL - 1 ) rbodies2twigs(bodies,twigs);         //  Put recv bodies into twig vector
      startTimer("Cells2twigs  ");                              //  Start timer
      cells2twigs(cells,twigs,l==LEVEL-1);                      //  Put cells into twig vector
      stopTimer("Cells2twigs  ",printNow);                      //  Stop timer 
      startTimer("Send2twigs   ");                              //  Start timer
      send2twigs(bodies,twigs,offTwigs);                        //  Put send buffer (sticks) into twig vector
      stopTimer("Send2twigs   ",printNow);                      //  Stop timer 
      startTimer("Recv2twigs   ");                              //  Start timer
      recv2twigs(bodies,twigs);                                 //  Put recv buffer into twig vector
      stopTimer("Recv2twigs   ",printNow);                      //  Stop timer 
#ifdef DEBUG
      if( l == LEVEL - 1 ) {                                    //  If at last level
        complex SUM = 0;                                        //   Initialize accumulator
        for(C_iter C=twigs.begin(); C!=twigs.end(); ++C) {      //   Loop over twigs
          if( C->NDLEAF == 0 ) SUM += C->M[0];                  //    Add multipoles of empty twigs
        }                                                       //   End loop over twigs
        print("Before recv   : ",0);                            //   Print identifier
        print(SUM);                                             //   Print sum of multipoles
      }                                                         //  Endif for last level
#endif
      zipTwigs(twigs,cells,sticks,l==LEVEL-1);                  //  Zip two groups of twigs that overlap
#ifdef DEBUG
      if( l == LEVEL - 1 ) {                                    //  If at last level
        complex SUM = 0;                                        //   Initialize accumulator
        for(C_iter C=twigs.begin(); C!=twigs.end(); ++C) {      //   Loop over twigs
          SUM += C->M[0];                                       //    Add multipoles
        }                                                       //   End loop over twigs
        print("After merge   : ",0);                            //   Print identifier
        print(SUM);                                             //   Print sum of multipoles
        print("sticks.size() : ",0);                            //   Print identifier
        print(sticks.size());                                   //   Print size of stick vector
      }                                                         //  Endif for last level
#endif
      if( l == LEVEL - 1 ) reindexBodies(bodies,twigs,cells,sticks);// Re-index bodies
      twigs2cells(twigs,cells,sticks);                          //  Turn twigs to cells
      startTimer("Sticks2send  ");                              //  Start timer
      sticks2send(sticks,offTwigs);                             //  Turn sticks to send buffer
      stopTimer("Sticks2send  ",printNow);                      //  Stop timer 
    }                                                           // End loop over levels of N-D hypercube communication
#endif

#ifdef DEBUG
    print("M[0] @ root   : ",0);                                // Print identifier
    print((cells.end()-1)->M[0]);                               // Print monopole of root (should be 1 for test)
    print("bodies.size() : ",0);                                // Print identifier
    print(bodies.size());                                       // Print size of body vector
#endif
    sendCells.clear();                                          // Clear send buffer
    recvCells.clear();                                          // Clear recv buffer
  }

//! Remove cells that belong to current process
  void eraseLocalTree(Cells &cells) {
    int level = int(log(MPISIZE-1) / M_LN2 / 3) + 1;            // Level of process root cell
    if( MPISIZE == 1 ) level = 0;                               // Account for serial case
    int off = ((1 << 3 * level) - 1) / 7;                       // Levelwise offset of ICELL
    int size = (1 << 3 * level) / MPISIZE;                      // Number of cells to remove
    unsigned begin = MPIRANK * size + off;                      // Begin index of cells to remove
    unsigned end = (MPIRANK + 1) * size + off;                  // End index of cells to remove
    for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {        // Loop over cells
      int nchild = 0;                                           //  Initialize child cell counter
      for( int c=0; c!=C->NCHILD; ++c ) {                       //  Loop over child cells
        C_iter CC = cells.begin()+C->CHILD+c;                   //   Iterator of child cell
        if( CC->ICELL < begin || end <= CC->ICELL ) {           //   If child cell is not within the removal range
          C_iter CH = cells.begin()+C->CHILD+nchild;            //    New iterator of child cell
          *CH = *CC;                                            //    Copy data of child cell
          nchild++;                                             //    Increment child cell counter
        }                                                       //   Endif for removal range
      }                                                         //  End loop over child cells
      C->NCHILD = nchild;                                       //  Update number of child cells
    }                                                           // End loop over cells
  }
};

#endif
