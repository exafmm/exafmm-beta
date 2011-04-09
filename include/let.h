#ifndef let_h
#define let_h
#include "partition.h"

class LocalEssentialTree : public Partition {                   // Handles all the communication in this code
private:
  JBodies sendBodies;                                           // Send buffer for bodies
  JBodies recvBodies;                                           // Receive buffer for bodies
  JCells  sendCells;                                            // Send buffer for cells
  JCells  recvCells;                                            // Receive buffer for cells

private:
  void getOtherDomain(vect &xmin, vect &xmax, int l) {          // Get boundries of domains on other processes
    startTimer("Get domain   ");                                //  Start timer
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
    stopTimer("Get domain   ",printNow);                        //  Stop timer & print
  }

  real getDistance(C_iter C, vect xmin, vect xmax) {            // Get disatnce to other domain
    vect dist;                                                  // Distance vector
    for( int d=0; d!=3; ++d )                                   // Loop over dimensions
      dist[d] = (C->X[d] + Xperiodic[d] > xmax[d])*             //  Calculate the distance between cell C and
                (C->X[d] + Xperiodic[d] - xmax[d])+             //  the nearest point in domain [xmin,xmax]^3
                (C->X[d] + Xperiodic[d] < xmin[d])*             //  Take the differnece from xmin or xmax
                (C->X[d] + Xperiodic[d] - xmin[d]);             //  or 0 if between xmin and xmax
    real R = std::sqrt(norm(dist));                             // Scalar distance
    return R;
  }

  void getLET(C_iter C0, C_iter C, vect xmin, vect xmax) {      // Determine which cells to send
    int level = int(log(SIZE-1) / M_LN2 / 3) + 1;               // Level of local root cell
    for( int i=0; i!=C->NCHILD; i++ ) {                         // Loop over child cells
      C_iter CC = C0+C->CHILD[i];                               //  Iterator for child cell
      bool divide = false;                                      //  Initialize logical for dividing
      if( IMAGES == 0 ) {                                       //  If free boundary condition
        Xperiodic = 0;                                          //   Set periodic coordinate offset
        real R = getDistance(CC,xmin,xmax);                     //   Get distance to other domain
        divide |= CLET * CC->R > THETA * R;                     //   If the cell seems too close and not twig
      } else {                                                  //  If periodic boundary condition
        for( int ix=-1; ix<=1; ++ix ) {                         //   Loop over x periodic direction
          for( int iy=-1; iy<=1; ++iy ) {                       //    Loop over y periodic direction
            for( int iz=-1; iz<=1; ++iz ) {                     //     Loop over z periodic direction
              Xperiodic[0] = ix * 2 * R0;                       //      Coordinate offset for x periodic direction
              Xperiodic[1] = iy * 2 * R0;                       //      Coordinate offset for y periodic direction
              Xperiodic[2] = iz * 2 * R0;                       //      Coordinate offset for z periodic direction
              real R = getDistance(CC,xmin,xmax);               //      Get distance to other domain
              divide |= CLET * CC->R > THETA * R;               //      If the cell seems too close and not twig
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

  void commCellsAlltoall(int l) {                               // Communicate cells by one-to-one MPI_Alltoallv
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

  void commCellsScatter(int l) {                                // Communicate cells by scattering from leftover proc
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

  void rbodies2twigs(Bodies &bodies, Cells &twigs) {            // Turn recv bodies to twigs
    for( JB_iter JB=recvBodies.begin(); JB!=recvBodies.end(); ++JB ) {// Loop over recv bodies
      Body body;                                                //  Body structure
      body.ICELL = JB->ICELL;                                   //  Set index of body
      body.X     = JB->X;                                       //  Set position of body
      body.SRC   = JB->SRC;                                     //  Set source values of body
      bodies.push_back(body);                                   //  Push body into bodies vector
    }                                                           // End loop over recv bodies
    buffer.resize(bodies.size());                               // Resize sort buffer
    sortBodies(bodies,buffer,false);                            // Sort bodies in descending order
    bodies2twigs(bodies,twigs);                                 // Turn bodies to twigs
  }

  void cells2twigs(Cells &cells, Cells &twigs, bool last) {     // Turn cells to twigs
    while( !cells.empty() ) {                                   // While cell vector is not empty
      if( cells.back().NCHILD == 0 ) {                          //  If cell has no child
        if( cells.back().NLEAF == 0 || !last ) {                //   If cell has no leaf or is not last iteration
          cells.back().NLEAF = 0;                               //    Set number of leafs to 0
          twigs.push_back(cells.back());                        //    Push cell into twig vector
        }                                                       //   Endif for no leaf
      }                                                         //  Endif for no child
      cells.pop_back();                                         //  Pop last element from cell vector
    }                                                           // End while for cell vector
  }

  void send2twigs(Bodies &bodies, Cells &twigs, int offTwigs) { // Turn send buffer to twigs
    for( JC_iter JC=sendCells.begin(); JC!=sendCells.begin()+offTwigs; ++JC ) {// Loop over send buffer
      Cell cell;                                                //  Cell structure
      cell.ICELL = JC->ICELL;                                   //  Set index of cell
      cell.M     = JC->M;                                       //  Set multipole of cell
      cell.NLEAF = cell.NCHILD = 0;                             //  Set number of leafs and children
      cell.LEAF  = bodies.end();                                //  Set pointer to first leaf
      getCenter(cell);                                          //  Set center and radius
      twigs.push_back(cell);                                    //  Push cell into twig vector
    }                                                           // End loop over send buffer
    sendCells.clear();                                          // Clear send buffer
  }

  void recv2twigs(Bodies &bodies, Cells &twigs) {               // Turn recv buffer to twigs
    for( JC_iter JC=recvCells.begin(); JC!=recvCells.end(); ++JC ) {// Loop over recv buffer
      Cell cell;                                                //  Cell structure
      cell.ICELL = JC->ICELL;                                   //  Set index of cell
      cell.M     = JC->M;                                       //  Set multipole of cell
      cell.NLEAF = cell.NCHILD = 0;                             //  Set number of leafs and children
      cell.LEAF  = bodies.end();                                //  Set pointer to first leaf
      getCenter(cell);                                          //  Set center and radius
      twigs.push_back(cell);                                    //  Push cell into twig vector
    }                                                           // End loop over recv buffer
  }

  void zipTwigs(Cells &twigs, Cells &cells, Cells &sticks, bool last) {// Zip two groups of twigs that overlap
    sortCells(twigs);                                           // Sort twigs in ascending order
    bigint index = -1;                                          // Initialize index counter
    while( !twigs.empty() ) {                                   // While twig vector is not empty
      if( twigs.back().ICELL != index ) {                       //  If twig's index is different from previous
        cells.push_back(twigs.back());                          //   Push twig into cell vector
        index = twigs.back().ICELL;                             //   Update index counter
      } else if ( twigs.back().NLEAF == 0 || !last ) {          //  Elseif twig-twig collision
        cells.back().M += twigs.back().M;                       //   Accumulate the multipole
      } else if ( cells.back().NLEAF == 0 ) {                   //  Elseif twig-body collision
        coef M;                                                 //   Multipole for temporary storage
        M = cells.back().M;                                     //   Save multipoles from cells
        cells.back() = twigs.back();                            //   Copy twigs to cells
        cells.back().M = M;                                     //   Copy back multipoles to cells
        twigs.back().M = M - twigs.back().M;                    //   Take the difference of the two
        if( std::abs(twigs.back().M[0]/M[0]) > 1e-6 ) {         //   If the difference is non-zero
          sticks.push_back(twigs.back());                       //    Save this difference in the sticks vector
        }                                                       //   Endif for non-zero difference
      } else {                                                  //  Else body-body collision (don't do anything)
      }                                                         //  Endif for collision type
      twigs.pop_back();                                         //  Pop last element from twig vector
    }                                                           // End while for twig vector
    sortCells(cells);                                           // Sort cells in ascending order
    twigs = cells;                                              // Copy cells to twigs
    cells.clear();                                              // Clear cells
  }

  void reindexBodies(Bodies &bodies, Cells &twigs, Cells &cells ,Cells &sticks) {// Re-index bodies
    while( !twigs.empty() ) {                                   // While twig vector is not empty
      if( twigs.back().NLEAF == 0 ) {                           //  If twig has no leafs
        cells.push_back(twigs.back());                          //   Push twig into cell vector
      }                                                         //  Endif for no leafs
      twigs.pop_back();                                         //  Pop last element from twig vector
    }                                                           // End while for twig vector
    BottomUp::setIndex(bodies,-1,0,0,true);                     // Set index of bodies
    buffer.resize(bodies.size());                               // Resize sort buffer
//    sortBodies(bodies,buffer);                                  // Sort bodies in ascending order
//    BottomUp::grow(bodies);                                     // Grow tree structure
    sortBodies(bodies,buffer);                                  // Sort bodies in ascending order
    bodies2twigs(bodies,twigs);                                 // Turn bodies to twigs
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
    sortCells(cells);                                           // Sort cells in ascending order
    twigs = cells;                                              // Copy cells to twigs
    cells.clear();                                              // Clear cells
  }

  void sticks2send(Cells &sticks, int &offTwigs) {              // Turn sticks to send buffer
    while( !sticks.empty() ) {                                  // While stick vector is not empty
      JCell cell;                                               //  Cell structure
      cell.ICELL = sticks.back().ICELL;                         //  Set index of cell
      cell.M     = sticks.back().M;                             //  Set multipole of cell
      sendCells.push_back(cell);                                //  Push cell into send buffer
      sticks.pop_back();                                        //  Pop last element of stick vector
    }                                                           // End while for stick vector
    offTwigs = sendCells.size();                                // Keep track of current send buffer size
  }

public:
  LocalEssentialTree() : Partition() {}                         // Constructor
  ~LocalEssentialTree() {}                                      // Destructor

  void commBodies(Cells &cells) {                               // Communicate bodies in the LET
    startTimer("Gather bounds");                                // Start timer
    MPI_Datatype MPI_TYPE = getType(XMIN[LEVEL][0]);            // Get MPI data type
    std::vector<vect> xmin(SIZE);                               // Buffer for gathering XMIN
    std::vector<vect> xmax(SIZE);                               // Buffer for gathering XMAX
    MPI_Allgather(&XMIN[LEVEL][0],3,MPI_TYPE,                   // Gather XMIN
                  &xmin[0][0],3,MPI_TYPE,MPI_COMM_WORLD);
    MPI_Allgather(&XMAX[LEVEL][0],3,MPI_TYPE,                   // Gather XMAX
                  &xmax[0][0],3,MPI_TYPE,MPI_COMM_WORLD);
    stopTimer("Gather bounds",printNow);                        // Stop timer & print
    startTimer("Get send rank");                                // Start timer
    std::vector<int> srnks,scnts;                               // Ranks to send to, and their send counts
    std::vector<C_iter> scells;                                 // Vector of cell iterators for cells to send
    int oldsize = 0;                                            // Per rank offset of the number of cells to send
    C_iter C0 = cells.begin();                                  // Set cell begin iterator
    for( int irank=0; irank!=SIZE; ++irank ) {                  // Loop over ranks
      int ic = 0;                                               //  Initialize neighbor dimension counter
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimensions
        if(xmin[irank][d] < 2 * XMAX[LEVEL][d] - XMIN[LEVEL][d] &&//   If the two domains are touching or overlapping
           xmax[irank][d] > 2 * XMIN[LEVEL][d] - XMAX[LEVEL][d]) {//   in all dimensions, they are neighboring domains
          ic++;                                                 //    Increment neighbor dimension counter
        }                                                       //   Endif for overlapping domains
      }                                                         //  End loop over dimensions
      ic = 3;
      if( ic == 3 && irank != RANK ) {                          //  If ranks are neighbors in all dimensions
        for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {    //   Loop over cells
          if( C->NCHILD == 0 ) {                                //   If cell is a twig
            bool send = false;                                  //    Initialize logical for sending
            if( IMAGES == 0 ) {                                 //    If free boundary condition
              Xperiodic = 0;                                    //     Set periodic coordinate offset
              real R = getDistance(C,xmin[irank],xmax[irank]);  //     Get distance to other domain
              send |= CLET * C->R > THETA * R;                  //      If the cell seems close enough for P2P
            } else {                                            //     If periodic boundary condition
              for( int ix=-1; ix<=1; ++ix ) {                   //      Loop over x periodic direction
                for( int iy=-1; iy<=1; ++iy ) {                 //       Loop over y periodic direction
                  for( int iz=-1; iz<=1; ++iz ) {               //        Loop over z periodic direction
                    Xperiodic[0] = ix * 2 * R0;                 //         Coordinate offset for x periodic direction
                    Xperiodic[1] = iy * 2 * R0;                 //         Coordinate offset for y periodic direction
                    Xperiodic[2] = iz * 2 * R0;                 //         Coordinate offset for z periodic direction
                    real R = getDistance(C,xmin[irank],xmax[irank]);//     Get distance to other domain
                    send |= CLET * C->R > THETA * R;            //         If the cell seems close enough for P2P
                  }                                             //        End loop over z periodic direction
                }                                               //       End loop over y periodic direction
              }                                                 //      End loop over x periodic direction
            }                                                   //     Endif for periodic boundary condition
            if( send ) {                                        //     If the cell seems close enough for P2P
              scells.push_back(C);                              //      Add cell iterator to scells
            }                                                   //     Endif for cell distance
          }                                                     //    Endif for twigs
        }                                                       //   End loop over cells
        srnks.push_back(irank);                                 //   Add current rank to srnks
        scnts.push_back(scells.size()-oldsize);                 //   Add current cell count to scnts
        oldsize = scells.size();                                //   Set new offset for cell count
      }                                                         //  Endif for neighbor ranks
    }                                                           // End loop over ranks
    stopTimer("Get send rank",printNow);                        // Stop timer & print

    startTimer("Send count   ");                                // Start timer
    int ic = 0, ssize = 0;                                      // Initialize counter and offset for scells
    std::vector<MPI_Request> reqs(2*srnks.size());              // Vector of MPI requests
    std::vector<int>         scnt(SIZE);                        // Vector of send counts
    std::vector<int>         sdsp(SIZE);                        // Vector of send displacements
    std::vector<int>         rcnt(SIZE);                        // Vector of recv counts
    std::vector<int>         rdsp(SIZE);                        // Vector of recv displacements
    for( int i=0; i!=int(srnks.size()); ++i ) {                 // Loop over ranks to send to & recv from
      int irank = srnks[i];                                     //  Rank to send to & recv from
      for( int c=0; c!=scnts[i]; ++c,++ic ) {                   //  Loop over cells to send to that rank
        C_iter C = scells[ic];                                  //   Set cell iterator
        for( B_iter B=C->LEAF; B!=C->LEAF+C->NLEAF; ++B ) {     //   Loop over bodies in that cell
          JBody body;                                           //    Set compact body type for sending
          body.ICELL = B->ICELL;                                //    Set cell index of compact body type
          body.X     = B->X;                                    //    Set position of compact body type
          body.SRC   = B->SRC;                                  //    Set source values of compact body type
          sendBodies.push_back(body);                           //    Push it into the send buffer
        }                                                       //   End loop over bodies
      }                                                         //  End loop over cells
      scnt[irank] = sendBodies.size()-ssize;                    //  Set send count of current rank
      sdsp[irank] = ssize;
      ssize += scnt[irank];                                     //  Increment offset for vector scells
    }                                                           // End loop over ranks
    MPI_Alltoall(&scnt[0],1,MPI_INT,&rcnt[0],1,MPI_INT,MPI_COMM_WORLD);
    stopTimer("Send count   ",printNow);                        // Stop timer & print

    startTimer("Send bodies  ");                                // Start timer
    int rsize = 0;                                              // Initialize total recv count
    for( int i=0; i!=SIZE; ++i ) {                              // Loop over ranks to recv from
      rdsp[i] = rsize;                                          //  Set recv displacements
      rsize += rcnt[i];                                         //  Accumulate recv counts
    }                                                           // End loop over ranks to recv from
    recvBodies.resize(rsize);                                   // Receive buffer for bodies
    int bytes = sizeof(sendBodies[0]);                          // Byte size of jbody structure
    for( int i=0; i!=SIZE; ++i ) {                              // Loop over ranks to recv from
      scnt[i] *= bytes;                                         //  Send as bytes
      sdsp[i] *= bytes;                                         //  Send as bytes
      rcnt[i] *= bytes;                                         //  Recv as bytes
      rdsp[i] *= bytes;                                         //  Recv as bytes
    }                                                           // End loop over ranks to recv from
    MPI_Alltoallv(&sendBodies[0],&scnt[0],&sdsp[0],MPI_BYTE,
                  &recvBodies[0],&rcnt[0],&rdsp[0],MPI_BYTE,MPI_COMM_WORLD);
    sendBodies.clear();                                         // Clear send buffer for bodies
    stopTimer("Send bodies  ",printNow);                        // Stop timer & print
  }

  void checkNumCells(int l) {                                   // Only works with octsection
    int maxLevel = int(log(SIZE-1) / M_LN2 / 3) + 1;
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
    if( numCellsExpect != numCells && RANK == 0) std::cout << numCells << " " << numCellsExpect << std::endl;
  }

  void checkSumMass(Cells &cells) {
    double localMass = 0;
    for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {
      if( C->NCHILD == 0 ) {
        localMass += C->M[0].real();
      }
    }
    double globalMass;
    MPI_Allreduce(&localMass,&globalMass,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    print("localMass : ",0);
    print(localMass);
    print("globalMass : ",0);
    print(globalMass,0);
    print("\n",0);
  }

  void commCells(Bodies &bodies, Cells &cells) {                // Communicate cell in the LET
    int offTwigs = 0;                                           // Initialize offset of twigs
    vect xmin = 0, xmax = 0;                                    // Initialize domain boundaries
    Cells twigs,sticks;                                         // Twigs and sticks are special types of cells

    for( int l=0; l!=LEVEL; ++l ) {                             // Loop over levels of N-D hypercube communication
      getOtherDomain(xmin,xmax,l+1);                            //  Get boundries of domains on other processes
      startTimer("Get LET      ");                              //  Start timer
      getLET(cells.begin(),cells.end()-1,xmin,xmax);            //  Determine which cells to send
#ifdef DEBUG
      checkNumCells(LEVEL-l-1);
      checkSumMass(cells);
#endif
      stopTimer("Get LET      ",printNow);                      //  Stop timer & print
      startTimer("Alltoall     ");                              //  Start timer
      commCellsAlltoall(l);                                     //  Communicate cells by one-to-one MPI_Alltoallv
      if( nprocs[l][0] % 2 == 1 && nprocs[l][0] != 1 && nprocs[l+1][0] <= nprocs[l+1][1] ) {// If scatter is necessary
        commCellsScatter(l);                                    //   Communicate cells by scattering from leftover proc
      }                                                         //  Endif for odd number of procs
      stopTimer("Alltoall     ",printNow);                      //  Stop timer & print
      startTimer("Bodies2twigs ");                              //  Start timer
      if( l == LEVEL - 1 ) rbodies2twigs(bodies,twigs);         //  Put recv bodies into twig vector
      stopTimer("Bodies2twigs ",printNow);                      //  Stop timer & print
      startTimer("Cells2twigs  ");                              //  Start timer
      cells2twigs(cells,twigs,l==LEVEL-1);                      //  Put cells into twig vector
      stopTimer("Cells2twigs  ",printNow);                      //  Stop timer & print
      startTimer("Send2twigs   ");                              //  Start timer
      send2twigs(bodies,twigs,offTwigs);                        //  Put send buffer (sticks) into twig vector
      stopTimer("Send2twigs   ",printNow);                      //  Stop timer & print
      startTimer("Recv2twigs   ");                              //  Start timer
      recv2twigs(bodies,twigs);                                 //  Put recv buffer into twig vector
      stopTimer("Recv2twigs   ",printNow);                      //  Stop timer & print
#ifdef DEBUG
      if( l == LEVEL - 1 ) {                                    //  If at last level
        complex SUM = 0;                                        //   Initialize accumulator
        for(C_iter C=twigs.begin(); C!=twigs.end(); ++C) {      //   Loop over twigs
          if( C->NLEAF == 0 ) SUM += C->M[0];                   //    Add multipoles of empty twigs
        }                                                       //   End loop over twigs
        print("Before recv   : ",0);                            //   Print identifier
        print(SUM);                                             //   Print sum of multipoles
      }                                                         //  Endif for last level
#endif
      startTimer("Ziptwigs     ");                              //  Start timer
      zipTwigs(twigs,cells,sticks,l==LEVEL-1);                  //  Zip two groups of twigs that overlap
      stopTimer("Ziptwigs     ",printNow);                      //  Stop timer & print
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
      startTimer("Reindex      ");                              // Start timer
      if( l == LEVEL - 1 ) reindexBodies(bodies,twigs,cells,sticks);// Re-index bodies
      stopTimer("Reindex      ",printNow);                      //  Stop timer & print
      twigs2cells(twigs,cells,sticks);                          //  Turn twigs to cells
      startTimer("Sticks2send  ");                              //  Start timer
      sticks2send(sticks,offTwigs);                             //  Turn sticks to send buffer
      stopTimer("Sticks2send  ",printNow);                      //  Stop timer & print
    }                                                           // End loop over levels of N-D hypercube communication

#ifdef DEBUG
    print("M[0] @ root   : ",0);                                // Print identifier
    print((cells.end()-1)->M[0]);                               // Print monopole of root (should be 1 for test)
    print("bodies.size() : ",0);                                // Print identifier
    print(bodies.size());                                       // Print size of body vector
#endif
    sendCells.clear();                                          // Clear send buffer
    recvCells.clear();                                          // Clear recv buffer
  }
};

#endif
