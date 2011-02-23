#ifndef let_h
#define let_h
#include "partition.h"

class LocalEssentialTree : public Partition {
private:
  JBodies sendBodies;                                           // Send buffer for bodies
  JBodies recvBodies;                                           // Receive buffer for bodies
  JCells  sendCells;                                            // Send buffer for cells
  JCells  recvCells;                                            // Receive buffer for cells

private:
  void getOtherDomain(vect &xmin, vect &xmax, int l) {          // Get boundries of domains on other processes
    int MPI_TYPE = getType(XMIN[l][0]);                         // Get MPI data type
    vect send[2],recv[2];                                       // Send and recv buffer
    MPI_Request req;                                            // MPI requests
    send[0] = send[1] = XMIN[l];                                // Set XMIN into send buffer
    recv[0] = recv[1] = 0;                                      // Initialize recv buffer
    MPI_Alltoall(send,3,MPI_TYPE,recv,3,MPI_TYPE,MPI_COMM[2]);  // Communicate XMIN
    xmin = recv[1-key[2]];                                      // Copy xmin from recv buffer
    send[0] = send[1] = XMAX[l];                                // Set XMAX into send buffer
    recv[0] = recv[1] = 0;                                      // Initialize recv buffer
    MPI_Alltoall(send,3,MPI_TYPE,recv,3,MPI_TYPE,MPI_COMM[2]);  // Communicate XMAX
    xmax = recv[1-key[2]];                                      // Copy xmax from recv buffer
    if( oldnprocs % 2 == 1 && nprocs[0] >= nprocs[1] ) {        // If right half of odd group
      int isend = (key[0] + 1            ) % nprocs[0];         //  Send to next rank (wrapped)
      int irecv = (key[0] - 1 + nprocs[0]) % nprocs[0];         //  Recv from previous rank (wrapped)
      send[0] = xmin;                                           //  Set xmin in send buffer
      MPI_Isend(send,3,MPI_TYPE,isend,0,MPI_COMM[0],&req);      //  Send to next rank
      MPI_Irecv(recv,3,MPI_TYPE,irecv,0,MPI_COMM[0],&req);      //  Recv from previous rank
      MPI_Wait(&req,MPI_STATUS_IGNORE);                         //  Wait for recv to finish
      if( color[0] != color[1] ) xmin = recv[0];                //  Update only if leftover process of odd group
      send[0] = xmax;                                           //  Set xmax in send buffer
      MPI_Isend(send,3,MPI_TYPE,isend,0,MPI_COMM[0],&req);      //  Send to next rank
      MPI_Irecv(recv,3,MPI_TYPE,irecv,0,MPI_COMM[0],&req);      //  Recv from previous rank
      MPI_Wait(&req,MPI_STATUS_IGNORE);                         //  Wait for recv to finish
      if( color[0] != color[1] ) xmax = recv[0];                //  Update only if leftover process of odd group
    }                                                           // Endif for right half of odd group
    if( oldnprocs == 1 ) {                                      // If previous group had one process
      xmin = XMIN[l];                                           //  Use own XMIN value for xmin
      xmax = XMAX[l];                                           //  Use own XMAX value for xmax
    }                                                           // Endif for isolated process
  }

  real getDistance(C_iter C, vect xmin, vect xmax) {            // Get disatnce to other domain
    vect dist;                                                  // Distance vector
    for( int d=0; d!=3; ++d )                                   // Loop over dimensions
      dist[d] = (C->X[d] > xmax[d])*                            //  Calculate the distance between cell C and
                (C->X[d] - xmax[d])+                            //  the nearest point in domain [xmin,xmax]^3
                (C->X[d] < xmin[d])*                            //  Take the differnece from xmin or xmax
                (C->X[d] - xmin[d]);                            //  or 0 if between xmin and xmax
    real R = std::sqrt(norm(dist));                             // Scalar distance
    return R;
  }

  void getLET(C_iter C0, C_iter C, vect xmin, vect xmax) {      // Determine which cells to send
    for( int i=0; i!=C->NCHILD; i++ ) {                         // Loop over child cells
      C_iter CC = C0+C->CHILD[i];                               //  Iterator for child cell
      real R = getDistance(CC,xmin,xmax);                       //  Get distance to other domain
      if( 3 * CC->R > THETA * R && CC->NCHILD != 0 ) {          //  If the cell seems too close and not twig
        getLET(C0,CC,xmin,xmax);                                //   Traverse the tree further
      } else {                                                  //  If the cell if far or a twig
        JCell cell;                                             //   Set compact cell type for sending
        cell.I = CC->I;                                         //   Set index of compact cell type
        cell.M = CC->M;                                         //   Set Multipoles of compact cell type
        sendCells.push_back(cell);                              //   Push cell into send buffer vector
      }                                                         //  Endif for interaction
    }                                                           // End loop over child cells
    if( C->I == 0 && C->NCHILD == 0 ) {                         // If the root cell has no children
      JCell cell;                                               //  Set compact cell type for sending
      cell.I = C->I;                                            //  Set index of compact cell type
      cell.M = C->M;                                            //  Set Multipoles of compact cell type
      sendCells.push_back(cell);                                //  Push cell into send buffer vector
    }                                                           // Endif for root cells children
  }

  void commCellsAlltoall() {                                    // Communicate cells by one-to-one MPI_Alltoallv
    int const bytes = sizeof(sendCells[0]);                     // Byte size of JCell structure
    int rcnt[2], scnt[2] = {0, 0};                              // Recv count, send count
    scnt[1-key[2]] = sendCells.size()*bytes;                    // Set send count to size of send buffer * bytes
    MPI_Alltoall(scnt,1,MPI_INT,rcnt,1,MPI_INT,MPI_COMM[2]);    // Communicate the send count to get recv count
    int sdsp[2] = {0, scnt[0]};                                 // Send displacement
    int rdsp[2] = {0, rcnt[0]};                                 // Recv displacement
    if( color[0] != color[1] ) {                                // If leftover process of odd group
      rcnt[1-key[2]] = 0;                                       //  It won't have a pair for this communication
    }                                                           // Endif for leftover process of odd group
    recvCells.resize(rcnt[1-key[2]]/bytes);                     // Resize recv buffer
    MPI_Alltoallv(&sendCells[0],scnt,sdsp,MPI_BYTE,             // Communicate cells
                  &recvCells[0],rcnt,rdsp,MPI_BYTE,MPI_COMM[2]);// MPI_COMM[2] is for the one-to-one pair
  }

  void commCellsScatter() {                                     // Communicate cells by scattering from leftover proc
    int const bytes = sizeof(sendCells[0]);                     // Byte size of JCell structure
    int numScatter = nprocs[1] - 1;                             // Number of processes to scatter to
    int oldSize = recvCells.size();                             // Size of recv buffer before communication
    int *scnt = new int [nprocs[1]];                            // Send count
    int *sdsp = new int [nprocs[1]];                            // Send displacement
    int rcnt;                                                   // Recv count
    if( key[1] == numScatter ) {                                // If this is the leftover proc to scatter from
      sdsp[0] = 0;                                              //  Initialize send displacement
      for(int i=0; i!=numScatter; ++i ) {                       //  Loop over processes to send to
        int begin = 0, end = sendCells.size();                  //   Set begin and end of range to send
        splitRange(begin,end,i,numScatter);                     //   Split range into numScatter and get i-th range
        scnt[i] = end - begin;                                  //   Set send count based on range
        sdsp[i+1] = sdsp[i] + scnt[i];                          //   Set send displacement based on send count
      }                                                         //  End loop over processes to send to
      scnt[numScatter] = 0;                                     //  Send count to self should be 0
    }                                                           // Endif for leftover proc
    MPI_Scatter(scnt,1,MPI_INT,&rcnt,1,MPI_INT,numScatter,MPI_COMM[1]);// Scatter the send count to get recv count

    recvCells.resize(oldSize+rcnt);                             // Resize recv buffer based on recv count
    for(int i=0; i!= nprocs[1]; ++i ) {                         // Loop over group of processes
      scnt[i] *= bytes;                                         //  Multiply send count by byte size of data
      sdsp[i] *= bytes;                                         //  Multiply send displacement by byte size of data
    }                                                           // End loop over group of processes
    rcnt *= bytes;                                              // Multiply recv count by byte size of data
    MPI_Scatterv(&sendCells[0],      scnt,sdsp,MPI_BYTE,        // Communicate cells via MPI_Scatterv
                 &recvCells[oldSize],rcnt,     MPI_BYTE,        // Offset recv buffer by oldSize
                 numScatter,MPI_COMM[1]);
    delete[] scnt;                                              // Delete send count
    delete[] sdsp;                                              // Delete send displacement
  }

  void rbodies2twigs(Bodies &bodies, Cells &twigs) {            // Turn recv bodies to twigs
    for( JB_iter JB=recvBodies.begin(); JB!=recvBodies.end(); ++JB ) {// Loop over recv bodies
      Body body;                                                //  Body structure
      body.I    = JB->I;                                        //  Set index of body
      body.pos  = JB->pos;                                      //  Set position of body
      body.scal = JB->scal;                                     //  Set mass/charge of body
      bodies.push_back(body);                                   //  Push body into bodies vector
    }                                                           // End loop over recv bodies
    buffer.resize(bodies.size());                               // Resize sort buffer
    sort(bodies,buffer,false);                                  // Sort bodies in descending order
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
      cell.I = JC->I;                                           //  Set index of cell
      cell.M = JC->M;                                           //  Set multipole of cell
      cell.NLEAF = cell.NCHILD = 0;                             //  Set number of leafs and children
      cell.LEAF = bodies.end();                                 //  Set pointer to first leaf
      getCenter(cell);                                          //  Set center and radius
      twigs.push_back(cell);                                    //  Push cell into twig vector
    }                                                           // End loop over send buffer
    sendCells.clear();                                          // Clear send buffer
  }

  void recv2twigs(Bodies &bodies, Cells &twigs) {               // Turn recv buffer to twigs
    for( JC_iter JC=recvCells.begin(); JC!=recvCells.end(); ++JC ) {// Loop over recv buffer
      Cell cell;                                                //  Cell structure
      cell.I = JC->I;                                           //  Set index of cell
      cell.M = JC->M;                                           //  Set multipole of cell
      cell.NLEAF = cell.NCHILD = 0;                             //  Set number of leafs and children
      cell.LEAF = bodies.end();                                 //  Set pointer to first leaf
      getCenter(cell);                                          //  Set center and radius
      twigs.push_back(cell);                                    //  Push cell into twig vector
    }                                                           // End loop over recv buffer
  }

  void zipTwigs(Cells &twigs, Cells &cells, Cells &sticks, bool last) {// Zip two groups of twigs that overlap
    sortCells(twigs);                                           // Sort twigs in ascending order
    bigint index = -1;                                          // Initialize index counter
    while( !twigs.empty() ) {                                   // While twig vector is not empty
      if( twigs.back().I != index ) {                           //  If twig's index is different from previous
        cells.push_back(twigs.back());                          //   Push twig into cell vector
        index = twigs.back().I;                                 //   Update index counter
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
    BottomUp::setIndex(bodies,0,0,0,true);                      // Set index of bodies
    buffer.resize(bodies.size());                               // Resize sort buffer
    sort(bodies,buffer);                                        // Sort bodies in ascending order
    BottomUp::grow(bodies);                                     // Grow tree structure
    sort(bodies,buffer);                                        // Sort bodies in ascending order
    bodies2twigs(bodies,twigs);                                 // Turn bodies to twigs
    for( C_iter C=twigs.begin(); C!=twigs.end(); ++C ) {        // Loop over cells
      if( sticks.size() > 0 ) {                                 //  If stick vector is not empty
        if( C->I == sticks.back().I ) {                         //   If twig's index is equal to stick's index
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
      cell.I = sticks.back().I;                                 //  Set index of cell
      cell.M = sticks.back().M;                                 //  Set multipole of cell
      sendCells.push_back(cell);                                //  Push cell into send buffer
      sticks.pop_back();                                        //  Pop last element of stick vector
    }                                                           // End while for stick vector
    offTwigs = sendCells.size();                                // Keep track of current send buffer size
  }

public:
  LocalEssentialTree() : Partition() {}                         // Constructor
  ~LocalEssentialTree() {}                                      // Destructor

  void commBodies(Cells &cells) {                               // Communicate bodies in the LET
    int MPI_TYPE = getType(XMIN[LEVEL][0]);                     // Get MPI data type
    std::vector<vect> xmin(SIZE);                               // Buffer for gathering XMIN
    std::vector<vect> xmax(SIZE);                               // Buffer for gathering XMAX
    MPI_Allgather(&XMIN[LEVEL][0],3,MPI_TYPE,                   // Gather XMIN
                  &xmin[0][0],3,MPI_TYPE,MPI_COMM_WORLD);
    MPI_Allgather(&XMAX[LEVEL][0],3,MPI_TYPE,                   // Gather XMAX
                  &xmax[0][0],3,MPI_TYPE,MPI_COMM_WORLD);
    std::vector<int> srnks,scnts;                               // Ranks to send to, and their send counts
    std::vector<C_iter> scells;                                 // Vector of cell iterators for cells to send
    int oldsize = 0;                                            // Per rank offset of the number of cells to send
    C_iter C0 = cells.begin();                                  // Set cell begin iterator
    for( int irank=0; irank!=SIZE; ++irank ) {                  // Loop over ranks
      int ic = 0;                                               //  Initialize neighbor dimension counter
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimensions
        if(xmin[irank][d] < XMAX[LEVEL][d] + C0->R &&           //   If the two domains are touching or overlapping
           XMIN[LEVEL][d] < xmax[irank][d] + C0->R) {           //   in all dimensions, they are neighboring domains
          ic++;                                                 //    Increment neighbor dimension counter
        }                                                       //   Endif for overlapping domains
      }                                                         //  End loop over dimensions
      if( ic == 3 && irank != RANK ) {                          //  If ranks are neighbors in all dimensions
        for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {    //   Loop over cells
          if( C->NCHILD == 0 ) {                                //    If cell is a twig
            real R = getDistance(C,xmin[irank],xmax[irank]);    //     Get distance to other domain
            if( 3 * C->R > THETA * R ) {                        //     If the cell seems close enough for P2P
              scells.push_back(C);                              //      Add cell iterator to scells
            }                                                   //     Endif for cell distance
          }                                                     //    Endif for twigs
        }                                                       //   End loop over cells
        srnks.push_back(irank);                                 //   Add current rank to srnks
        scnts.push_back(scells.size()-oldsize);                 //   Add current cell count to scnts
        oldsize = scells.size();                                //   Set new offset for cell count
      }                                                         //  Endif for neighbor ranks
    }                                                           // End loop over ranks

    int ic = 0, ssize = 0;                                      // Initialize counter and offset for scells
    std::vector<MPI_Request> reqs(2*srnks.size());              // Vector of MPI requests
    std::vector<int>         scnt(srnks.size());                // Vector of send counts
    std::vector<int>         rcnt(srnks.size());                // Vector of recv counts
    for( int i=0; i!=int(srnks.size()); ++i ) {                 // Loop over ranks to send to & recv from
      int irank = srnks[i];                                     //  Rank to send to & recv from
      for( int c=0; c!=scnts[i]; ++c,++ic ) {                   //  Loop over cells to send to that rank
        C_iter C = scells[ic];                                  //   Set cell iterator
        for( B_iter B=C->LEAF; B!=C->LEAF+C->NLEAF; ++B ) {     //   Loop over bodies in that cell
          JBody body;                                           //    Set compact body type for sending
          body.I    = B->I;                                     //    Set cell index of compact body type
          body.pos  = B->pos;                                   //    Set position of compact body type
          body.scal = B->scal;                                  //    Set mass/charge of compact body type
          sendBodies.push_back(body);                           //    Push it into the send buffer
        }                                                       //   End loop over bodies
      }                                                         //  End loop over cells
      scnt[i] = sendBodies.size()-ssize;                        //  Set send count of current rank
      ssize += scnt[i];                                         //  Increment offset for vector scells
      MPI_Isend(&scnt[i],1,MPI_INT,irank,0,                     //  Send the send count
                MPI_COMM_WORLD,&reqs[i]);
      MPI_Irecv(&rcnt[i],1,MPI_INT,irank,MPI_ANY_TAG,           //  Receive the recv count
                MPI_COMM_WORLD,&reqs[i+srnks.size()]);
    }                                                           // End loop over ranks
    MPI_Waitall(2*srnks.size(),&reqs[0],MPI_STATUSES_IGNORE);   // Wait for all communication to finish

    int rsize = 0;                                              // Initialize total recv count
    for( int i=0; i!=int(srnks.size()); ++i )                   // Loop over ranks to recv from
      rsize += rcnt[i];                                         //  Accumulate recv counts
    recvBodies.resize(rsize);                                   // Receive buffer for bodies
    int bytes = sizeof(sendBodies[0]);                          // Byte size of jbody structure
    rsize = ssize = 0;                                          // Initialize offset for sendBodies & recvBodies
    for( int i=0; i!=int(srnks.size()); ++i ) {                 // Loop over ranks to send to & recv from
      int irank = srnks[i];                                     // Rank to send to & recv from
      MPI_Isend(&sendBodies[ssize],scnt[i]*bytes,MPI_BYTE,irank,0,// Send bodies
                MPI_COMM_WORLD,&reqs[i]);
      MPI_Irecv(&recvBodies[rsize],rcnt[i]*bytes,MPI_BYTE,irank,MPI_ANY_TAG,// Receive bodies
                MPI_COMM_WORLD,&reqs[i+srnks.size()]);
      ssize += scnt[i];                                         // Increment offset for sendBodies
      rsize += rcnt[i];                                         // Increment offset for recvBodies
    }                                                           // End loop over ranks
    MPI_Waitall(2*srnks.size(),&reqs[0],MPI_STATUSES_IGNORE);   // Wait for all communication to finish
    sendBodies.clear();                                         // Clear send buffer for bodies
  }

  void commCells(Bodies &bodies, Cells &cells) {                // Communicate cell in the LET
    int offTwigs = 0;                                           // Initialize offset of twigs
    vect xmin = 0, xmax = 0;                                    // Initialize domain boundaries
    Cells twigs,sticks;                                         // Twigs and sticks are special types of cells
    nprocs[0] = nprocs[1] = SIZE;                               // Initialize number of processes in groups
    offset[0] = offset[1] = 0;                                  // Initialize offset of body in groups
     color[0] =  color[1] =  color[2] = 0;                      // Initialize color of communicators
       key[0] =    key[1] =    key[2] = 0;                      // Initialize key of communicators

    for( int l=0; l!=LEVEL; ++l ) {                             // Loop over levels of N-D hypercube communication
      bisectionGetComm(l);                                      //  Split the MPI communicator for that level
      getOtherDomain(xmin,xmax,l+1);                            //  Get boundries of domains on other processes
      getLET(cells.begin(),cells.end()-1,xmin,xmax);            //  Determine which cells to send
      commCellsAlltoall();                                      //  Communicate cells by one-to-one MPI_Alltoallv
      if( oldnprocs % 2 == 1 && oldnprocs != 1 && nprocs[0] <= nprocs[1] ) {// If scatter is necessary
        commCellsScatter();                                     //   Communicate cells by scattering from leftover proc
      }                                                         //  Endif for odd number of procs
      if( l == LEVEL - 1 ) rbodies2twigs(bodies,twigs);         //  Put recv bodies into twig vector
      cells2twigs(cells,twigs,l==LEVEL-1);                      //  Put cells into twig vector
      send2twigs(bodies,twigs,offTwigs);                        //  Put send buffer (sticks) into twig vector
      recv2twigs(bodies,twigs);                                 //  Put recv buffer into twig vector
#ifdef DEBUG
      if( l == LEVEL - 1 ) {                                    //  If at last level
        real SUM = 0;                                           //   Initialize accumulator
        for(C_iter C=twigs.begin(); C!=twigs.end(); ++C) {      //   Loop over twigs
          if( C->NLEAF == 0 ) SUM += C->M[0];                   //    Add multipoles of empty twigs
        }                                                       //   End loop over twigs
        print("Before recv   : ",0);                            //   Print identifier
        print(SUM);                                             //   Print sum of multipoles
      }                                                         //  Endif for last level
#endif
      zipTwigs(twigs,cells,sticks,l==LEVEL-1);                  //  Zip two groups of twigs that overlap
#ifdef DEBUG
      if( l == LEVEL - 1 ) {                                    //  If at last level
        real SUM = 0;                                           //   Initialize accumulator
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
      sticks2send(sticks,offTwigs);                             //  Turn sticks to send buffer
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
