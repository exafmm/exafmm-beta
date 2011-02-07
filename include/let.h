#ifndef let_h
#define let_h
#include "partition.h"
#include "tree.h"
#include "construct.h"

class LocalEssentialTree : public TreeConstructor,
                           public Partition {
private:
  JBodies sendBodies;                                           // Send buffer for bodies
  JBodies recvBodies;                                           // Receive buffer for bodies
  JCells  sendCells;                                            // Send buffer for cells
  JCells  recvCells;                                            // Receive buffer for cells
public:
  LocalEssentialTree(Bodies &b) : TreeStructure(b),             // Constructor
                                  TreeConstructor(b),
                                  Partition(b) {}
  ~LocalEssentialTree() {}                                      // Destructor

  void commBodies() {                                           // Communicate bodies in LET
    int MPI_TYPE = getType(XMIN[LEVEL-1][0]);                   // Get MPI data type
    std::vector<vect> xmin(SIZE);                               // Buffer for gathering XMIN
    std::vector<vect> xmax(SIZE);                               // Buffer for gathering XMAX
    MPI_Allgather(&XMIN[LEVEL-1][0],3,MPI_TYPE,                 // Gather XMIN
                  &xmin[0][0],3,MPI_TYPE,MPI_COMM_WORLD);
    MPI_Allgather(&XMAX[LEVEL-1][0],3,MPI_TYPE,                 // Gather XMAX
                  &xmax[0][0],3,MPI_TYPE,MPI_COMM_WORLD);
    std::vector<int> srnks,scnts;                               // Ranks to send to, and their send counts
    std::vector<C_iter> scells;                                 // Vector of cell iterators for cells to send
    int oldsize = 0;                                            // Per rank offset of the number of cells to send
    for( int irank=0; irank!=SIZE; ++irank ) {                  // Loop over all ranks
      int ic = 0;                                               //  Initialize neighbor dimension counter
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimensions
        if(xmin[irank][d] < XMAX[LEVEL-1][d] +                  // If the two domains are touching or overlapping
           EPS * std::abs(XMAX[LEVEL-1][d]) &&                  // they are neighbors in this dimension,
           XMIN[LEVEL-1][d] < xmax[irank][d] +                  // and if they are neighbors in all dimensions
           EPS * std::abs(xmax[irank][d]))                      // this is by definition, a neighboring domain
          ic++;                                                 //    Increment neighbor dimension counter
      }                                                         //  End loop over dimensions
      if( ic == 3 && irank != RANK ) {                          //  If ranks are neighbors in all dimensions
        C_iter C(C0);                                           //   Initialize cell iterator
        while( C->NCHILD == 0 ) {                               //   Loop through the twig cells
          vect dist;                                            //    Distance vector
          for( int d=0; d!=3; ++d )                             //    Loop over dimensions
            dist[d] = (C->X[d] > xmax[irank][d])*               //     Calculate the distance between cell C and
                      (C->X[d] - xmax[irank][d])+               //     the nearest point in domain [xmin,xmax]^3
                      (C->X[d] < xmin[irank][d])*
                      (C->X[d] - xmin[irank][d]);
          real R = std::sqrt(norm(dist));                       //    Scalar distance
          if( C0->R + C->R > THETA * R )                        //    If the cell seems close enough for P2P
            scells.push_back(C);                                //     Add cell iterator to scells
          ++C;                                                  //    Increment cell iterator
        }                                                       //   End while loop over twig cells
        srnks.push_back(irank);                                 //   Add current rank to srnks
        scnts.push_back(scells.size()-oldsize);                 //   Add current cell count to scnts
        oldsize = scells.size();                                //   Set new offset for cell count
      }                                                         //  Endif for neighbor ranks
    }                                                           // End loop over all ranks

    int ic = 0;                                                 // Initialize counter for scells
    std::vector<MPI_Request> reqs(2*srnks.size());              // Vector of MPI requests
    std::vector<int>         scnt(srnks.size());                // Vector of send counts
    std::vector<int>         rcnt(srnks.size());                // Vector of receive counts
    for( int i=0; i!=int(srnks.size()); ++i ) {                 // Loop over all ranks to send to & receive from
      int irank = srnks[i];                                     //  Rank to send to & receive from
      for( int c=0; c!=scnts[i]; ++c,++ic ) {                   //  Loop over all cells to send to that rank
        C_iter C = scells[ic];                                  //   Set cell iterator
        for( B_iter B=C->LEAF; B!=C->LEAF+C->NLEAF; ++B ) {     //   Loop over all bodies in that cell
          JBody body;                                           //    Set compact body type for sending
          body.pos  = B->pos;                                   //    Set position of compact body type
          body.scal = B->scal;                                  //    Set mass/charge of compact body type
          sendBodies.push_back(body);                           //    Push it into the send buffer
        }                                                       //   End loop over bodies
      }                                                         //  End loop over cells
      scnt[i] = sendBodies.size();                              //  Set send count of current rank
      MPI_Isend(&scnt[i],1,MPI_INT,irank,0,                     //  Send the send count
                MPI_COMM_WORLD,&reqs[i]);
      MPI_Irecv(&rcnt[i],1,MPI_INT,irank,MPI_ANY_TAG,           //  Receive the recv count
                MPI_COMM_WORLD,&reqs[i+srnks.size()]);
    }                                                           // End loop over ranks
    MPI_Waitall(2*srnks.size(),&reqs[0],MPI_STATUSES_IGNORE);   // Wait for all communication to finish
    int rsize = 0;                                              // Initialize total receive count
    for( int i=0; i!=int(srnks.size()); ++i )                   // Loop over all ranks to receive from
      rsize += rcnt[i];                                         //  Accumulate receive counts

    int bytes = sizeof(sendBodies[0]);                          // Byte size of jbody structure
    recvBodies.resize(rsize);                                   // Receive buffer for bodies
    for( int i=0; i!=int(srnks.size()); ++i ) {                 // Loop over all ranks to send to & receive from
      int irank = srnks[i];                                     // Rank to send to & receive from
      MPI_Isend(&sendBodies[0],scnt[i]*bytes,MPI_BYTE,irank,0,  // Send bodies
                MPI_COMM_WORLD,&reqs[i]);
      MPI_Irecv(&recvBodies[0],rcnt[i]*bytes,MPI_BYTE,irank,MPI_ANY_TAG,// Receive bodies
                MPI_COMM_WORLD,&reqs[i+srnks.size()]);
    }                                                           // End loop over ranks
    MPI_Waitall(2*srnks.size(),&reqs[0],MPI_STATUSES_IGNORE);   // Wait for all communication to finish
    sendBodies.clear();                                         // Clear send buffer for bodies
    recvBodies.clear();                                         // Clear receive buffer for bodies
  }

  void getOtherDomain(vect &xmin, vect &xmax, int l) {
    int MPI_TYPE = getType(XMIN[l][0]);
    vect send[2],recv[2];
    MPI_Request req;
    send[0] = send[1] = XMIN[l];
    recv[0] = recv[1] = 0;
    MPI_Alltoall(send,3,MPI_TYPE,recv,3,MPI_TYPE,MPI_COMM[2]);
    xmin = recv[1-key[2]];
    send[0] = send[1] = XMAX[l];
    recv[0] = recv[1] = 0;
    MPI_Alltoall(send,3,MPI_TYPE,recv,3,MPI_TYPE,MPI_COMM[2]);
    xmax = recv[1-key[2]];
    if( oldnprocs % 2 == 1 && nprocs[0] >= nprocs[1] ) {
      int isend = (key[0] + 1            ) % nprocs[0];
      int irecv = (key[0] - 1 + nprocs[0]) % nprocs[0];
      send[0] = xmin;
      MPI_Isend(send,3,MPI_TYPE,isend,0,MPI_COMM[0],&req);
      MPI_Irecv(recv,3,MPI_TYPE,irecv,0,MPI_COMM[0],&req);
      MPI_Wait(&req,MPI_STATUS_IGNORE);
      if( color[0] != color[1] ) xmin = recv[0];
      send[0] = xmax;
      MPI_Isend(send,3,MPI_TYPE,isend,0,MPI_COMM[0],&req);
      MPI_Irecv(recv,3,MPI_TYPE,irecv,0,MPI_COMM[0],&req);
      MPI_Wait(&req,MPI_STATUS_IGNORE);
      if( color[0] != color[1] ) xmax = recv[0];
    }
    if( oldnprocs == 1 ) {
      xmin = XMIN[l];
      xmax = XMAX[l];
    }
  }

  void getLET(C_iter C, vect xmin, vect xmax) {
    for( int i=0; i!=C->NCHILD; i++ ) {                         // Loop over child cells
      C_iter CC = C0+C->CHILD[i];                               //  Iterator for child cell
      vect dist;                                                //  Distance vector
      for( int d=0; d!=3; ++d )                                 //  Loop over dimensions
        dist[d] = (CC->X[d] > xmax[d])*                         //   Calculate the distance between cell C and
                  (CC->X[d] - xmax[d])+                         //   the nearest point in domain [xmin,xmax]^3
                  (CC->X[d] < xmin[d])*
                  (CC->X[d] - xmin[d]);
      real R = std::sqrt(norm(dist));                           //  Scalar distance
      if( C0->R + CC->R > THETA * R ) {                         //  If the cell seems too close for interaction
        if( CC->NCHILD != 0 )                                   //   If the child cell is not a twig
          getLET(CC,xmin,xmax);                                 //    Traverse the tree further
      } else {                                                  //  If the cell is far enough for interaction
        JCell cell;                                             //   Set compact cell type for sending
        cell.I = CC->I;                                         //   Set index of compact cell type
        cell.M = CC->M;                                         //   Set Multipoles of compact cell type
        sendCells.push_back(cell);                              //   Add cell iterator to send buffer
      }                                                         //  Endif for interaction
    }                                                           // End loop over child cells
  }

  void commCellsAlltoall() {
    int const bytes = sizeof(sendCells[0]);
    int rcnt[2],scnt[2]={0, 0};
    scnt[1-key[2]] = sendCells.size()*bytes;
    MPI_Alltoall(scnt,1,MPI_INT,rcnt,1,MPI_INT,MPI_COMM[2]);
    int sdsp[2] = {0, scnt[0]};
    int rdsp[2] = {0, rcnt[0]};
    if( color[0] != color[1] )
      rcnt[1-key[2]] = 0;

    recvCells.resize(rcnt[1-key[2]]/bytes);
    MPI_Alltoallv(&sendCells[0],scnt,sdsp,MPI_BYTE,
                  &recvCells[0],rcnt,rdsp,MPI_BYTE,MPI_COMM[2]);
  }

  void commCellsScatter() {
    int const bytes = sizeof(sendCells[0]);
    int numScatter = nprocs[1] - 1;
    int oldSize = recvCells.size();
    int *scnt = new int [nprocs[1]];
    int *sdsp = new int [nprocs[1]];
    int rcnt;
    if( key[1] == numScatter ) {
      sdsp[0] = 0;
      for(int i=0; i!=numScatter-1; ++i ) {
        scnt[i] = sendCells.size() / numScatter;
        sdsp[i+1] = sdsp[i] + scnt[i];
      }
      scnt[numScatter-1] = sendCells.size() - (sendCells.size() / numScatter) * (numScatter - 1);
      sdsp[numScatter] = sdsp[numScatter-1] + scnt[numScatter-1];
      scnt[numScatter] = 0;
    }
    MPI_Scatter(scnt,1,MPI_INT,&rcnt,1,MPI_INT,numScatter,MPI_COMM[1]);

    recvCells.resize(oldSize+rcnt);
    for(int i=0; i!= nprocs[1]; ++i ) {
      scnt[i] *= bytes;
      sdsp[i] *= bytes;
    }
    rcnt *= bytes;
    MPI_Scatterv(&sendCells[0],      scnt,sdsp,MPI_BYTE,
                 &recvCells[oldSize],rcnt,     MPI_BYTE,
                 numScatter,MPI_COMM[1]);
    delete[] scnt;
    delete[] sdsp;
  }

  void getTwigs() {
    for( int c=0; c!=int(cells.size()); ++c ) {
      if( cells[c].NCHILD != 0 ) {
        cells.erase(cells.begin()+c);
        c--;
      }
    }
  }

  void unique(int begin, int end) {
    int c_old = begin;
    for( int c=begin; c!=end; ++c ) {
      if( cells[c].I != cells[c_old].I ) {
        c_old = c;
      } else if( c != c_old ) {
        cells[c_old].M += cells[c].M;
        if( cells[c].NLEAF != 0 ) {
          cells[c_old].NLEAF += cells[c].NLEAF;
          cells[c_old].LEAF = cells[c].LEAF;
        }
        cells.erase(cells.begin()+c);
        c--;
        end--;
      }
    }
  }

  void graft() {
    getTwigs();

    BI_iter CI = Icell.begin();
    for( JC_iter JC=recvCells.begin(); JC!=recvCells.end(); ++JC,++CI ) *CI = JC->I;
    JCells jbuffer;
    jbuffer.resize(recvCells.size());
    sort(Icell,recvCells,jbuffer,false);

    for( JC_iter JC=recvCells.begin(); JC!=recvCells.end(); ++JC ) {
      Cell cell;
      cell.I = JC->I;
      cell.M = JC->M;
      cell.NLEAF = cell.NCHILD = 0;
      cell.LEAF = TreeStructure::bodies.end();
      getCenter(cell);
      cells.push_back(cell);
    }

    Cells buffer;
    sortCells(buffer,0,cells.size());
    unique(0,cells.size());

    int begin=0, end=0;
    int level = getLevel(cells[0].I);
    Cells twigs = cells;
    cells.clear();
    for( C_iter C=twigs.begin(); C!=twigs.end(); ++C ) {
      while( getLevel(C->I) != level ) {
        unique(begin,end);
        linkParent(buffer,begin,end);
        level--;
      }
      cells.push_back(*C);
      end++;
    }
    for( int l=level; l>0; --l )                                // Once all the twigs are done, do the rest
      unique(begin,end);
      linkParent(buffer,begin,end);                             //  Form parent-child mutual link

  }

  void commCells() {
    vect xmin=0, xmax=0;
    nprocs[0] = nprocs[1] = SIZE;                               // Initialize number of processes in groups
    offset[0] = offset[1] = 0;                                  // Initialize offset of body in groups
     color[0] =  color[1] =  color[2] = 0;                      // Initialize color of communicators
       key[0] =    key[1] =    key[2] = 0;                      // Initialize key of communicators
    LEVEL = 1;

    for( int l=0; l!=LEVEL; ++l ) {
      multisectionGetComm(l);

      getOtherDomain(xmin,xmax,l+1);
      getLET(cells.end()-1,xmin,xmax);

      commCellsAlltoall();
      if( oldnprocs % 2 == 1 && oldnprocs != 1 && nprocs[0] <= nprocs[1] )
        commCellsScatter();
      graft();
    }
    sendCells.clear();
    recvCells.clear();
  }

};
#endif
