#ifndef parallelfmm_h
#define parallelfmm_h
#include "partition.h"

//! Handles all the communication of local essential trees
class ParallelFMM : public Partition {
private:
  int IRANK;                                                    //!< MPI rank loop counter
  Cells sendCells;                                              //!< Send buffer for cells
  Cells recvCells;                                              //!< Receive buffer for cells
  int * sendCellCount;                                          //!< Send count
  int * sendCellDispl;                                          //!< Send displacement
  int * recvCellCount;                                          //!< Receive count
  int * recvCellDispl;                                          //!< Receive displacement

public:
  using Partition::alltoall;
  using Partition::alltoallv;

private:
//! Get distance to other domain
  real_t getDistance(C_iter C) {
    vec3 dX;                                                    // Distance vector
    for (int d=0; d<3; d++) {                                   // Loop over dimensions
      dX[d] = (C->X[d] + Xperiodic[d] > localXmax[d])*          //  Calculate the distance between cell C and
              (C->X[d] + Xperiodic[d] - localXmax[d])+          //  the nearest point in domain [xmin,xmax]^3
              (C->X[d] + Xperiodic[d] < localXmin[d])*          //  Take the differnece from xmin or xmax
              (C->X[d] + Xperiodic[d] - localXmin[d]);          //  or 0 if between xmin and xmax
    }                                                           // End loop over dimensions
    real_t R2 = norm(dX);                                       // Distance squared
    return R2;                                                  // Return distance squared
  }

//! Add cells to send buffer
  void addSendCell(C_iter C, int &iparent, int &icell) {
    Cell cell(*C);                                              // Initialize send cell
    cell.NCHILD = cell.NCBODY = cell.NDBODY = 0;                // Reset counters
    cell.PARENT = iparent;                                      // Index of parent
    sendCells.push_back(cell);                                  // Push to send cell vector
    icell++;                                                    // Increment cell counter
    C_iter Cparent = sendCells.begin() + sendCellDispl[IRANK] + iparent;// Get parent iterator
    if (Cparent->NCHILD == 0) Cparent->CHILD = icell;           // Index of parent's first child
    Cparent->NCHILD++;                                          // Increment parent's child counter
  }

//! Add bodies to send buffer
  void addSendBody(C_iter C, int &ibody, int icell) {
    C_iter Csend = sendCells.begin() + sendCellDispl[IRANK] + icell;// Get send cell iterator
    Csend->NCBODY = C->NCBODY;                                  // Set number of bodies
    Csend->NDBODY = ibody;                                      // Set body index per rank
    for (B_iter B=C->BODY; B!=C->BODY+C->NCBODY; B++) {         // Loop over bodies in cell
      sendBodies.push_back(*B);                                 //  Push to send body vector
      sendBodies.back().IPROC = IRANK;                          //  Set current rank
    }                                                           // End loop over bodies in cell
    ibody+=C->NCBODY;                                           // Increment body counter
  }

//! Determine which cells to send
  void traverseLET(CellQueue cellQueue) {
    int ibody = 0;                                              // Current send body's offset
    int icell = 0;                                              // Current send cell's offset
    int iparent = 0;                                            // Parent send cell's offset
    int level = int(logf(MPISIZE-1) / M_LN2 / 3) + 1;           // Level of local root cell
    if (MPISIZE == 1) level = 0;                                // Account for serial case
    while (!cellQueue.empty()) {                                // While traversal queue is not empty
      C_iter C = cellQueue.front();                             //  Get front item in traversal queue
      cellQueue.pop();                                          //  Pop item from traversal queue
      for (C_iter CC=Cj0+C->CHILD; CC!=Cj0+C->CHILD+C->NCHILD; CC++) {// Loop over child cells
        addSendCell(CC, iparent, icell);                        //   Add cells to send
        if (CC->NCHILD == 0) {                                  //   If cell is twig
          addSendBody(CC, ibody, icell);                        //    Add bodies to send
        } else {                                                //   If cell is not twig
          bool divide = false;                                  //    Initialize logical for dividing
          if (IMAGES == 0) {                                    //    If free boundary condition
            Xperiodic = 0;                                      //     Set periodic coordinate offset
            real_t R2 = getDistance(CC);                        //     Get distance to other domain
            divide |= 4 * CC->RCRIT * CC->RCRIT > R2;           //     Divide if the cell seems too close
          } else {                                              //    If periodic boundary condition
            for (int ix=-1; ix<=1; ix++) {                      //     Loop over x periodic direction
              for (int iy=-1; iy<=1; iy++) {                    //      Loop over y periodic direction
                for (int iz=-1; iz<=1; iz++) {                  //       Loop over z periodic direction
                  Xperiodic[0] = ix * 2 * globalRadius;         //        Coordinate offset for x periodic direction
                  Xperiodic[1] = iy * 2 * globalRadius;         //        Coordinate offset for y periodic direction
                  Xperiodic[2] = iz * 2 * globalRadius;         //        Coordinate offset for z periodic direction
                  real_t R2 = getDistance(CC);                  //        Get distance to other domain
                  divide |= 4 * CC->RCRIT * CC->RCRIT > R2;     //        Divide if cell seems too close
                }                                               //       End loop over z periodic direction
              }                                                 //      End loop over y periodic direction
            }                                                   //     End loop over x periodic direction
          }                                                     //    Endif for periodic boundary condition
          divide |= CC->R > (globalRadius / (1 << level));      //    Divide if cell is larger than local root cell
          if (!divide) {                                        //    If cell does not have to be divided
            CC->NCHILD = 0;                                     //     Cut off child links
          }                                                     //    Endif for cell division
        }                                                       //   Endif for twig
        cellQueue.push(CC);                                     //   Push cell to traveral queue
      }                                                         //  End loop over child cells
      iparent++;                                                //  Increment parent send cell's offset
    }                                                           // End while loop for traversal queue
  }

//! Exchange send count for cells
  void alltoall(Cells) {
    MPI_Alltoall(sendCellCount, 1, MPI_INT,                     // Communicate send count to get receive count
                 recvCellCount, 1, MPI_INT, MPI_COMM_WORLD);
    recvCellDispl[0] = 0;                                       // Initialize receive displacements
    for (int irank=0; irank<MPISIZE-1; irank++) {               // Loop over ranks
      recvCellDispl[irank+1] = recvCellDispl[irank] + recvCellCount[irank];//  Set receive displacement
    }                                                           // End loop over ranks
  }

//! Exchange cells
  void alltoallv(Cells &cells) {
    int word = sizeof(cells[0]) / 4;                            // Word size of body structure
    recvCells.resize(recvCellDispl[MPISIZE-1]+recvCellCount[MPISIZE-1]);// Resize receive buffer
    for (int irank=0; irank<MPISIZE; irank++) {                 // Loop over ranks
      sendCellCount[irank] *= word;                             //  Multiply send count by word size of data
      sendCellDispl[irank] *= word;                             //  Multiply send displacement by word size of data
      recvCellCount[irank] *= word;                             //  Multiply receive count by word size of data
      recvCellDispl[irank] *= word;                             //  Multiply receive displacement by word size of data
    }                                                           // End loop over ranks
    MPI_Alltoallv(&cells[0], sendCellCount, sendCellDispl, MPI_INT,// Communicate cells
                  &recvCells[0], recvCellCount, recvCellDispl, MPI_INT, MPI_COMM_WORLD);
    for (int irank=0; irank<MPISIZE; irank++) {                 // Loop over ranks
      sendCellCount[irank] /= word;                             //  Divide send count by word size of data
      sendCellDispl[irank] /= word;                             //  Divide send displacement by word size of data
      recvCellCount[irank] /= word;                             //  Divide receive count by word size of data
      recvCellDispl[irank] /= word;                             //  Divide receive displacement by word size of data
    }                                                           // End loop over ranks
  }

public:
//! Constructor
  ParallelFMM() {
    sendCellCount = new int [MPISIZE];                          // Allocate send count
    sendCellDispl = new int [MPISIZE];                          // Allocate send displacement
    recvCellCount = new int [MPISIZE];                          // Allocate receive count
    recvCellDispl = new int [MPISIZE];                          // Allocate receive displacement
  }
//! Destructor
  ~ParallelFMM() {
    delete[] sendCellCount;                                     // Deallocate send count
    delete[] sendCellDispl;                                     // Deallocate send displacement
    delete[] recvCellCount;                                     // Deallocate receive count
    delete[] recvCellDispl;                                     // Deallocate receive displacement
  }

//! Set local essential tree to send to each process
  void setLET(Cells &cells) {
    startTimer("Set LET");                                      // Start timer
    sendBodies.clear();                                         // Clear send buffer for bodies
    sendCells.clear();                                          // Clear send buffer for cells
    sendCellDispl[0] = 0;                                       // Initialize displacement vector
    for (IRANK=0; IRANK<MPISIZE; IRANK++) {                     // Loop over ranks
      if (IRANK != 0) sendCellDispl[IRANK] = sendCellDispl[IRANK-1] + sendCellCount[IRANK-1];// Update displacement
      if (IRANK != MPIRANK) {                                   //  If not current rank
        recvCells = cells;                                      //   Use recvCells as temporary storage
        Cj0 = recvCells.begin();                                //   Set cells begin iterator
        localXmin = rankXmin[IRANK];                            //   Set local Xmin for IRANK
        localXmax = rankXmax[IRANK];                            //   Set local Xmax for IRANK
        Cell cell(*Cj0);                                        //   Send root cell
        cell.NCHILD = cell.NCBODY = cell.NDBODY = 0;            //   Reset link to children and bodies
        sendCells.push_back(cell);                              //   Push it into send buffer
        CellQueue cellQueue;                                    //   Traversal queue
        cellQueue.push(Cj0);                                    //   Push root to traversal queue
        traverseLET(cellQueue);                                 //   Traverse tree to get LET
      }                                                         //  Endif for current rank
      sendCellCount[IRANK] = sendCells.size() - sendCellDispl[IRANK];// Send count for IRANK
    }                                                           // End loop over ranks
    stopTimer("Set LET",printNow);                              // Stop timer
  }

//! Get local essential tree from rank "irank".
  void getLET(Cells &cells, int irank) {
    startTimer("Get LET");                                      // Start timer
    for( int i=recvCellCount[irank]-1; i>=0; i-- ) {            // Loop over receive cells
      C_iter C = recvCells.begin() + recvCellDispl[irank] + i;  //  Iterator of receive cell
      if (C->NCBODY != 0) {                                     //  If cell has bodies
        C->BODY = recvBodies.begin() + recvBodyDispl[irank] + C->NDBODY;// Iterator of first body
        C->NDBODY = C->NCBODY;                                  //   Initialize number of bodies
      }                                                         //  End if for bodies
      if (i != 0) {                                             //  If cell is not root
        C_iter Cparent = recvCells.begin() + recvCellDispl[irank] + C->PARENT;// Iterator of parent cell
        Cparent->NDBODY += C->NDBODY;                           //   Accululate number of bodies
      }                                                         //  End if for root cell
    }                                                           // End loop over receive cells
    cells.resize(recvCellCount[irank]);                         // Resize cell vector for LET
    cells.assign(recvCells.begin()+recvCellDispl[irank],recvCells.begin()+recvCellDispl[irank]+recvCellCount[irank]);
    stopTimer("Get LET",printNow);                              // Stop timer
  }

//! Send bodies
  void commBodies() {
    startTimer("Comm bodies");                                  // Start timer
    alltoall(sendBodies);                                       // Send body count
    alltoallv(sendBodies);                                      // Send bodies
    stopTimer("Comm bodies",printNow);                          // Stop timer
  }

//! Send cells
  void commCells() {
    startTimer("Comm cells");                                   // Start timer
    alltoall(sendCells);                                        // Send cell count
    alltoallv(sendCells);                                       // Senc cells
    stopTimer("Comm cells",printNow);                           // Stop timer
  }
};

#endif
