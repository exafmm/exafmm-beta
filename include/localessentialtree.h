#ifndef localessentialtree_h
#define localessentialtree_h
#include <queue>
#include "partition.h"

//! Handles all the communication of local essential trees
class LocalEssentialTree : public Partition {
 private:
  typedef std::queue<C_iter> CellQueue;                         //!< Queue of cell iterators
  int IRANK;                                                    //!< MPI rank loop counter
  int IMAGES;                                                   //!< Number of periodic image sublevels
  fvec3 localXmin;                                              //!< Local Xmin for a given rank
  fvec3 localXmax;                                              //!< Local Xmax for a given rank
  fvec3 * allLocalXmin;                                         //!< Array for local Xmin for all ranks
  fvec3 * allLocalXmax;                                         //!< Array for local Xmax for all ranks
  Cells sendCells;                                              //!< Send buffer for cells
  Cells recvCells;                                              //!< Receive buffer for cells
  C_iter C0;                                                    //!< Iterator of first cell
  int * sendCellCount;                                          //!< Send count
  int * sendCellDispl;                                          //!< Send displacement
  int * recvCellCount;                                          //!< Receive count
  int * recvCellDispl;                                          //!< Receive displacement

 public:
  using Partition::alltoall;
  using Partition::alltoallv;

 private:
  void allgatherBounds(Bounds local) {
    fvec3 Xmin, Xmax;
    for (int d=0; d<3; d++) {                                   // Loop over dimensions
      Xmin[d] = local.Xmin[d];                                  //  Convert Xmin to float
      Xmax[d] = local.Xmax[d];                                  //  Convert Xmax to float
    }                                                           // End loop over dimensions
    MPI_Allgather(Xmin, 3, MPI_FLOAT, &allLocalXmin[0], 3, MPI_FLOAT, MPI_COMM_WORLD);// Gather all domain bounds
    MPI_Allgather(Xmax, 3, MPI_FLOAT, &allLocalXmax[0], 3, MPI_FLOAT, MPI_COMM_WORLD);// Gather all domain bounds
  }

//! Get distance to other domain
  real_t getDistance(C_iter C, vec3 Xperiodic) {
    vec3 dX;                                                    // Distance vector
    for (int d=0; d<3; d++) {                                   // Loop over dimensions
      dX[d] = (C->X[d] + Xperiodic[d] > localXmax[d])*          //  Calculate the distance between cell C and
              (C->X[d] + Xperiodic[d] - localXmax[d])+          //  the nearest point in domain [xmin,xmax]^3
              (C->X[d] + Xperiodic[d] < localXmin[d])*          //  Take the differnece from xmin or xmax
              (C->X[d] + Xperiodic[d] - localXmin[d]);          //  or 0 if between xmin and xmax
    }                                                           // End loop over dimensions
    return norm(dX);                                            // Return distance squared
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
  void traverseLET(CellQueue cellQueue, real_t cycle) {
    int ibody = 0;                                              // Current send body's offset
    int icell = 0;                                              // Current send cell's offset
    int iparent = 0;                                            // Parent send cell's offset
    int level = int(logf(MPISIZE-1) / M_LN2 / 3) + 1;           // Level of local root cell
    if (MPISIZE == 1) level = 0;                                // Account for serial case
    while (!cellQueue.empty()) {                                // While traversal queue is not empty
      C_iter C = cellQueue.front();                             //  Get front item in traversal queue
      cellQueue.pop();                                          //  Pop item from traversal queue
      for (C_iter CC=C0+C->CHILD; CC!=C0+C->CHILD+C->NCHILD; CC++) {// Loop over child cells
        addSendCell(CC, iparent, icell);                        //   Add cells to send
        if (CC->NCHILD == 0) {                                  //   If cell is twig
          addSendBody(CC, ibody, icell);                        //    Add bodies to send
        } else {                                                //   If cell is not twig
          bool divide = false;                                  //    Initialize logical for dividing
          vec3 Xperiodic = 0;                                   //    Periodic coordinate offset
          if (IMAGES == 0) {                                    //    If free boundary condition
            real_t R2 = getDistance(CC, Xperiodic);             //     Get distance to other domain
            divide |= 4 * CC->RCRIT * CC->RCRIT > R2;           //     Divide if the cell seems too close
          } else {                                              //    If periodic boundary condition
            for (int ix=-1; ix<=1; ix++) {                      //     Loop over x periodic direction
              for (int iy=-1; iy<=1; iy++) {                    //      Loop over y periodic direction
                for (int iz=-1; iz<=1; iz++) {                  //       Loop over z periodic direction
                  Xperiodic[0] = ix * cycle;                    //        Coordinate offset for x periodic direction
                  Xperiodic[1] = iy * cycle;                    //        Coordinate offset for y periodic direction
                  Xperiodic[2] = iz * cycle;                    //        Coordinate offset for z periodic direction
                  real_t R2 = getDistance(CC, Xperiodic);       //        Get distance to other domain
                  divide |= 4 * CC->RCRIT * CC->RCRIT > R2;     //        Divide if cell seems too close
                }                                               //       End loop over z periodic direction
              }                                                 //      End loop over y periodic direction
            }                                                   //     End loop over x periodic direction
          }                                                     //    Endif for periodic boundary condition
          divide |= CC->R > (cycle / (1 << (level+1)));         //    Divide if cell is larger than local root cell
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
  LocalEssentialTree(int images) : IMAGES(images) {
    allLocalXmin = new fvec3 [MPISIZE];                         // Allocate array for minimum of local domains
    allLocalXmax = new fvec3 [MPISIZE];                         // Allocate array for maximum of local domains
    sendCellCount = new int [MPISIZE];                          // Allocate send count
    sendCellDispl = new int [MPISIZE];                          // Allocate send displacement
    recvCellCount = new int [MPISIZE];                          // Allocate receive count
    recvCellDispl = new int [MPISIZE];                          // Allocate receive displacement
  }
//! Destructor
  ~LocalEssentialTree() {
    delete[] allLocalXmin;                                      // Deallocate array for minimum of local domains
    delete[] allLocalXmax;                                      // Deallocate array for maximum of local domains
    delete[] sendCellCount;                                     // Deallocate send count
    delete[] sendCellDispl;                                     // Deallocate send displacement
    delete[] recvCellCount;                                     // Deallocate receive count
    delete[] recvCellDispl;                                     // Deallocate receive displacement
  }

//! Set local essential tree to send to each process
  void setLET(Cells &cells, Bounds bounds, real_t cycle) {
    startTimer("Set LET");                                      // Start timer
    allgatherBounds(bounds);                                    // Gather local bounds from all ranks
    sendBodies.clear();                                         // Clear send buffer for bodies
    sendCells.clear();                                          // Clear send buffer for cells
    sendCellDispl[0] = 0;                                       // Initialize displacement vector
    for (IRANK=0; IRANK<MPISIZE; IRANK++) {                     // Loop over ranks
      if (IRANK != 0) sendCellDispl[IRANK] = sendCellDispl[IRANK-1] + sendCellCount[IRANK-1];// Update displacement
      if (IRANK != MPIRANK) {                                   //  If not current rank
        recvCells = cells;                                      //   Use recvCells as temporary storage
        C0 = recvCells.begin();                                 //   Set cells begin iterator
        localXmin = allLocalXmin[IRANK];                        //   Set local Xmin for IRANK
        localXmax = allLocalXmax[IRANK];                        //   Set local Xmax for IRANK
        Cell cell(*C0);                                         //   Send root cell
        cell.NCHILD = cell.NCBODY = cell.NDBODY = 0;            //   Reset link to children and bodies
        sendCells.push_back(cell);                              //   Push it into send buffer
        CellQueue cellQueue;                                    //   Traversal queue
        cellQueue.push(C0);                                     //   Push root to traversal queue
        traverseLET(cellQueue,cycle);                           //   Traverse tree to get LET
      }                                                         //  Endif for current rank
      sendCellCount[IRANK] = sendCells.size() - sendCellDispl[IRANK];// Send count for IRANK
    }                                                           // End loop over ranks
    stopTimer("Set LET",printNow);                              // Stop timer
  }

//! Get local essential tree from rank "irank".
  void getLET(Cells &cells, int irank) {
    startTimer("Get LET");                                      // Start timer
    for (int i=recvCellCount[irank]-1; i>=0; i--) {             // Loop over receive cells
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
  Bodies commBodies() {
    startTimer("Comm bodies");                                  // Start timer
    alltoall(sendBodies);                                       // Send body count
    alltoallv(sendBodies);                                      // Send bodies
    stopTimer("Comm bodies",printNow);                          // Stop timer
    return recvBodies;                                          // Return received bodies
  }

//! Send bodies
  Bodies commBodies(Bodies bodies) {
    startTimer("Comm bodies");                                  // Start timer
    alltoall(bodies);                                           // Send body count
    alltoallv(bodies);                                          // Send bodies
    stopTimer("Comm bodies",printNow);                          // Stop timer
    return recvBodies;                                          // Return received bodies
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
