#ifndef tree_mpi_h
#define tree_mpi_h
#include "kernel.h"
#include "logger.h"

//! Handles all the communication of local essential trees
class TreeMPI {
protected:
  const int mpirank;                                            //!< Rank of MPI communicator
  const int mpisize;                                            //!< Size of MPI communicator
  const int images;                                             //!< Number of periodic image sublevels
  float (* allBoundsXmin)[3];                                   //!< Array for local Xmin for all ranks
  float (* allBoundsXmax)[3];                                   //!< Array for local Xmax for all ranks
  Bodies sendBodies;                                            //!< Send buffer for bodies
  Bodies recvBodies;                                            //!< Receive buffer for bodies
  Cells sendCells;                                              //!< Send buffer for cells
  Cells recvCells;                                              //!< Receive buffer for cells
  int * sendBodyCount;                                          //!< Send count
  int * sendBodyDispl;                                          //!< Send displacement
  int * recvBodyCount;                                          //!< Receive count
  int * recvBodyDispl;                                          //!< Receive displacement
  int * sendCellCount;                                          //!< Send count
  int * sendCellDispl;                                          //!< Send displacement
  int * recvCellCount;                                          //!< Receive count
  int * recvCellDispl;                                          //!< Receive displacement

private:
  //! Exchange send count for bodies
  void alltoall(Bodies & bodies) {
    for (int i=0; i<mpisize; i++) {                             // Loop over ranks
      sendBodyCount[i] = 0;                                     //  Initialize send counts
    }                                                           // End loop over ranks
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      assert(0 <= B->IRANK && B->IRANK < mpisize);              //  Check bounds for process ID
      sendBodyCount[B->IRANK]++;                                //  Fill send count bucket
      B->IRANK = mpirank;                                       //  Tag for sending back to original rank
    }                                                           // End loop over bodies
    MPI_Alltoall(sendBodyCount, 1, MPI_INT,                     // Communicate send count to get receive count
                 recvBodyCount, 1, MPI_INT, MPI_COMM_WORLD);
    sendBodyDispl[0] = recvBodyDispl[0] = 0;                    // Initialize send/receive displacements
    for (int irank=0; irank<mpisize-1; irank++) {               // Loop over ranks
      sendBodyDispl[irank+1] = sendBodyDispl[irank] + sendBodyCount[irank];//  Set send displacement
      recvBodyDispl[irank+1] = recvBodyDispl[irank] + recvBodyCount[irank];//  Set receive displacement
    }                                                           // End loop over ranks
  }

  //! Exchange bodies
  void alltoallv(Bodies & bodies) {
    assert( (sizeof(bodies[0]) & 3) == 0 );                     // Body structure must be 4 Byte aligned
    int word = sizeof(bodies[0]) / 4;                           // Word size of body structure
    recvBodies.resize(recvBodyDispl[mpisize-1]+recvBodyCount[mpisize-1]);// Resize receive buffer
    for (int irank=0; irank<mpisize; irank++) {                 // Loop over ranks
      sendBodyCount[irank] *= word;                             //  Multiply send count by word size of data
      sendBodyDispl[irank] *= word;                             //  Multiply send displacement by word size of data
      recvBodyCount[irank] *= word;                             //  Multiply receive count by word size of data
      recvBodyDispl[irank] *= word;                             //  Multiply receive displacement by word size of data
    }                                                           // End loop over ranks
    MPI_Alltoallv((int*)&bodies[0], sendBodyCount, sendBodyDispl, MPI_INT,// Communicate bodies
                  (int*)&recvBodies[0], recvBodyCount, recvBodyDispl, MPI_INT, MPI_COMM_WORLD);
    for (int irank=0; irank<mpisize; irank++) {                 // Loop over ranks
      sendBodyCount[irank] /= word;                             //  Divide send count by word size of data
      sendBodyDispl[irank] /= word;                             //  Divide send displacement by word size of data
      recvBodyCount[irank] /= word;                             //  Divide receive count by word size of data
      recvBodyDispl[irank] /= word;                             //  Divide receive displacement by word size of data
    }                                                           // End loop over ranks
  }

  //! Exchange send count for cells
  void alltoall(Cells) {
    MPI_Alltoall(sendCellCount, 1, MPI_INT,                     // Communicate send count to get receive count
                 recvCellCount, 1, MPI_INT, MPI_COMM_WORLD);
    recvCellDispl[0] = 0;                                       // Initialize receive displacements
    for (int irank=0; irank<mpisize-1; irank++) {               // Loop over ranks
      recvCellDispl[irank+1] = recvCellDispl[irank] + recvCellCount[irank];//  Set receive displacement
    }                                                           // End loop over ranks
  }

  //! Exchange cells
  void alltoallv(Cells & cells) {
    assert( (sizeof(cells[0]) & 3) == 0 );                      // Cell structure must be 4 Byte aligned
    int word = sizeof(cells[0]) / 4;                            // Word size of body structure
    recvCells.resize(recvCellDispl[mpisize-1]+recvCellCount[mpisize-1]);// Resize receive buffer
    for (int irank=0; irank<mpisize; irank++) {                 // Loop over ranks
      sendCellCount[irank] *= word;                             //  Multiply send count by word size of data
      sendCellDispl[irank] *= word;                             //  Multiply send displacement by word size of data
      recvCellCount[irank] *= word;                             //  Multiply receive count by word size of data
      recvCellDispl[irank] *= word;                             //  Multiply receive displacement by word size of data
    }                                                           // End loop over ranks
    MPI_Alltoallv((int*)&cells[0], sendCellCount, sendCellDispl, MPI_INT,// Communicate cells
                  (int*)&recvCells[0], recvCellCount, recvCellDispl, MPI_INT, MPI_COMM_WORLD);
    for (int irank=0; irank<mpisize; irank++) {                 // Loop over ranks
      sendCellCount[irank] /= word;                             //  Divide send count by word size of data
      sendCellDispl[irank] /= word;                             //  Divide send displacement by word size of data
      recvCellCount[irank] /= word;                             //  Divide receive count by word size of data
      recvCellDispl[irank] /= word;                             //  Divide receive displacement by word size of data
    }                                                           // End loop over ranks
  }

protected:
  //! Get distance to other domain
  real_t getDistance(C_iter C, Bounds bounds, vec3 Xperiodic) {
    vec3 dX;                                                    // Distance vector
    for (int d=0; d<3; d++) {                                   // Loop over dimensions
      dX[d] = (C->X[d] + Xperiodic[d] > bounds.Xmax[d]) *       //  Calculate the distance between cell C and
	(C->X[d] + Xperiodic[d] - bounds.Xmax[d]) +             //  the nearest point in domain [xmin,xmax]^3
	(C->X[d] + Xperiodic[d] < bounds.Xmin[d]) *             //  Take the differnece from xmin or xmax
	(C->X[d] + Xperiodic[d] - bounds.Xmin[d]);              //  or 0 if between xmin and xmax
    }                                                           // End loop over dimensions
    return norm(dX);                                            // Return distance squared
  }

  //! Add cells to send buffer
  void addSendCell(C_iter C, int & irank, int & icell, int & iparent, bool copyData) {
    if (copyData) {                                             // If copying data to send cells
      Cell cell(*C);                                            //  Initialize send cell
      cell.NCHILD = cell.NBODY = 0;                             //  Reset counters
      cell.IPARENT = iparent;                                   //  Index of parent
      sendCells[sendCellDispl[irank]+icell] = cell;             //  Copy cell to send buffer
      C_iter Cparent = sendCells.begin() + sendCellDispl[irank] + iparent;// Get parent iterator
      if (Cparent->NCHILD == 0) Cparent->ICHILD = icell;        //  Index of parent's first child
      Cparent->NCHILD++;                                        //  Increment parent's child counter
    }                                                           // End if for copying data to send cells
    icell++;                                                    // Increment cell counter
  }

  //! Add bodies to send buffer
  void addSendBody(C_iter C, int & irank, int & ibody, int icell, bool copyData) {
    if (copyData) {                                             // If copying data to send bodies
      C_iter Csend = sendCells.begin() + sendCellDispl[irank] + icell; // Send cell iterator
      Csend->NBODY = C->NBODY;                                  //  Number of bodies
      Csend->IBODY = ibody;                                     //  Body index per rank
      B_iter Bsend = sendBodies.begin() + sendBodyDispl[irank] + ibody; // Send body iterator
      for (B_iter B=C->BODY; B!=C->BODY+C->NBODY; B++,Bsend++) {//  Loop over bodies in cell
	*Bsend = *B;                                            //   Copy body to send buffer
	Bsend->IRANK = irank;                                   //   Assign destination rank
      }                                                         //  End loop over bodies in cell
    }                                                           // End if for copying data to send bodies
    ibody += C->NBODY;                                          // Increment body counter
  }

  //! Determine which cells to send
  void traverseLET(C_iter C, C_iter C0, Bounds bounds, real_t cycle,
		   int & irank, int & ibody, int & icell, int iparent, bool copyData) {
    int level = int(logf(mpisize-1) / M_LN2 / 3) + 1;           // Level of local root cell
    if (mpisize == 1) level = 0;                                // Account for serial case
    bool divide[8] = {0, 0, 0, 0, 0, 0, 0, 0};                  // Initialize divide flag
    int icells[8] = {0, 0, 0, 0, 0, 0, 0, 0};                   // Initialize icell array
    int cc = 0;                                                 // Initialize child index
    for (C_iter CC=C0+C->ICHILD; CC!=C0+C->ICHILD+C->NCHILD; CC++,cc++) { // Loop over child cells
      icells[cc] = icell;                                       //  Store cell index
      addSendCell(CC, irank, icell, iparent, copyData);         //  Add cells to send
      if (CC->NCHILD == 0) {                                    //  If cell is leaf
	addSendBody(CC, irank, ibody, icell-1, copyData);       //   Add bodies to send
      } else {                                                  //  If cell is not leaf
	vec3 Xperiodic = 0;                                     //   Periodic coordinate offset
	if (images == 0) {                                      //   If free boundary condition
	  real_t R2 = getDistance(CC, bounds, Xperiodic);       //    Get distance to other domain
	  divide[cc] |= 4 * CC->R * CC->R > R2;                 //    Divide if the cell seems too close
	} else {                                                //   If periodic boundary condition
	  for (int ix=-1; ix<=1; ix++) {                        //    Loop over x periodic direction
	    for (int iy=-1; iy<=1; iy++) {                      //     Loop over y periodic direction
	      for (int iz=-1; iz<=1; iz++) {                    //      Loop over z periodic direction
		Xperiodic[0] = ix * cycle;                      //       Coordinate offset for x periodic direction
		Xperiodic[1] = iy * cycle;                      //       Coordinate offset for y periodic direction
		Xperiodic[2] = iz * cycle;                      //       Coordinate offset for z periodic direction
		real_t R2 = getDistance(CC, bounds, Xperiodic); //       Get distance to other domain
		divide[cc] |= 4 * CC->R * CC->R > R2;           //       Divide if cell seems too close
	      }                                                 //      End loop over z periodic direction
	    }                                                   //     End loop over y periodic direction
	  }                                                     //    End loop over x periodic direction
	}                                                       //   Endif for periodic boundary condition
	divide[cc] |= CC->R > (cycle / (1 << (level+1)));       //   Divide if cell is larger than local root cell
      }                                                         //  Endif for leaf
    }                                                           // End loop over child cells
    cc = 0;                                                     // Initialize child index
    for (C_iter CC=C0+C->ICHILD; CC!=C0+C->ICHILD+C->NCHILD; CC++,cc++) { // Loop over child cells
      if (divide[cc]) {                                         //  If cell must be divided further
	iparent = icells[cc];                                   //   Parent cell index
	traverseLET(CC, C0, bounds, cycle, irank, ibody, icell, iparent, copyData);// Recursively traverse tree to set LET
      }                                                         //  End if for cell division
    }                                                           // End loop over child cells
  }

public:
  //! Constructor
  TreeMPI(int _mpirank, int _mpisize, int _images) :
    mpirank(_mpirank), mpisize(_mpisize), images(_images) {     // Initialize variables
    allBoundsXmin = new float [mpisize][3];                     // Allocate array for minimum of local domains
    allBoundsXmax = new float [mpisize][3];                     // Allocate array for maximum of local domains
    sendBodyCount = new int [mpisize];                          // Allocate send count
    sendBodyDispl = new int [mpisize];                          // Allocate send displacement
    recvBodyCount = new int [mpisize];                          // Allocate receive count
    recvBodyDispl = new int [mpisize];                          // Allocate receive displacement
    sendCellCount = new int [mpisize];                          // Allocate send count
    sendCellDispl = new int [mpisize];                          // Allocate send displacement
    recvCellCount = new int [mpisize];                          // Allocate receive count
    recvCellDispl = new int [mpisize];                          // Allocate receive displacement
  }
  //! Destructor
  ~TreeMPI() {
    delete[] allBoundsXmin;                                     // Deallocate array for minimum of local domains
    delete[] allBoundsXmax;                                     // Deallocate array for maximum of local domains
    delete[] sendBodyCount;                                     // Deallocate send count
    delete[] sendBodyDispl;                                     // Deallocate send displacement
    delete[] recvBodyCount;                                     // Deallocate receive count
    delete[] recvBodyDispl;                                     // Deallocate receive displacement
    delete[] sendCellCount;                                     // Deallocate send count
    delete[] sendCellDispl;                                     // Deallocate send displacement
    delete[] recvCellCount;                                     // Deallocate receive count
    delete[] recvCellDispl;                                     // Deallocate receive displacement
  }

  //! Allgather bounds from all ranks
  void allgatherBounds(Bounds bounds) {
    float Xmin[3], Xmax[3];
    for (int d=0; d<3; d++) {                                   // Loop over dimensions
      Xmin[d] = bounds.Xmin[d];                                 //  Convert Xmin to float
      Xmax[d] = bounds.Xmax[d];                                 //  Convert Xmax to float
    }                                                           // End loop over dimensions
    MPI_Allgather(Xmin, 3, MPI_FLOAT, allBoundsXmin[0], 3, MPI_FLOAT, MPI_COMM_WORLD);// Gather all domain bounds
    MPI_Allgather(Xmax, 3, MPI_FLOAT, allBoundsXmax[0], 3, MPI_FLOAT, MPI_COMM_WORLD);// Gather all domain bounds
  }

  //! Set local essential tree to send to each process
  void setLET(Cells & cells, real_t cycle) {
    logger::startTimer("Set LET size");                         // Start timer
    C_iter C0 = cells.begin();                                  // Set cells begin iterator
    Bounds bounds;                                              // Bounds of local subdomain
    sendBodyDispl[0] = 0;                                       // Initialize body displacement vector
    sendCellDispl[0] = 0;                                       // Initialize cell displacement vector
    for (int irank=0; irank<mpisize; irank++) {                 // Loop over ranks
      if (irank != 0) sendBodyDispl[irank] = sendBodyDispl[irank-1] + sendBodyCount[irank-1];// Update body displacement
      if (irank != 0) sendCellDispl[irank] = sendCellDispl[irank-1] + sendCellCount[irank-1];// Update cell displacement
      sendBodyCount[irank] = 0;                                 //  Initialize send body count for current rank
      sendCellCount[irank] = 0;                                 //  Initialize send cell count for current rank
      if (irank != mpirank && !cells.empty()) {                 //  If not current rank and cell vector is not empty
	int ibody = 0;                                          //   Initialize send body's offset
	int icell = 1;                                          //   Initialize send cell's offset
	for (int d=0; d<3; d++) {                               //   Loop over dimensions
	  bounds.Xmin[d] = allBoundsXmin[irank][d];             //    Local Xmin for irank
	  bounds.Xmax[d] = allBoundsXmax[irank][d];             //    Local Xmax for irank
	}                                                       //   End loop over dimensions
	traverseLET(C0, C0, bounds, cycle, irank, ibody, icell, 0, false); // Traverse tree to set LET
	sendBodyCount[irank] = ibody;                           //   Send body count for current rank
	sendCellCount[irank] = icell;                           //   Send cell count for current rank
      }                                                         //  Endif for current rank
    }                                                           // End loop over ranks
    logger::stopTimer("Set LET size");                          // Stop timer
    logger::startTimer("Set LET");                              // Start timer
    int numSendBodies = sendBodyDispl[mpisize-1] + sendBodyCount[mpisize-1];// Total number of send bodies
    int numSendCells = sendCellDispl[mpisize-1] + sendCellCount[mpisize-1];// Total number of send cells
    sendBodies.resize(numSendBodies);                           // Clear send buffer for bodies
    sendCells.resize(numSendCells);                             // Clear send buffer for cells
    for (int irank=0; irank<mpisize; irank++) {                 // Loop over ranks
      if (irank != mpirank && !cells.empty()) {                 //  If not current rank and cell vector is not empty
	int ibody = 0;                                          //   Reinitialize send body's offset
	int icell = 0;                                          //   Reinitialize send cell's offset
	for (int d=0; d<3; d++) {                               //   Loop over dimensions
	  bounds.Xmin[d] = allBoundsXmin[irank][d];             //   Local Xmin for irank
	  bounds.Xmax[d] = allBoundsXmax[irank][d];             //   Local Xmax for irank
	}                                                       //   End loop over dimensions
	C_iter Csend = sendCells.begin() + sendCellDispl[irank];//   Send cell iterator
	*Csend = *C0;                                           //   Copy cell to send buffer
	Csend->NCHILD = Csend->NBODY = 0;                       //   Reset link to children and bodies
	icell++;                                                //   Increment send cell counter
	if (C0->NCHILD == 0) {                                  //   If root cell is leaf
	  addSendBody(C0, irank, ibody, icell-1, true);         //    Add bodies to send
	}                                                       //   End if for root cell leaf
        traverseLET(C0, C0, bounds, cycle, irank, ibody, icell, 0, true); // Traverse tree to set LET
      }                                                         //  Endif for current rank
    }                                                           // End loop over ranks
    logger::stopTimer("Set LET");                               // Stop timer
  }

  //! Get local essential tree from irank
  void getLET(Cells & cells, int irank) {
    std::stringstream event;                                    // Event name
    event << "Get LET from rank " << irank;                     // Create event name based on irank
    logger::startTimer(event.str());                            // Start timer
    for (int i=0; i<recvCellCount[irank]; i++) {                // Loop over receive cells
      C_iter C = recvCells.begin() + recvCellDispl[irank] + i;  //  Iterator of receive cell
      if (C->NBODY != 0) {                                      //  If cell has bodies
        C->BODY = recvBodies.begin() + recvBodyDispl[irank] + C->IBODY;// Iterator of first body
      }                                                         //  End if for bodies
    }                                                           // End loop over receive cells
    cells.resize(recvCellCount[irank]);                         // Resize cell vector for LET
    cells.assign(recvCells.begin()+recvCellDispl[irank],        // Assign receive cells to vector
		 recvCells.begin()+recvCellDispl[irank]+recvCellCount[irank]);
    logger::stopTimer(event.str());                             // Stop timer
  }

  //! Link LET with received bodies and calcualte NBODY
  void linkLET() {
    logger::startTimer("Link LET");                             // Start timer
    for (int irank=0; irank<mpisize; irank++) {                 // Loop over ranks
      for (int i=0; i<recvCellCount[irank]; i++) {              //  Loop over receive cells
        C_iter C = recvCells.begin() + recvCellDispl[irank] + i;//   Iterator of receive cell
        if (C->NBODY != 0) {                                    //   If cell has bodies
          C->BODY = recvBodies.begin() + recvBodyDispl[irank] + C->IBODY;// Iterator of first body
        }                                                       //   End if for bodies
      }                                                         //  End loop over receive cells
    }                                                           // End loop over ranks
    logger::stopTimer("Link LET");                              // End timer
  }

  //! Copy remote root cells to body structs (for building global tree)
  Bodies root2body() {
    logger::startTimer("Root to body");                         // Start timer
    Bodies bodies;                                              // Bodies to contain remote root coordinates
    bodies.reserve(mpisize-1);                                  // Reserve size of body vector
    for (int irank=0; irank<mpisize; irank++) {                 // Loop over ranks
      if (irank != mpirank) {                                   //  If not current rank
	C_iter C0 = recvCells.begin() + recvCellDispl[irank];   //   Root cell iterator for irank
	Body body;                                              //   Body to contain remote root coordinates
	body.X = C0->X;                                         //   Copy remote root coordinates
	body.IBODY = recvCellDispl[irank];                      //   Copy remote root displacement in vector
	bodies.push_back(body);                                 //   Push this root cell to body vector
      }                                                         //  End if for not current rank
    }                                                           // End loop over ranks
    logger::stopTimer("Root to body");                          // Stop timer
    return bodies;                                              // Return body vector
  }

  //! Graft remote trees to global tree
  void attachRoot(Cells & cells) {
    logger::startTimer("Attach root");                          // Start timer
    int globalCells = cells.size();                             // Number of global cells
    cells.insert(cells.end(), recvCells.begin(), recvCells.end()); // Join LET cell vectors
    for (C_iter C=cells.begin(); C!=cells.begin()+globalCells; C++) { // Loop over global cells
      if (C->NCHILD==0) {                                       // If leaf cell
	int offset = globalCells + C->BODY->IBODY;              //  Offset of received root cell index
	C_iter C0 = cells.begin() + offset;                     //  Root cell iterator
	C0->IPARENT = C->IPARENT;                               //  Link remote root to global leaf
	*C = *C0;                                               //  Copy remote root to global leaf
	C->ICHILD += offset;                                    //  Add offset to child index
      } else {                                                  // If not leaf cell
	C->BODY = recvBodies.end();                             //  Use BODY as flag to indicate non-leaf global cell
      }                                                         // End if for leaf cell
    }                                                           // End loop over global cells
    for (int irank=0; irank<mpisize; irank++) {                 // Loop over ranks
      if (irank != mpirank && recvCellCount[irank] > 0) {       //  If not current rank
	C_iter C0 = cells.begin() + globalCells + recvCellDispl[irank];// Root cell iterator for irank
	for (C_iter C=C0+1; C!=C0+recvCellCount[irank]; C++) {  //    Loop over cells received from irank
	  int offset = globalCells + recvCellDispl[irank];      //     Offset of received root cell index
	  C->IPARENT += offset;                                 //     Add offset to parent index
	  C->ICHILD += offset;                                  //     Add offset to child index
	}                                                       //    End loop over cells received from irank
      }                                                         //  End if for not current rank
    }                                                           // End loop over ranks
    C_iter C0 = cells.begin();                                  // Root cell iterator of entire LET
    for (int i=globalCells-1; i>=0; i--) {                      // Loop over global cells bottom up
      C_iter C = cells.begin() + i;                             //  Iterator of current cell
      if (C->BODY == recvBodies.end()) {                        //  If non-leaf global cell
	vec3 Xmin = C->X, Xmax = C->X;                          //   Initialize Xmin, Xmax
	for (C_iter CC=C0+C->ICHILD; CC!=C0+C->ICHILD+C->NCHILD; CC++) { // Loop over child cells
	  Xmin = min(CC->X-CC->R, Xmin);                        //    Update Xmin
          Xmax = max(CC->X+CC->R, Xmax);                        //    Update Xmax
	}                                                       //   End loop over child cells
	C->X = (Xmax + Xmin) / 2;                               //   Calculate center of domain
	for (int d=0; d<3; d++) {                               //   Loop over dimensions
	  C->R = std::max(C->X[d] - Xmin[d], C->R);             //    Calculate min distance from center
	  C->R = std::max(Xmax[d] - C->X[d], C->R);             //    Calculate max distance from center
	}                                                       //   End loop over dimensions
	C->M = 0;                                               //   Reset multipoles
	kernel::M2M(C, C0);                                     //   M2M kernel
      }                                                         //  End if for non-leaf global cell
    }                                                           // End loop over global cells bottom up
    logger::stopTimer("Attach root");                           // Stop timer
  }

  //! Send bodies
  Bodies commBodies(Bodies bodies) {
    logger::startTimer("Comm partition");                       // Start timer
    alltoall(bodies);                                           // Send body count
    alltoallv(bodies);                                          // Send bodies
    logger::stopTimer("Comm partition");                        // Stop timer
    return recvBodies;                                          // Return received bodies
  }

  //! Send bodies
  Bodies commBodies() {
    logger::startTimer("Comm LET bodies");                      // Start timer
    alltoall(sendBodies);                                       // Send body count
    alltoallv(sendBodies);                                      // Send bodies
    logger::stopTimer("Comm LET bodies");                       // Stop timer
    return recvBodies;                                          // Return received bodies
  }

  //! Send cells
  void commCells() {
    logger::startTimer("Comm LET cells");                       // Start timer
    alltoall(sendCells);                                        // Send cell count
    alltoallv(sendCells);                                       // Senc cells
    logger::stopTimer("Comm LET cells");                        // Stop timer
  }

  //! Copy recvBodies to bodies
  Bodies getRecvBodies() {
    return recvBodies;                                          // Return recvBodies
  }

  //! Send bodies to next rank (round robin)
  void shiftBodies(Bodies & bodies) {
    int newSize;                                                // New number of bodies
    int oldSize = bodies.size();                                // Current number of bodies
    const int word = sizeof(bodies[0]) / 4;                     // Word size of body structure
    const int isend = (mpirank + 1          ) % mpisize;        // Send to next rank (wrap around)
    const int irecv = (mpirank - 1 + mpisize) % mpisize;        // Receive from previous rank (wrap around)
    MPI_Request sreq,rreq;                                      // Send, receive request handles

    MPI_Isend(&oldSize, 1, MPI_INT, irecv, 0, MPI_COMM_WORLD, &sreq);// Send current number of bodies
    MPI_Irecv(&newSize, 1, MPI_INT, isend, 0, MPI_COMM_WORLD, &rreq);// Receive new number of bodies
    MPI_Wait(&sreq, MPI_STATUS_IGNORE);                         // Wait for send to complete
    MPI_Wait(&rreq, MPI_STATUS_IGNORE);                         // Wait for receive to complete

    recvBodies.resize(newSize);                                 // Resize buffer to new number of bodies
    MPI_Isend((int*)&bodies[0], oldSize*word, MPI_INT, irecv,   // Send bodies to next rank
              1, MPI_COMM_WORLD, &sreq);
    MPI_Irecv((int*)&recvBodies[0], newSize*word, MPI_INT, isend,// Receive bodies from previous rank
              1, MPI_COMM_WORLD, &rreq);
    MPI_Wait(&sreq, MPI_STATUS_IGNORE);                         // Wait for send to complete
    MPI_Wait(&rreq, MPI_STATUS_IGNORE);                         // Wait for receive to complete
    bodies = recvBodies;                                        // Copy bodies from buffer
  }

  //! Allgather bodies
  Bodies allgatherBodies(Bodies & bodies) {
    const int word = sizeof(bodies[0]) / 4;                     // Word size of body structure
    sendBodyCount[0] = bodies.size();                           // Determine send count
    MPI_Allgather(sendBodyCount, 1, MPI_INT,                    // Allgather number of bodies
                  recvBodyCount, 1, MPI_INT, MPI_COMM_WORLD);
    recvBodyDispl[0] = 0;                                       // Initialize receive displacement
    for (int irank=0; irank<mpisize-1; irank++) {               // Loop over ranks
      recvBodyDispl[irank+1] = recvBodyDispl[irank] + recvBodyCount[irank];// Set receive displacement
    }                                                           // End loop over ranks
    recvBodies.resize(recvBodyDispl[mpisize-1]+recvBodyCount[mpisize-1]);// Resize receive buffer
    for (int irank=0; irank<mpisize; irank++) {                 // Loop over ranks
      recvBodyCount[irank] *= word;                             //  Multiply receive count by word size of data
      recvBodyDispl[irank] *= word;                             //  Multiply receive displacement by word size of data
    }                                                           // End loop over ranks
    MPI_Allgatherv((int*)&bodies[0], sendBodyCount[0]*word, MPI_INT,// Allgather bodies
                   (int*)&recvBodies[0], recvBodyCount, recvBodyDispl, MPI_INT, MPI_COMM_WORLD);
    return recvBodies;                                          // Return bodies
  }
};
#endif
