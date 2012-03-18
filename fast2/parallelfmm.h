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
  int IRANK;                                                    //!< MPI rank loop counter
  vect thisXMIN;                                                //!< XMIN for a given rank
  vect thisXMAX;                                                //!< XMAX for a given rank
  B_iter Bj0;                                                   //!< Source bodies begin iterator
  Cells sendCells;                                              //!< Send buffer for cells
  Cells recvCells;                                              //!< Receive buffer for cells
  int *sendCellCount;                                           //!< Send count
  int *sendCellDispl;                                           //!< Send displacement
  int *recvCellCount;                                           //!< Receive count
  int *recvCellDispl;                                           //!< Receive displacement

public:
  using Kernel<equation>::printNow;                             //!< Switch to print timings
  using Kernel<equation>::startTimer;                           //!< Start timer for given event
  using Kernel<equation>::stopTimer;                            //!< Stop timer for given event
  using Kernel<equation>::X0;                                   //!< Center of root cell
  using Kernel<equation>::R0;                                   //!< Radius of root cel
  using Kernel<equation>::Cj0;                                  //!< Begin iterator for source cells
  using Evaluator<equation>::TOPDOWN;                           //!< Flag for top down tree construction
  using Partition<equation>::print;                             //!< Print in MPI
  using Partition<equation>::XMIN;                              //!< Minimum position vector of bodies
  using Partition<equation>::XMAX;                              //!< Maximum position vector of bodies
  using Partition<equation>::sendBodies;                        //!< Send buffer for bodies
  using Partition<equation>::recvBodies;                        //!< Receive buffer for bodies
  using Partition<equation>::sendBodyCount;                     //!< Send count
  using Partition<equation>::sendBodyDispl;                     //!< Send displacement
  using Partition<equation>::recvBodyCount;                     //!< Receive count
  using Partition<equation>::recvBodyDispl;                     //!< Receive displacement
  using Partition<equation>::alltoall;                          //!< Exchange send count for bodies
  using Partition<equation>::alltoallv;                         //!< Exchange bodies

private:
//! Get disatnce to other domain
  real getDistance(C_iter C) {
    vect dist;                                                  // Distance vector
    for( int d=0; d!=3; ++d ) {                                 // Loop over dimensions
      dist[d] = (C->X[d] + Xperiodic[d] > thisXMAX[d])*         //  Calculate the distance between cell C and
                (C->X[d] + Xperiodic[d] - thisXMAX[d])+         //  the nearest point in domain [xmin,xmax]^3
                (C->X[d] + Xperiodic[d] < thisXMIN[d])*         //  Take the differnece from xmin or xmax
                (C->X[d] + Xperiodic[d] - thisXMIN[d]);         //  or 0 if between xmin and xmax
    }                                                           // End loop over dimensions
    real R = std::sqrt(norm(dist));                             // Scalar distance
    return R;                                                   // Return scalar distance
  }

//! Add cells to send buffer
  void addSendCell(C_iter C, int &icell) {
    Cell cell(*C);
    cell.NCHILD = cell.NCLEAF = cell.NDLEAF = 0;
    cell.PARENT = icell - sendCellDispl[IRANK];
    sendCells.push_back(cell);
    C_iter Cparent = sendCells.begin() + icell;
    icell = sendCells.size() - 1;
    if( Cparent->NCHILD == 0 ) Cparent->CHILD = icell - sendCellDispl[IRANK];
    Cparent->NCHILD++;
  }

//! Add bodies to send buffer
  void addSendBody(C_iter C, int &ibody, int icell) {
    C_iter Csend = sendCells.begin() + icell;
    Csend->ILEAF = ibody;
    Csend->NCLEAF = C->NCLEAF;
    Csend->NDLEAF = C->NDLEAF;
    for( B_iter B=C->LEAF; B!=C->LEAF+C->NCLEAF; ++B ) {
      sendBodies.push_back(*B);
      sendBodies.back().IPROC = IRANK;
      ibody++;
    }
  }

//! Determine which cells to send
  void traverseLET(C_iter C, int &ibody, int icell) {
    int level = int(log(MPISIZE-1) / M_LN2 / 3) + 1;            // Level of local root cell
    if( MPISIZE == 1 ) level = 0;                               // Account for serial case
    for( C_iter CC=Cj0+C->CHILD; CC!=Cj0+C->CHILD+C->NCHILD; ++CC ) {// Loop over child cells
      addSendCell(CC,icell);                                    //  Add cells to send
      if( CC->NCHILD == 0 ) {                                   //  If cell is twig
        addSendBody(CC,ibody,icell);                            //   Add bodies to send
      } else {                                                  //  If cell is not twig
        bool divide = false;                                    //   Initialize logical for dividing
        if( IMAGES == 0 ) {                                     //   If free boundary condition
          Xperiodic = 0;                                        //    Set periodic coordinate offset
          real R = getDistance(CC);                             //    Get distance to other domain
          divide |= 2 * CC->R > THETA * R;                      //    Divide if the cell seems too close
        } else {                                                //   If periodic boundary condition
          for( int ix=-1; ix<=1; ++ix ) {                       //    Loop over x periodic direction
            for( int iy=-1; iy<=1; ++iy ) {                     //     Loop over y periodic direction
              for( int iz=-1; iz<=1; ++iz ) {                   //      Loop over z periodic direction
                Xperiodic[0] = ix * 2 * R0;                     //       Coordinate offset for x periodic direction
                Xperiodic[1] = iy * 2 * R0;                     //       Coordinate offset for y periodic direction
                Xperiodic[2] = iz * 2 * R0;                     //       Coordinate offset for z periodic direction
                real R = getDistance(CC);                       //       Get distance to other domain
                divide |= 2 * CC->R > THETA * R;                //       Divide if cell seems too close
              }                                                 //      End loop over z periodic direction
            }                                                   //     End loop over y periodic direction
          }                                                     //    End loop over x periodic direction
        }                                                       //   Endif for periodic boundary condition
        divide |= CC->R > (R0 / (1 << level));                  //   Divide if cell is larger than local root cell
        if( divide ) {                                          //   If cell has to be divided
          traverseLET(CC,ibody,icell);                          //    Traverse the tree further
        }                                                       //   Endif for cell division
      }                                                         //  Endif for twig
    }                                                           // End loop over child cells
  }

//! Exchange send count for cells
  void alltoall(Cells &cells) {
    MPI_Alltoall(sendCellCount,1,MPI_INT,                       // Communicate send count to get receive count
                 recvCellCount,1,MPI_INT,MPI_COMM_WORLD);
    recvCellDispl[0] = 0;                                       // Initialize receive displacements
    for( int irank=0; irank!=MPISIZE-1; ++irank ) {             // Loop over ranks
      recvCellDispl[irank+1] = recvCellDispl[irank] + recvCellCount[irank];//  Set receive displacement
    }                                                           // End loop over ranks
  }

//! Exchange cells
  void alltoallv(Cells &cells) {
    int word = sizeof(cells[0]) / 4;                           // Word size of body structure
    recvCells.resize(recvCellDispl[MPISIZE-1]+recvCellCount[MPISIZE-1]);// Resize receive buffer
    for( int irank=0; irank!=MPISIZE; ++irank ) {               // Loop over ranks
      sendCellCount[irank] *= word;                             //  Multiply send count by word size of data
      sendCellDispl[irank] *= word;                             //  Multiply send displacement by word size of data
      recvCellCount[irank] *= word;                             //  Multiply receive count by word size of data
      recvCellDispl[irank] *= word;                             //  Multiply receive displacement by word size of data
    }                                                           // End loop over ranks
    MPI_Alltoallv(&cells[0],sendCellCount,sendCellDispl,MPI_INT,// Communicate cells
                  &recvCells[0],recvCellCount,recvCellDispl,MPI_INT,MPI_COMM_WORLD);
  }

public:
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

//! Get local essential tree to send to each process
  void getLET(Bodies &bodies, Cells &cells) {
    startTimer("Get LET");                                      // Start timer
    Bj0 = bodies.begin();                                       // Set bodies begin iterator
    Cj0 = cells.begin();                                        // Set cells begin iterator
    C_iter Root;                                                // Root cell
    if( TOPDOWN ) {                                             // If tree was constructed top down
      Root = cells.begin();                                     //  The first cell is the root cell
    } else {                                                    // If tree was constructed bottom up
      Root = cells.end() - 1;                                   //  The last cell is the root cell
    }                                                           // Endif for tree construction
    sendCells.reserve(cells.size()*27);                         // Reserve space for send buffer
    sendCellDispl[0] = 0;                                       // Initialize displacement vector
    for( IRANK=0; IRANK!=MPISIZE; ++IRANK ) {                   // Loop over ranks
      if( IRANK != 0 ) sendCellDispl[IRANK] = sendCellDispl[IRANK-1] + sendCellCount[IRANK-1];// Update displacement
      if( IRANK != MPIRANK ) {                                  //  If not current rank
        thisXMIN = XMIN[IRANK];                                 //   Set XMIN for IRANK
        thisXMAX = XMAX[IRANK];                                 //   Set XMAX for IRANK
        Cell cell(*Root);                                       //   Send root cell
        cell.NCHILD = cell.NCLEAF = cell.NDLEAF = 0;            //   Reset link to children and leafs
        sendCells.push_back(cell);                              //   Push it into send buffer
        int ibody = sendBodyDispl[IRANK];                       //   Current send bodies offset
        int icell = sendCellDispl[IRANK];                       //   Current send cells offset
        traverseLET(Root,ibody,icell);                          //   Traverse tree to get LET
      }                                                         //  Endif for current rank
      sendCellCount[IRANK] = sendCells.size() - sendCellDispl[IRANK];// Send count for IRANK
    }                                                           // End loop over ranks
    stopTimer("Get LET",printNow);                              // Stop timer
    std::cout << MPIRANK << " " << cells.size() << " " << sendCells.size() << std::endl;
    std::cout << MPIRANK << " " << bodies.size() << " " << sendBodies.size() << std::endl;
  }

//! Send bodies
  void commBodies() {
    startTimer("Comm bodies");                                  // Start timer
    alltoall(sendBodies);
    alltoallv(sendBodies);
    stopTimer("Comm bodies",printNow);                          // Stop timer
  }

//! Send cells
  void commCells() {

  }

};
#endif
