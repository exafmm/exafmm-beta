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
#ifndef partition_h
#define partition_h
#include "../include/mympi.h"
#include "serialfmm.h"

//! Handles all the partitioning of domains
template<Equation equation>
class Partition : public MyMPI, public SerialFMM<equation> {
public:
  using Kernel<equation>::printNow;                             //!< Switch to print timings
  using Kernel<equation>::startTimer;                           //!< Start timer for given event
  using Kernel<equation>::stopTimer;                            //!< Stop timer for given event
  using Kernel<equation>::sortBodies;                           //!< Sort bodies according to cell index
  using Kernel<equation>::setX0;                                //!< Set center of root cell
  using Kernel<equation>::setR0;                                //!< Set radius of root cell
  using Evaluator<equation>::MAXLEVEL;                          //!< Max depth of tree

public:
//! Send bodies to next rank (round robin)
  void shiftBodies(Bodies &bodies) {
    int newSize;                                                // New number of bodies
    int oldSize = bodies.size();                                // Current number of bodies
    const int bytes = sizeof(bodies[0]);                        // Byte size of body structure
    const int isend = (MPIRANK + 1          ) % MPISIZE;        // Send to next rank (wrap around)
    const int irecv = (MPIRANK - 1 + MPISIZE) % MPISIZE;        // Receive from previous rank (wrap around)
    MPI_Request sreq,rreq;                                      // Send, recv request handles

    MPI_Isend(&oldSize,1,MPI_INT,irecv,0,MPI_COMM_WORLD,&sreq); // Send current number of bodies
    MPI_Irecv(&newSize,1,MPI_INT,isend,0,MPI_COMM_WORLD,&rreq); // Receive new number of bodies
    MPI_Wait(&sreq,MPI_STATUS_IGNORE);                          // Wait for send to complete
    MPI_Wait(&rreq,MPI_STATUS_IGNORE);                          // Wait for recv to complete

    Bodies buffer(newSize);                                     // Resize buffer to new number of bodies
    MPI_Isend(&bodies[0],oldSize*bytes,MPI_BYTE,irecv,          // Send bodies to next rank
              1,MPI_COMM_WORLD,&sreq);
    MPI_Irecv(&buffer[0],newSize*bytes,MPI_BYTE,isend,          // Receive bodies from previous rank
              1,MPI_COMM_WORLD,&rreq);
    MPI_Wait(&sreq,MPI_STATUS_IGNORE);                          // Wait for send to complete
    MPI_Wait(&rreq,MPI_STATUS_IGNORE);                          // Wait for recv to complete
    bodies = buffer;                                            // Copy bodies from buffer
  }

//! Partition by recursive octsection
  void octsection(Bodies &bodies) {
    startTimer("Partition");                                    // Start timer
    int byte = sizeof(bodies[0]);                               // Byte size of body structure
    MAXLEVEL = int(log(MPISIZE-1) / M_LN2 / 3) + 1;             // Level of local root cell
    if( MPISIZE == 1 ) MAXLEVEL = 0;                            // For serial execution local root cell is root cell
    setX0(0);
    setR0(M_PI);
    BottomUp<equation>::setIndex(bodies);                       // Set index of bodies for that level
    Bodies buffer = bodies;
    stopTimer("Partition");                                     // Stop timer 
    sortBodies(bodies,buffer);                                  // Sort bodies in ascending order
    startTimer("Partition");                                    // Start timer
    int *scnt = new int [MPISIZE];                              // Send count
    int *sdsp = new int [MPISIZE];                              // Send displacement
    int *rcnt = new int [MPISIZE];                              // Recv count
    int *rdsp = new int [MPISIZE];                              // Recv displacement
    for( int i=0; i!=MPISIZE; ++i ) {                           // Loop over ranks
      scnt[i] = 0;                                              //  Initialize send counts
    }                                                           // End loop over ranks
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      int irank = B->ICELL / (int(pow(8,MAXLEVEL)) / MPISIZE);  //  Get rank which the cell belongs to
      scnt[irank]++;                                            //  Fill send count bucket
    }                                                           // End loop over bodies
    MPI_Alltoall(scnt,1,MPI_INT,rcnt,1,MPI_INT,MPI_COMM_WORLD); // Communicate send count to get recv count
    sdsp[0] = rdsp[0] = 0;                                      // Initialize send/recv displacements
    for( int irank=0; irank!=MPISIZE-1; ++irank ) {             // Loop over ranks
      sdsp[irank+1] = sdsp[irank] + scnt[irank];                //  Set send displacement based on send count
      rdsp[irank+1] = rdsp[irank] + rcnt[irank];                //  Set recv displacement based on recv count
    }                                                           // End loop over ranks
    buffer.resize(rdsp[MPISIZE-1]+rcnt[MPISIZE-1]);             // Resize recv buffer
    for( int irank=0; irank!=MPISIZE; ++irank ) {               // Loop over ranks
      scnt[irank] *= byte;                                      //  Multiply send count by byte size of data
      sdsp[irank] *= byte;                                      //  Multiply send displacement by byte size of data
      rcnt[irank] *= byte;                                      //  Multiply recv count by byte size of data
      rdsp[irank] *= byte;                                      //  Multiply recv displacement by byte size of data
    }                                                           // End loop over ranks
    MPI_Alltoallv(&bodies[0],scnt,sdsp,MPI_BYTE,&buffer[0],rcnt,rdsp,MPI_BYTE,MPI_COMM_WORLD);// Communicat bodies
    bodies = buffer;                                            // Copy recv buffer to bodies
    delete[] scnt;                                              // Delete send count
    delete[] sdsp;                                              // Delete send displacement
    delete[] rcnt;                                              // Delete recv count
    delete[] rdsp;                                              // Delete recv displacement
    stopTimer("Partition",printNow);                            // Stop timer 
  }

//! Send bodies back to where they came from
  void unpartition(Bodies &bodies) {
    startTimer("Unpartition");                                  // Start timer
    int byte = sizeof(bodies[0]);                               // Byte size of body structure
    int *scnt = new int [MPISIZE];                              // Send count
    int *sdsp = new int [MPISIZE];                              // Send displacement
    int *rcnt = new int [MPISIZE];                              // Recv count
    int *rdsp = new int [MPISIZE];                              // Recv displacement
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      B->ICELL = B->IPROC;                                      //  Copy process rank to cell index for sorting
    }                                                           // End loop over bodies
    Bodies buffer = bodies;                                     // Resize sort buffer
    stopTimer("Unpartition");                                   // Stop timer 
    sortBodies(bodies,buffer);                                  // Sort bodies in ascending order
    startTimer("Unpartition");                                  // Start timer
    for( int i=0; i!=MPISIZE; ++i ) {                           // Loop over ranks
      scnt[i] = 0;                                              //  Initialize send counts
    }                                                           // End loop over ranks
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      int irank = B->IPROC;                                     //  Get rank which the body belongs to
      scnt[irank]++;                                            //  Fill send count bucket
    }                                                           // End loop over bodies
    MPI_Alltoall(scnt,1,MPI_INT,rcnt,1,MPI_INT,MPI_COMM_WORLD); // Communicate send count to get recv count
    sdsp[0] = rdsp[0] = 0;                                      // Initialize send/recv displacements
    for( int irank=0; irank!=MPISIZE-1; ++irank ) {             // Loop over ranks
      sdsp[irank+1] = sdsp[irank] + scnt[irank];                //  Set send displacement based on send count
      rdsp[irank+1] = rdsp[irank] + rcnt[irank];                //  Set recv displacement based on recv count
    }                                                           // End loop over ranks
    buffer.resize(rdsp[MPISIZE-1]+rcnt[MPISIZE-1]);             // Resize recv buffer
    for( int irank=0; irank!=MPISIZE; ++irank ) {               // Loop over ranks
      scnt[irank] *= byte;                                      //  Multiply send count by byte size of data
      sdsp[irank] *= byte;                                      //  Multiply send displacement by byte size of data
      rcnt[irank] *= byte;                                      //  Multiply recv count by byte size of data
      rdsp[irank] *= byte;                                      //  Multiply recv displacement by byte size of data
    }                                                           // End loop over ranks
    MPI_Alltoallv(&bodies[0],scnt,sdsp,MPI_BYTE,&buffer[0],rcnt,rdsp,MPI_BYTE,MPI_COMM_WORLD);// Communicat bodies
    bodies = buffer;                                            // Copy recv buffer to bodies
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      assert(B->IPROC == MPIRANK);                              //  Make sure bodies are in the right rank
    }                                                           // End loop over bodies
    delete[] scnt;                                              // Delete send count
    delete[] sdsp;                                              // Delete send displacement
    delete[] rcnt;                                              // Delete recv count
    delete[] rdsp;                                              // Delete recv displacement
    stopTimer("Unpartition",printNow);                          // Stop timer 
  }
};

#endif
