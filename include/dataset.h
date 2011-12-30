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
#ifndef dataset_h
#define dataset_h
#include "kernel.h"

//! Contains all the different datasets
template<Equation equation>
class Dataset : public Kernel<equation> {
private:
  long filePosition;                                            //!< Position of file stream

public:
//! Constructor
  Dataset() : filePosition(0) {}
//! Destructor
  ~Dataset() {}

//! Initialize source values
  void initSource(Bodies &bodies) {
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      B->IBODY = B-bodies.begin();                              //  Tag body with initial index
      B->IPROC = MPIRANK;                                       //  Tag body with initial MPI rank
      B->SRC = 1. / bodies.size() / MPISIZE;                    //  Initialize mass/charge
    }                                                           // End loop over bodies
  }

//! Initialize target values
  void initTarget(Bodies &bodies, bool IeqJ=true) {
    srand48(0);                                                 // Set seed for random number generator
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      B->IBODY = B-bodies.begin();                              //  Tag body with initial index
      B->IPROC = MPIRANK;                                       //  Tag body with initial MPI rank
      B->TRG = 0;                                               //  Clear previous target values (IeqJ is dummy)
      if( EPS2 != 0 ) B->TRG[0] = -B->SRC / std::sqrt(EPS2) * IeqJ;//  Initialize potential (0 if I != J)
    }                                                           // End loop over bodies
  }

//! Read target values from file
  void readTarget(Bodies &bodies) {
    char fname[256];                                            // File name for saving direct calculation values
    sprintf(fname,"direct%4.4d",MPIRANK);                       // Set file name
    std::ifstream file(fname,std::ios::in | std::ios::binary);  // Open file
    file.seekg(filePosition);                                   // Set position in file
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      file >> B->TRG[0];                                        //  Read data for potential
      file >> B->TRG[1];                                        //  Read data for x acceleration
      file >> B->TRG[2];                                        //  Read data for y acceleration
      file >> B->TRG[3];                                        //  Read data for z acceleration
    }                                                           // End loop over bodies
    filePosition = file.tellg();                                // Get position in file
    file.close();                                               // Close file
  }

//! Write target values to file
  void writeTarget(Bodies &bodies) {
    char fname[256];                                            // File name for saving direct calculation values
    sprintf(fname,"direct%4.4d",MPIRANK);                       // Set file name
    std::ofstream file(fname,std::ios::out | std::ios::app | std::ios::binary);// Open file
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      file << B->TRG[0] << std::endl;                           //  Write data for potential
      file << B->TRG[1] << std::endl;                           //  Write data for x acceleration
      file << B->TRG[2] << std::endl;                           //  Write data for y acceleration
      file << B->TRG[3] << std::endl;                           //  Write data for z acceleration
    }                                                           // End loop over bodies
    file.close();                                               // Close file
  }

//! Evaluate relative L2 norm error
  void evalError(Bodies &bodies, Bodies &bodies2,
                 real &diff1, real &norm1, real &diff2, real &norm2, bool ewald=false) {
    real p = 0, p2 = 0;                                         // Total energy
    B_iter B2 = bodies2.begin();                                // Set iterator for bodies2
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B, ++B2 ) {// Loop over bodies & bodies2
#ifdef DEBUG
      std::cout << B->ICELL << " " << B->TRG[0] << " " << B2->TRG[0] << std::endl;// Compare every element
#endif
      if( ewald ) {                                             // If ewald method is used
        p += B->TRG[0] * B->SRC;                                //  total energy
        p2 += B2->TRG[0] * B2->SRC;                             //  total energy
      } else {                                                  // If normal Laplace kernel is used
        diff1 += (B->TRG[0] - B2->TRG[0]) * (B->TRG[0] - B2->TRG[0]);// Difference of potential
        norm1 += B2->TRG[0] * B2->TRG[0];                       //  Value of potential
      }                                                         // End if for Ewald method
      diff2 += (B->TRG[1] - B2->TRG[1]) * (B->TRG[1] - B2->TRG[1]);// Difference of x acceleration
      diff2 += (B->TRG[2] - B2->TRG[2]) * (B->TRG[2] - B2->TRG[2]);// Difference of y acceleration
      diff2 += (B->TRG[3] - B2->TRG[3]) * (B->TRG[3] - B2->TRG[3]);// Difference of z acceleration
      norm2 += B2->TRG[1] * B2->TRG[1];                         //  Value of x acceleration
      norm2 += B2->TRG[2] * B2->TRG[2];                         //  Value of y acceleration
      norm2 += B2->TRG[3] * B2->TRG[3];                         //  Value of z acceleration
    }                                                           //  End loop over bodies & bodies2
    if( ewald ) {                                               // If ewald method is used
      diff1 = (p - p2) * (p - p2);                              //  Difference of total energy
      norm1 = p2 * p2;                                          //  Value of total energy
    }                                                           // End if for Ewald method
  }

//! Print relative L2 norm error
  void printError(real diff1, real norm1, real diff2, real norm2) {
    std::cout << "Error (pot)   : " << std::sqrt(diff1/norm1) << std::endl;
    std::cout << "Error (acc)   : " << std::sqrt(diff2/norm2) << std::endl;
  }
};

template<>
class Dataset<VanDerWaals> : public Kernel<VanDerWaals> {
private:
  long filePosition;                                            //!< Position of file stream

public:
//! Constructor
  Dataset() : filePosition(0) {}
//! Destructor
  ~Dataset() {}

//! Initialize source values
  void initSource(Bodies &bodies) {
    THETA = .1;                                                 // Force opening angle to be small
    ATOMS = 16;                                                 // Set number of atoms
    RSCALE.resize(ATOMS*ATOMS);                                 // Resize rscale vector
    GSCALE.resize(ATOMS*ATOMS);                                 // Resize gscale vector
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      B->IBODY = B-bodies.begin();                              //  Tag body with initial index
      B->IPROC = MPIRANK;                                       //  Tag body with initial MPI rank
      B->SRC = drand48() * ATOMS;                               //  Initialize mass/charge
    }                                                           // End loop over bodies
    for( int i=0; i!=ATOMS; ++i ) {                             // Loop over atoms
      GSCALE[i*ATOMS+i] = drand48();                            //  Set VanDerWaals post scaling factor
      RSCALE[i*ATOMS+i] = drand48();                            //  Set VanDerWaals pre scaling factor
    }                                                           // End loop over atoms
    for( int i=0; i!=ATOMS; ++i ) {                             // Loop over target atoms
      for( int j=0; j!=ATOMS; ++j ) {                           //  Loop over source atoms
        if( i != j ) {                                          //   If target and source are different
          GSCALE[i*ATOMS+j] = std::sqrt(GSCALE[i*ATOMS+i] * GSCALE[j*ATOMS+j]);// Set post scaling factor
          RSCALE[i*ATOMS+j] = (std::sqrt(RSCALE[i*ATOMS+i]) + std::sqrt(RSCALE[j*ATOMS+j])) * 0.5;
          RSCALE[i*ATOMS+j] *= RSCALE[i*ATOMS+j];               //    Set pre scaling factor
        }                                                       //   End if for different target and source
      }                                                         //  End loop over source atoms
    }                                                           // End loop over target atoms
  }

//! Initialize target values
  void initTarget(Bodies &bodies, bool IeqJ=true) {
    srand48(0);                                                 // Set seed for random number generator
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      B->IBODY = B-bodies.begin();                              //  Tag body with initial index
      B->IPROC = MPIRANK;                                       //  Tag body with initial MPI rank
      B->TRG = 0 * IeqJ;                                        //  Clear previous target values
    }                                                           // End loop over bodies
  }

//! Read target values from file
  void readTarget(Bodies &bodies) {
    char fname[256];                                            // File name for saving direct calculation values
    sprintf(fname,"direct%4.4d",MPIRANK);                       // Set file name
    std::ifstream file(fname,std::ios::in | std::ios::binary);  // Open file
    file.seekg(filePosition);                                   // Set position in file
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      file >> B->TRG[0];                                        //  Read data for potential
      file >> B->TRG[1];                                        //  Read data for x acceleration
      file >> B->TRG[2];                                        //  Read data for y acceleration
      file >> B->TRG[3];                                        //  Read data for z acceleration
    }                                                           // End loop over bodies
    filePosition = file.tellg();                                // Get position in file
    file.close();                                               // Close file
  }

//! Write target values to file
  void writeTarget(Bodies &bodies) {
    char fname[256];                                            // File name for saving direct calculation values
    sprintf(fname,"direct%4.4d",MPIRANK);                       // Set file name
    std::ofstream file(fname,std::ios::out | std::ios::app | std::ios::binary);// Open file
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      file << B->TRG[0] << std::endl;                           //  Write data for potential
      file << B->TRG[1] << std::endl;                           //  Write data for x acceleration
      file << B->TRG[2] << std::endl;                           //  Write data for y acceleration
      file << B->TRG[3] << std::endl;                           //  Write data for z acceleration
    }                                                           // End loop over bodies
    file.close();                                               // Close file
  }

//! Evaluate relative L2 norm error
  void evalError(Bodies &bodies, Bodies &bodies2,
                 real &diff1, real &norm1, real &diff2, real &norm2) {
    B_iter B2 = bodies2.begin();                                // Set iterator for bodies2
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B, ++B2 ) {// Loop over bodies & bodies2
#ifdef DEBUG
      std::cout << B->ICELL << " " << B->TRG[0] << " " << B2->TRG[0] << std::endl;// Compare every element
#endif
      diff1 += (B->TRG[0] - B2->TRG[0]) * (B->TRG[0] - B2->TRG[0]);// Difference of potential
      norm1 += B2->TRG[0] * B2->TRG[0];                         //  Value of potential
      diff2 += (B->TRG[1] - B2->TRG[1]) * (B->TRG[1] - B2->TRG[1]);// Difference of x acceleration
      diff2 += (B->TRG[2] - B2->TRG[2]) * (B->TRG[2] - B2->TRG[2]);// Difference of y acceleration
      diff2 += (B->TRG[3] - B2->TRG[3]) * (B->TRG[3] - B2->TRG[3]);// Difference of z acceleration
      norm2 += B2->TRG[1] * B2->TRG[1];                         //  Value of x acceleration
      norm2 += B2->TRG[2] * B2->TRG[2];                         //  Value of y acceleration
      norm2 += B2->TRG[3] * B2->TRG[3];                         //  Value of z acceleration
    }                                                           // End loop over bodies & bodies2
  }

//! Print relative L2 norm error
  void printError(real diff1, real norm1, real diff2, real norm2) {
    std::cout << "Error (pot)   : " << std::sqrt(diff1/norm1) << std::endl;
    std::cout << "Error (acc)   : " << std::sqrt(diff2/norm2) << std::endl;
  }
};

#endif
