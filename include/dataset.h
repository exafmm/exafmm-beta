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
#include "types.h"

//! Contains all the different datasets
template<Equation equation>
class Dataset {
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
      B->SRC = 0;                                               //  Clear previous source values
      B->SRC[0] = 1. / bodies.size() / MPISIZE;                 //  Initialize mass/charge
    }                                                           // End loop over bodies
  }

//! Initialize target values
  void initTarget(Bodies &bodies, bool IeqJ=true) {
    srand48(0);                                                 // Set seed for random number generator
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      B->IBODY = B-bodies.begin();                              //  Tag body with initial index
      B->IPROC = MPIRANK;                                       //  Tag body with initial MPI rank
      B->TRG = 0;                                               //  Clear previous target values (IeqJ is dummy)
      B->TRG[0] = -B->SRC[0] / std::sqrt(EPS2) * IeqJ;          //  Initialize potential (0 if I != J)
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
    B_iter B2 = bodies2.begin();                                //  Set iterator for bodies2
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
    }                                                           //  End loop over bodies & bodies2
  }

//! Print relative L2 norm error
  void printError(real diff1, real norm1, real diff2, real norm2) {
    std::cout << "Error (pot)   : " << std::sqrt(diff1/norm1) << std::endl;
    std::cout << "Error (acc)   : " << std::sqrt(diff2/norm2) << std::endl;
  }
};

template<>
class Dataset<BiotSavart> {
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
      B->SRC[0] = (drand48() * 2 - 1)/ bodies.size() / MPISIZE; //  Initialize x vortex strength
      B->SRC[1] = (drand48() * 2 - 1)/ bodies.size() / MPISIZE; //  Initialize y vortex strength
      B->SRC[2] = (drand48() * 2 - 1)/ bodies.size() / MPISIZE; //  Initialize z vortex strength
      B->SRC[3] = powf(bodies.size() * MPISIZE,-1./3);          //  Initialize core radius
    }                                                           // End loop over bodies
  }

//! Initialize target values 
  void initTarget(Bodies &bodies, bool) {
    srand48(0);                                                 // Set seed for random number generator
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      B->IBODY = B-bodies.begin();                              //  Tag body with initial index
      B->IPROC = MPIRANK;                                       //  Tag body with initial MPI rank
      B->TRG = 0;                                               //  Clear previous target values (IeqJ is dummy)
    }                                                           // End loop over bodies
  }

//! Read target values from file
  void readTarget(Bodies &bodies) {
    char fname[256];                                            // File name for saving direct calculation values
    sprintf(fname,"direct%4.4d",MPIRANK);                       // Set file name
    std::ifstream file(fname,std::ios::in | std::ios::binary);  // Open file
    file.seekg(filePosition);                                   // Set position in file
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      file >> B->TRG[0];                                        //  Read data for x velocity
      file >> B->TRG[1];                                        //  Read data for y velocity
      file >> B->TRG[2];                                        //  Read data for z velocity
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
      file << B->TRG[0] << std::endl;                           //  Write data for x velocity
      file << B->TRG[1] << std::endl;                           //  Write data for y velocity
      file << B->TRG[2] << std::endl;                           //  Write data for z velocity
    }                                                           // End loop over bodies
    file.close();                                               // Close file
  }

//! Evaluate relative L2 norm error
  void evalError(Bodies &bodies, Bodies &bodies2,
                 real &diff1, real &norm1, real &diff2, real &norm2) {
    diff2 = norm2 = 0;                                          // Set unused values to 0
    B_iter B2 = bodies2.begin();                                // Set iterator for bodies2
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B, ++B2 ) {// Loop over bodies & bodies2
#ifdef DEBUG
      std::cout << B->ICELL << " " << B->TRG[0] << " " << B2->TRG[0] << std::endl;// Compare every element
#endif
      diff1 += (B->TRG[0] - B2->TRG[0]) * (B->TRG[0] - B2->TRG[0]);// Difference of x velocity
      diff1 += (B->TRG[1] - B2->TRG[1]) * (B->TRG[1] - B2->TRG[1]);// Difference of y velocity
      diff1 += (B->TRG[2] - B2->TRG[2]) * (B->TRG[2] - B2->TRG[2]);// Difference of z velocity
      norm1 += B2->TRG[0] * B2->TRG[0];                         //  Value of x velocity
      norm1 += B2->TRG[1] * B2->TRG[1];                         //  Value of y velocity
      norm1 += B2->TRG[2] * B2->TRG[2];                         //  Value of z velocity
    }                                                           // End loop over bodies & bodies2
  }

//! Print relative L2 norm error
  void printError(real diff1, real norm1, real diff2, real norm2) {
    vect dummy = diff2; dummy = norm2;                          //  Use the values so compiler does not complain
    std::cout << "Error         : " << std::sqrt(diff1/norm1) << std::endl;
  }
};

template<>
class Dataset<Stretching> {
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
      B->SRC[0] = (drand48() * 2 - 1) / bodies.size() / MPISIZE;//  Initialize x vortex strength
      B->SRC[1] = (drand48() * 2 - 1) / bodies.size() / MPISIZE;//  Initialize y vortex strength
      B->SRC[2] = (drand48() * 2 - 1) / bodies.size() / MPISIZE;//  Initialize z vortex strength
      B->SRC[3] = powf(bodies.size() * MPISIZE,-1./3);          //  Initialize core radius
    }                                                           // End loop over bodies
  }

//! Initialize target values
  void initTarget(Bodies &bodies, bool IeqJ=true) {
    srand48(0);                                                 // Set seed for random number generator
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      B->IBODY = B-bodies.begin();                              //  Tag body with initial index
      B->IPROC = MPIRANK;                                       //  Tag body with initial MPI rank
      B->TRG = 0;                                               //  Clear previous target values (IeqJ is dummy)
      if( !IeqJ ) {                                             //  If source and target are different
        B->SRC[0] = (drand48() * 2 - 1) / bodies.size() / MPISIZE;// Initialize x vortex strength
        B->SRC[1] = (drand48() * 2 - 1) / bodies.size() / MPISIZE;// Initialize y vortex strength
        B->SRC[2] = (drand48() * 2 - 1) / bodies.size() / MPISIZE;// Initialize z vortex strength
      }                                                         //  Endif for different source and target
    }                                                           // End loop over bodies
  }

//! Read target values from file
  void readTarget(Bodies &bodies) {
    char fname[256];                                            // File name for saving direct calculation values
    sprintf(fname,"direct%4.4d",MPIRANK);                       // Set file name
    std::ifstream file(fname,std::ios::in | std::ios::binary);  // Open file
    file.seekg(filePosition);                                   // Set position in file
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      file >> B->TRG[0];                                        //  Read data for change rate of x vortex strength
      file >> B->TRG[1];                                        //  Read data for change rate of y vortex strength
      file >> B->TRG[2];                                        //  Read data for change rate of z vortex strength
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
      file << B->TRG[0] << std::endl;                           //  Write data for change rate of x vortex strength
      file << B->TRG[1] << std::endl;                           //  Write data for change rate of y vortex strength
      file << B->TRG[2] << std::endl;                           //  Write data for change rate of z vortex strength
    }                                                           // End loop over bodies
    file.close();                                               // Close file
  }

//! Evaluate relative L2 norm error
  void evalError(Bodies &bodies, Bodies &bodies2,
                 real &diff1, real &norm1, real &diff2, real &norm2) {
    diff2 = norm2 = 0;                                          // Set unused values to 0
    B_iter B2 = bodies2.begin();                                // Set iterator for bodies2
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B, ++B2 ) {// Loop over bodies & bodies2
#ifdef DEBUG
      std::cout << B->ICELL << " " << B->TRG[0] << " " << B2->TRG[0] << std::endl;// Compare every element
#endif
      diff1 += (B->TRG[0] - B2->TRG[0]) * (B->TRG[0] - B2->TRG[0]);// Difference of x change rate of vortex strength
      diff1 += (B->TRG[1] - B2->TRG[1]) * (B->TRG[1] - B2->TRG[1]);// Difference of y change rate of vortex strength
      diff1 += (B->TRG[2] - B2->TRG[2]) * (B->TRG[2] - B2->TRG[2]);// Difference of z change rate of vortex strength
      norm1 += B2->TRG[0] * B2->TRG[0];                         //  Value of x change rate of vortex strength
      norm1 += B2->TRG[1] * B2->TRG[1];                         //  Value of y change rate of vortex strength
      norm1 += B2->TRG[2] * B2->TRG[2];                         //  Value of z change rate of vortex strength
    }                                                           // End loop over bodies & bodies2
  }

//! Print relative L2 norm error
  void printError(real diff1, real norm1, real diff2, real norm2) {
    vect dummy = diff2; dummy = norm2;                          // Use the values so compiler does not complain
    std::cout << "Error         : " << std::sqrt(diff1/norm1) << std::endl;
  }
};

template<>
class Dataset<Gaussian> {
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
      B->SRC = 0;                                               //  Clear previous source values
      B->SRC[0] = 1. / bodies.size() / MPISIZE;                 //  Initialize mass/charge
      B->SRC[3] = powf(bodies.size() * MPISIZE,-1./3);          //  Initialize core radius
    }                                                           // End loop over bodies
  }

//! Initialize target values
  void initTarget(Bodies &bodies, bool) {
    srand48(0);                                                 // Set seed for random number generator
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      B->IBODY = B-bodies.begin();                              //  Tag body with initial index
      B->IPROC = MPIRANK;                                       //  Tag body with initial MPI rank
      B->TRG = 0;                                               //  Clear previous target values (IeqJ is dummy)
    }                                                           // End loop over bodies
  }

//! Read target values from file
  void readTarget(Bodies &bodies) {
    char fname[256];                                            // File name for saving direct calculation values
    sprintf(fname,"direct%4.4d",MPIRANK);                       // Set file name
    std::ifstream file(fname,std::ios::in | std::ios::binary);  // Open file
    file.seekg(filePosition);                                   // Set position in file
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      file >> B->TRG[0];                                        //  Read data for value
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
      file << B->TRG[0] << std::endl;                           //  Write data for value
    }                                                           // End loop over bodies
    file.close();                                               // Close file
  }

//! Evaluate relative L2 norm error
  void evalError(Bodies &bodies, Bodies &bodies2,
                 real &diff1, real &norm1, real &diff2, real &norm2) {
    diff2 = norm2 = 0;                                          // Set unused values to 0
    B_iter B2 = bodies2.begin();                                // Set iterator for bodies2
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B, ++B2 ) {// Loop over bodies & bodies2
#ifdef DEBUG
      std::cout << B->ICELL << " " << B->TRG[0] << " " << B2->TRG[0] << std::endl;// Compare every element
#endif
      diff1 += (B->TRG[0] - B2->TRG[0]) * (B->TRG[0] - B2->TRG[0]);// Difference of potential
      norm1 += B2->TRG[0] * B2->TRG[0];                         //  Value of potential
    }                                                           // End loop over bodies & bodies2
  }

//! Print relative L2 norm error
  void printError(real diff1, real norm1, real diff2, real norm2) {
    vect dummy = diff2; dummy = norm2;                          // Use the values so compiler does not complain
    std::cout << "Error         : " << std::sqrt(diff1/norm1) << std::endl;
  }
};

template<>
class Dataset<CoulombVdW> {
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
      B->SRC = 0;                                               //  Clear previous source values
      B->SRC[0] = 1. / bodies.size() / MPISIZE;                 //  Initialize mass/charge
    }                                                           // End loop over bodies
  }

//! Initialize target values
  void initTarget(Bodies &bodies, bool IeqJ=true) {
    srand48(0);                                                 // Set seed for random number generator
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      B->IBODY = B-bodies.begin();                              //  Tag body with initial index
      B->IPROC = MPIRANK;                                       //  Tag body with initial MPI rank
      B->TRG = 0;                                               //  Clear previous target values (IeqJ is dummy)
      B->TRG[0] = -B->SRC[0] / std::sqrt(EPS2) * IeqJ;          //  Initialize potential (0 if I != J)
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
