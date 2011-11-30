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

class Dataset {                                                 // Contains all the different datasets
private:
  long filePosition;                                            // Position of file stream

public:
  std::string kernelName;                                       // Name of kernel

  Dataset() : filePosition(0) {}                                // Constructor
  ~Dataset() {}                                                 // Destructor

  void initSource(Bodies &bodies) {                             // Initialize source values
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      B->IBODY = B-bodies.begin();                              //  Tag body with initial index
      B->IPROC = MPIRANK;                                       //  Tag body with initial MPI rank
      B->SRC = 0;                                               //  Clear previous source values
      if( kernelName == "Laplace" ) {                           //  If Laplace kernel
        B->SRC[0] = 1. / bodies.size() / MPISIZE;               //   Initialize mass/charge
      } else if ( kernelName == "BiotSavart" ) {                //  If Biot Savart kernel
        B->SRC[0] = (drand48() * 2 - 1) / bodies.size() / MPISIZE;// Initialize x vortex strength
        B->SRC[1] = (drand48() * 2 - 1) / bodies.size() / MPISIZE;// Initialize y vortex strength
        B->SRC[2] = (drand48() * 2 - 1) / bodies.size() / MPISIZE;// Initialize z vortex strength
        B->SRC[3] = powf(bodies.size() * MPISIZE,-1./3);        //   Initialize core radius
      } else if ( kernelName == "Stretching" ) {                //  If Stretching kernel
        B->SRC[0] = (drand48() * 2 - 1) / bodies.size() / MPISIZE;// Initialize x vortex strength
        B->SRC[1] = (drand48() * 2 - 1) / bodies.size() / MPISIZE;// Initialize y vortex strength
        B->SRC[2] = (drand48() * 2 - 1) / bodies.size() / MPISIZE;// Initialize z vortex strength
        B->SRC[3] = powf(bodies.size() * MPISIZE,-1./3);        //   Initialize core radius
      } else if ( kernelName == "Gaussian" ) {                  //  If Gaussian kernel
        B->SRC[0] = 1. / bodies.size() / MPISIZE;               //   Initialize mass/charge
        B->SRC[3] = powf(bodies.size() * MPISIZE,-1./3);        //   Initialize core radius
      } else {                                                  //  If kernel is none of the above
        if(MPIRANK == 0) std::cout << "Invalid kernel type" << std::endl;// Invalid kernel type
        abort();                                                //   Abort execution
      }                                                         //  Endif for kernel type
    }                                                           // End loop over bodies
  }

  void initTarget(Bodies &bodies, bool IeqJ=true) {             // Initialize target values
    srand(1);                                                   // Set seed for random number generator
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
//      B->IBODY = B-bodies.begin();                              //  Tag body with initial index
      B->IPROC = MPIRANK;                                       //  Tag body with initial MPI rank
      B->TRG = 0 * IeqJ;                                        //  Clear previous target values (IeqJ is dummy)
      if( kernelName == "Laplace" ) {                           //  If Laplace kernel
        B->TRG = 0;                                             //   Initialize potential (0 if I != J)
      } else if ( kernelName == "BiotSavart" ) {                //  If Biot Savart kernel
      } else if ( kernelName == "Stretching" ) {                //  If Stretching kernel
        if( !IeqJ ) {                                           //   If source and target are different
          B->SRC[0] = (drand48() * 2 - 1) / bodies.size() / MPISIZE;// Initialize x vortex strength
          B->SRC[1] = (drand48() * 2 - 1) / bodies.size() / MPISIZE;// Initialize y vortex strength
          B->SRC[2] = (drand48() * 2 - 1) / bodies.size() / MPISIZE;// Initialize z vortex strength
        }                                                       //   Endif for different source and target
      } else if ( kernelName == "Gaussian" ) {                  //  If Gaussian kernel
      } else {                                                  //  If kernel is none of the above
        if(MPIRANK == 0) std::cout << "Invalid kernel type" << std::endl;// Invalid kernel type
        abort();                                                //   Abort execution
      }                                                         //  Endif for kernel type
    }                                                           // End loop over bodies
  }

  void random(Bodies &bodies, int seed=1, int numSplit=1) {     // Random distribution in [-1,1]^3 cube
    srand(seed);                                                // Set seed for random number generator
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      if( numSplit != 1 && B-bodies.begin() == int(seed*bodies.size()/numSplit) ) {// Mimic parallel dataset
        seed++;                                                 //   Mimic seed at next rank
        srand(seed);                                            //   Set seed for random number generator
      }                                                         //  Endif for mimicing parallel dataset
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        B->X[d] = drand48() * 2 * M_PI - M_PI;                  //   Initialize positions
      }                                                         //  End loop over dimension
    }                                                           // End loop over bodies
    initSource(bodies);                                         // Initialize source values
    initTarget(bodies);                                         // Initialize target values
  }

  void sphere(Bodies &bodies, int seed=1, int numSplit=1) {     // Random distribution on r = 1 sphere
    srand(seed);                                                // Set seed for random number generator
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      if( numSplit != 1 && B-bodies.begin() == int(seed*bodies.size()/numSplit) ) {// Mimic parallel dataset
        seed++;                                                 //   Mimic seed at next rank
        srand(seed);                                            //   Set seed for random number generator
      }                                                         //  Endif for mimicing parallel dataset
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        B->X[d] = drand48() * 2 - 1;                            //   Initialize positions
      }                                                         //  End loop over dimension
      real r = std::sqrt(norm(B->X));                           //  Distance from center
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        B->X[d] /= r * 1.1;                                     //   Normalize positions
      }                                                         //  End loop over dimension
    }                                                           // End loop over bodies
    initSource(bodies);                                         // Initialize source values
    initTarget(bodies);                                         // Initialize target values
  }

  void lattice(Bodies &bodies) {                                // Uniform distribution on [-1,1]^3 lattice (for debug)
    int level = int(log(bodies.size()*MPISIZE+1.)/M_LN2/3);     // Level of tree
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      int d = 0, l = 0;                                         //  Initialize dimension and level
      int index = MPIRANK * bodies.size() + (B-bodies.begin()); //  Set index of body iterator
      vec<3,int> nx = 0;                                        //  Initialize 3-D cell index
      while( index != 0 ) {                                     //  Deinterleave bits while index is nonzero
        nx[d] += (index % 2) * (1 << l);                        //   Add deinterleaved bit to 3-D cell index
        index >>= 1;                                            //   Right shift the bits
        d = (d+1) % 3;                                          //   Increment dimension
        if( d == 0 ) l++;                                       //   If dimension is 0 again, increment level
      }                                                         //  End while loop for deinterleaving bits
      for( d=0; d!=3; ++d ) {                                   //  Loop over dimensions
        B->X[d] = -1 + (2 * nx[d] + 1.) / (1 << level);         //   Calculate cell center from 3-D cell index
      }                                                         //  End loop over dimensions
    }                                                           // End loop over bodies
    initSource(bodies);                                         // Initialize source values
    initTarget(bodies);                                         // Initialize target values
  }

  void readTarget(Bodies &bodies) {                             // Read target values from file
    char fname[256];                                            // File name for saving direct calculation values
    sprintf(fname,"direct%4.4d",MPIRANK);                       // Set file name
    std::ifstream file(fname,std::ios::in | std::ios::binary);  // Open file
    file.seekg(filePosition);                                   // Set position in file
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      if( kernelName == "Laplace" ) {                           //  If Laplace kernel
        file >> B->TRG[0];                                      //   Read data for potential
        file >> B->TRG[1];                                      //   Read data for x acceleration
        file >> B->TRG[2];                                      //   Read data for y acceleration
        file >> B->TRG[3];                                      //   Read data for z acceleration
      } else if ( kernelName == "BiotSavart" ) {                //  If Biot Savart kernel
        file >> B->TRG[0];                                      //   Read data for x velocity
        file >> B->TRG[1];                                      //   Read data for y velocity
        file >> B->TRG[2];                                      //   Read data for z velocity
      } else if ( kernelName == "Stretching" ) {                //  If Stretching kernel
        file >> B->TRG[0];                                      //   Read data for change rate of x vortex strength
        file >> B->TRG[1];                                      //   Read data for change rate of y vortex strength
        file >> B->TRG[2];                                      //   Read data for change rate of z vortex strength
      } else if ( kernelName == "Gaussian" ) {                  //  If Gaussian kernel
        file >> B->TRG[0];                                      //   Read data for value
      } else {                                                  //  If kernel is none of the above
        if(MPIRANK == 0) std::cout << "Invalid kernel type" << std::endl;// Invalid kernel type
        abort();                                                //   Abort execution
      }                                                         //  Endif for kernel type
    }                                                           // End loop over bodies
    filePosition = file.tellg();                                // Get position in file
    file.close();                                               // Close file
  }

  void writeTarget(Bodies &bodies) {                            // Write target values to file
    char fname[256];                                            // File name for saving direct calculation values
    sprintf(fname,"direct%4.4d",MPIRANK);                       // Set file name
    std::ofstream file(fname,std::ios::out | std::ios::app | std::ios::binary);// Open file
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      if( kernelName == "Laplace" ) {                           //  If Laplace kernel
        file << B->TRG[0] << std::endl;                         //   Write data for potential
        file << B->TRG[1] << std::endl;                         //   Write data for x acceleration
        file << B->TRG[2] << std::endl;                         //   Write data for y acceleration
        file << B->TRG[3] << std::endl;                         //   Write data for z acceleration
      } else if ( kernelName == "BiotSavart" ) {                //  If Biot Savart kernel
        file << B->TRG[0] << std::endl;                         //   Write data for x velocity
        file << B->TRG[1] << std::endl;                         //   Write data for y velocity
        file << B->TRG[2] << std::endl;                         //   Write data for z velocity
      } else if ( kernelName == "Stretching" ) {                //  If Stretching kernel
        file << B->TRG[0] << std::endl;                         //   Write data for change rate of x vortex strength
        file << B->TRG[1] << std::endl;                         //   Write data for change rate of y vortex strength
        file << B->TRG[2] << std::endl;                         //   Write data for change rate of z vortex strength
      } else if ( kernelName == "Gaussian" ) {                  //  If Gaussian kernel
        file << B->TRG[0] << std::endl;                         //   Write data for value
      } else {                                                  //  If kernel is none of the above
        if(MPIRANK == 0) std::cout << "Invalid kernel type" << std::endl;// Invalid kernel type
        abort();                                                //   Abort execution
      }                                                         //  Endif for kernel type
    }                                                           // End loop over bodies
    file.close();                                               // Close file

  }

  void evalError(Bodies &bodies, Bodies &bodies2,               // Evaluate error
                 real &diff1, real &norm1, real &diff2, real &norm2) {
    if( kernelName == "Laplace" ) {                             // If Laplace kernel
      B_iter B2 = bodies2.begin();                              //  Set iterator for bodies2
      for( B_iter B=bodies.begin(); B!=bodies.end(); ++B, ++B2 ) {// Loop over bodies & bodies2
#ifdef DEBUG
        std::cout << B->IBODY << " " << B->TRG[0] << " " << B2->TRG[0] << std::endl;// Compare every element
#endif
        diff1 += (B->TRG[0] - B2->TRG[0]) * (B->TRG[0] - B2->TRG[0]);// Difference of potential
        norm1 += B2->TRG[0] * B2->TRG[0];                       //  Value of potential
        diff2 += (B->TRG[1] - B2->TRG[1]) * (B->TRG[1] - B2->TRG[1]);// Difference of x acceleration
        diff2 += (B->TRG[2] - B2->TRG[2]) * (B->TRG[2] - B2->TRG[2]);// Difference of y acceleration
        diff2 += (B->TRG[3] - B2->TRG[3]) * (B->TRG[3] - B2->TRG[3]);// Difference of z acceleration
        norm2 += B2->TRG[1] * B2->TRG[1];                       //  Value of x acceleration
        norm2 += B2->TRG[2] * B2->TRG[2];                       //  Value of y acceleration
        norm2 += B2->TRG[3] * B2->TRG[3];                       //  Value of z acceleration
      }                                                         //  End loop over bodies & bodies2
    } else if ( kernelName == "BiotSavart" ) {                  // If Biot Savart kernel
      diff2 = norm2 = 0;                                        //  Set unused values to 0
      B_iter B2 = bodies2.begin();                              //  Set iterator for bodies2
      for( B_iter B=bodies.begin(); B!=bodies.end(); ++B, ++B2 ) {// Loop over bodies & bodies2
#ifdef DEBUG
        std::cout << B->ICELL << " " << B->TRG[0] << " " << B2->TRG[0] << std::endl;// Compare every element
#endif
        diff1 += (B->TRG[0] - B2->TRG[0]) * (B->TRG[0] - B2->TRG[0]);// Difference of x velocity
        diff1 += (B->TRG[1] - B2->TRG[1]) * (B->TRG[1] - B2->TRG[1]);// Difference of y velocity
        diff1 += (B->TRG[2] - B2->TRG[2]) * (B->TRG[2] - B2->TRG[2]);// Difference of z velocity
        norm1 += B2->TRG[0] * B2->TRG[0];                       //  Value of x velocity
        norm1 += B2->TRG[1] * B2->TRG[1];                       //  Value of y velocity
        norm1 += B2->TRG[2] * B2->TRG[2];                       //  Value of z velocity
      }                                                         //  End loop over bodies & bodies2
    } else if ( kernelName == "Stretching" ) {                  // If Stretching kernel
      diff2 = norm2 = 0;                                        //  Set unused values to 0
      B_iter B2 = bodies2.begin();                              //  Set iterator for bodies2
      for( B_iter B=bodies.begin(); B!=bodies.end(); ++B, ++B2 ) {// Loop over bodies & bodies2
#ifdef DEBUG
        std::cout << B->ICELL << " " << B->TRG[0] << " " << B2->TRG[0] << std::endl;// Compare every element
#endif
        diff1 += (B->TRG[0] - B2->TRG[0]) * (B->TRG[0] - B2->TRG[0]);// Difference of x change rate of vortex strength
        diff1 += (B->TRG[1] - B2->TRG[1]) * (B->TRG[1] - B2->TRG[1]);// Difference of y change rate of vortex strength
        diff1 += (B->TRG[2] - B2->TRG[2]) * (B->TRG[2] - B2->TRG[2]);// Difference of z change rate of vortex strength
        norm1 += B2->TRG[0] * B2->TRG[0];                       //  Value of x change rate of vortex strength
        norm1 += B2->TRG[1] * B2->TRG[1];                       //  Value of y change rate of vortex strength
        norm1 += B2->TRG[2] * B2->TRG[2];                       //  Value of z change rate of vortex strength
      }                                                         // End loop over bodies & bodies2
    } else if ( kernelName == "Gaussian" ) {                    // If Gaussian kernel
      diff2 = norm2 = 0;                                        //  Set unused values to 0
      B_iter B2 = bodies2.begin();                              //  Set iterator for bodies2
      for( B_iter B=bodies.begin(); B!=bodies.end(); ++B, ++B2 ) {// Loop over bodies & bodies2
#ifdef DEBUG
        std::cout << B->ICELL << " " << B->TRG[0] << " " << B2->TRG[0] << std::endl;// Compare every element
#endif
        diff1 += (B->TRG[0] - B2->TRG[0]) * (B->TRG[0] - B2->TRG[0]);// Difference of potential
        norm1 += B2->TRG[0] * B2->TRG[0];                       //  Value of potential
      }                                                         //  End loop over bodies & bodies2
    } else {                                                    // If kernel is none of the above
      if(MPIRANK == 0) std::cout << "Invalid kernel type" << std::endl;// Invalid kernel type
      abort();                                                  //  Abort execution
    }                                                           // Endif for kernel type
  }

  void printError(real diff1, real norm1, real diff2, real norm2) {// Print relative L2 norm error
    if( kernelName == "Laplace" ) {                             // If Laplace kernel
      std::cout << "Error (pot)   : " << std::sqrt(diff1/norm1) << std::endl;
      std::cout << "Error (acc)   : " << std::sqrt(diff2/norm2) << std::endl;
    } else if ( kernelName == "BiotSavart" ) {                  // If Biot Savart kernel
      vect dummy = diff2; dummy = norm2;                        //  Use the values so compiler does not complain
      std::cout << "Error         : " << std::sqrt(diff1/norm1) << std::endl;
    } else if ( kernelName == "Stretching" ) {                  // If Stretching kernel
      vect dummy = diff2; dummy = norm2;                        //  Use the values so compiler does not complain
      std::cout << "Error         : " << std::sqrt(diff1/norm1) << std::endl;
    } else if ( kernelName == "Gaussian" ) {                    // If Gaussian kernel
      vect dummy = diff2; dummy = norm2;                        //  Use the values so compiler does not complain
      std::cout << "Error         : " << std::sqrt(diff1/norm1) << std::endl;
    } else {                                                    // If kernel is none of the above
      if(MPIRANK == 0) std::cout << "Invalid kernel type" << std::endl;// Invalid kernel type
      abort();                                                  //  Abort execution
    }                                                           // Endif for kernel type
  }
};

#endif
