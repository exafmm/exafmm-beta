#ifndef dataset_h
#define dataset_h
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "types.h"

class Dataset {                                                 // Contains all the different datasets
private:
  long filePosition;                                            // Position of file stream

private:
//! Initialize source values
  void initSource(Bodies &bodies, int mpisize) {
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
#if 0
      B->SRC = (drand48() - .5) / bodies.size() / mpisize;      //  Initialize mass/charge
#else
      B->SRC = 1. / bodies.size() / mpisize;                    //  Initialize mass/charge
#endif
    }                                                           // End loop over bodies
  }

//! Uniform distribution on [-1,1]^3 lattice
  void lattice(Bodies &bodies, int mpirank, int mpisize) {
    int level = int(log(bodies.size()*mpisize+1.)/M_LN2/3);     // Level of tree
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      int d = 0, l = 0;                                         //  Initialize dimension and level
      int index = mpirank * bodies.size() + (B-bodies.begin()); //  Set index of body iterator
      vec<3,int> nx = 0;                                        //  Initialize 3-D cell index
      while (index != 0) {                                      //  Deinterleave bits while index is nonzero
        nx[d] += (index & 1) * (1 << l);                        //   Add deinterleaved bit to 3-D cell index
        index >>= 1;                                            //   Right shift the bits
        d = (d+1) % 3;                                          //   Increment dimension
        if( d == 0 ) l++;                                       //   If dimension is 0 again, increment level
      }                                                         //  End while loop for deinterleaving bits
      for (d=0; d<3; d++) {                                     //  Loop over dimensions
        B->X[d] = -1 + (2 * nx[d] + 1.) / (1 << level);         //   Calculate cell center from 3-D cell index
      }                                                         //  End loop over dimensions
    }                                                           // End loop over bodies
  }

//! Random distribution in [-1,1]^3 cube
  void cube(Bodies &bodies, int seed, int numSplit) {
    srand48(seed);                                              // Set seed for random number generator
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      if (numSplit != 1 && B-bodies.begin() == int(seed*bodies.size()/numSplit)) {// Mimic parallel dataset
        seed++;                                                 //   Mimic seed at next rank
        srand48(seed);                                          //   Set seed for random number generator
      }                                                         //  Endif for mimicing parallel dataset
      for (int d=0; d<3; d++) {                                 //  Loop over dimension
        B->X[d] = drand48() * 2 * M_PI - M_PI;                  //   Initialize positions
      }                                                         //  End loop over dimension
    }                                                           // End loop over bodies
  }

//! Random distribution on r = 1 sphere
  void sphere(Bodies &bodies, int seed, int numSplit) {
    srand48(seed);                                              // Set seed for random number generator
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      if (numSplit != 1 && B-bodies.begin() == int(seed*bodies.size()/numSplit)) {// Mimic parallel dataset
        seed++;                                                 //   Mimic seed at next rank
        srand48(seed);                                          //   Set seed for random number generator
      }                                                         //  Endif for mimicing parallel dataset
      for (int d=0; d<3; d++) {                                 //  Loop over dimension
        B->X[d] = drand48() * 2 - 1;                            //   Initialize positions
      }                                                         //  End loop over dimension
      real_t r = std::sqrt(norm(B->X));                         //  Distance from center
      for (int d=0; d<3; d++) {                                 //  Loop over dimension
        B->X[d] /= r * 1.1;                                     //   Normalize positions
      }                                                         //  End loop over dimension
    }                                                           // End loop over bodies
  }

//! Plummer distribution in a r = M_PI/2 sphere
  void plummer(Bodies &bodies, int seed, int numSplit) {
    srand48(seed);                                              // Set seed for random number generator
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      if (numSplit != 1 && B-bodies.begin() == int(seed*bodies.size()/numSplit)) { // Mimic parallel dataset
        seed++;                                                 //   Mimic seed at next rank
        srand48(seed);                                          //   Set seed for random number generator
      }                                                         //  Endif for mimicing parallel dataset
      for (int d=0; d<3; d++) B->X[d] = drand48() * 2 - 1;      //  Generate random number between [-1,1]
      real_t R2 = powf(drand48()*.999, -2.0/3.0) - 1;           //  Generate distribution radius squared
      real_t rscale = 0.015 * M_PI / std::sqrt(norm(B->X) * R2);//  Scaling to fit in M_PI box
      for (int d=0; d<3; d++) B->X[d] *= rscale;                //  Rescale particle coordinates
    }                                                           // End loop over bodies
  }

public:
//! Constructor
  Dataset() : filePosition(0) {}
//! Destructor
  ~Dataset() {}

//! Initialize target values
  void initTarget(Bodies &bodies) {
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      B->TRG = 0;                                               //  Clear target values
    }                                                           // End loop over bodies
  }

//! Initialize dsitribution, source & target value of bodies
  void initBodies(Bodies &bodies, const char * distribution, int mpirank=0, int mpisize=1) {
    switch (distribution[0]) {                                  // Switch between data distribution type
    case 'l':                                                   // Case for lattice
      lattice(bodies,mpirank,mpisize);                          //  Uniform distribution on [-1,1]^3 lattice
      break;                                                    // End case for lattice
    case 'c':                                                   // Case for cube
      cube(bodies,mpirank,mpisize);                             //  Random distribution in [-1,1]^3 cube
      break;                                                    // End case for cube
    case 's':                                                   // Case for sphere
      sphere(bodies,mpirank,mpisize);                           //  Random distribution on surface of r = 1 sphere
      break;                                                    // End case for sphere
    case 'p':                                                   // Case plummer
      plummer(bodies,mpirank,mpisize);                          //  Plummer distribution in a r = M_PI/2 sphere
      break;                                                    // End case for plummer
    default:                                                    // If none of the above
      fprintf(stderr, "unknown data distribution %s\n", distribution);// Print error message
    }                                                           // End switch between data distribution type
    initSource(bodies,mpisize);                                 // Initialize source values
    initTarget(bodies);                                         // Initialize target values
  }

//! Read target values from file
  void readTarget(Bodies &bodies, int mpirank) {
    char fname[256];                                            // File name for saving direct calculation values
    sprintf(fname,"direct%4.4d",mpirank);                       // Set file name
    std::ifstream file(fname,std::ios::in | std::ios::binary);  // Open file
    file.seekg(filePosition);                                   // Set position in file
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      file >> B->TRG[0];                                        //  Read data for potential
      file >> B->TRG[1];                                        //  Read data for x acceleration
      file >> B->TRG[2];                                        //  Read data for y acceleration
      file >> B->TRG[3];                                        //  Read data for z acceleration
    }                                                           // End loop over bodies
    filePosition = file.tellg();                                // Get position in file
    file.close();                                               // Close file
  }

//! Write target values to file
  void writeTarget(Bodies &bodies, int mpirank) {
    char fname[256];                                            // File name for saving direct calculation values
    sprintf(fname,"direct%4.4d",mpirank);                       // Set file name
    std::ofstream file(fname,std::ios::out | std::ios::app | std::ios::binary);// Open file
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      file << B->TRG[0] << std::endl;                           //  Write data for potential
      file << B->TRG[1] << std::endl;                           //  Write data for x acceleration
      file << B->TRG[2] << std::endl;                           //  Write data for y acceleration
      file << B->TRG[3] << std::endl;                           //  Write data for z acceleration
    }                                                           // End loop over bodies
    file.close();                                               // Close file
  }

//! Evaluate relaitve L2 norm error
  void evalError(Bodies &bodies, Bodies &bodies2,
                 double &diff1, double &norm1, double &diff2, double &norm2) {
    B_iter B2 = bodies2.begin();                                // Set iterator for bodies2
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++, B2++) { // Loop over bodies & bodies2
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
  void printError(double diff1, double norm1, double diff2, double norm2) {
    std::cout << std::setw(20) << std::left
              << "Rel. L2 Error (pot)" << " : " << std::sqrt(diff1/norm1) << std::endl;
    std::cout << std::setw(20) << std::left
              << "Rel. L2 Error (acc)" << " : " << std::sqrt(diff2/norm2) << std::endl;
  }
};

#endif
