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
  //! Split range and return partial range
  void splitRange(int & begin, int & end, int iSplit, int numSplit) {
    assert(end > begin);                                        // Check that size > 0
    int size = end - begin;                                     // Size of range
    int increment = size / numSplit;                            // Increment of splitting
    int remainder = size % numSplit;                            // Remainder of splitting
    begin += iSplit * increment + std::min(iSplit,remainder);   // Increment the begin counter
    end = begin + increment;                                    // Increment the end counter
    if (remainder > iSplit) end++;                              // Adjust the end counter for remainder
  }

//! Uniform distribution on [-1,1]^2 lattice
  Bodies lattice(int numBodies, int mpirank, int mpisize) {
    int nx = int(sqrt(numBodies*mpisize));                      // Number of points in x direction
    int ny = nx;                                                // Number of points in y direction
    int begin = 0;                                              // Begin index in y direction
    int end = ny;                                               // End index in the y direction
    splitRange(begin, end, mpirank, mpisize);                   // Split range in y direction
    int numLattice = nx * (end - begin);                        // Total number of lattice points
    Bodies bodies(numLattice);                                  // Initialize bodies
    B_iter B = bodies.begin();                                  // Initialize body iterator
    for (int ix=0; ix<nx; ix++) {                               // Loop over x direction
      for (int iy=begin; iy<end; iy++, B++) {                   //  Loop over y direction
        B->X[0] = (ix / real_t(nx-1)) * 2 - 1;                  //   x coordinate
        B->X[1] = (iy / real_t(nx-1)) * 2 - 1;                  //   y coordinate
      }                                                         //  End loop over y direction
    }                                                           // End loop over x direction 
    return bodies;                                              // Return bodies
  }

//! Random distribution in [-1,1]^2 square
  Bodies square(int numBodies, int seed, int numSplit) {
    srand48(seed);                                              // Set seed for random number generator
    Bodies bodies(numBodies);                                   // Initialize bodies
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      if (numSplit != 1 && B-bodies.begin() == int((seed+1)*bodies.size()/numSplit)) {// Mimic parallel dataset
        seed++;                                                 //   Mimic seed at next rank
        srand48(seed);                                          //   Set seed for random number generator
      }                                                         //  Endif for mimicing parallel dataset
      for (int d=0; d<2; d++) {                                 //  Loop over dimension
        B->X[d] = drand48() * 2 * M_PI - M_PI;                  //   Initialize positions
      }                                                         //  End loop over dimension
    }                                                           // End loop over bodies
    return bodies;                                              // Return bodies
  }

//! Random distribution on r = 1 circle
  Bodies circle(int numBodies, int seed, int numSplit) {
    srand48(seed);                                              // Set seed for random number generator
    Bodies bodies(numBodies);                                   // Initialize bodies
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      if (numSplit != 1 && B-bodies.begin() == int((seed+1)*bodies.size()/numSplit)) {// Mimic parallel dataset
        seed++;                                                 //   Mimic seed at next rank
        srand48(seed);                                          //   Set seed for random number generator
      }                                                         //  Endif for mimicing parallel dataset
      for (int d=0; d<2; d++) {                                 //  Loop over dimension
        B->X[d] = drand48() * 2 - 1;                            //   Initialize positions
      }                                                         //  End loop over dimension
      real_t r = std::sqrt(norm(B->X));                         //  Distance from center
      for (int d=0; d<2; d++) {                                 //  Loop over dimension
        B->X[d] /= r * 1.1;                                     //   Normalize positions
      }                                                         //  End loop over dimension
    }                                                           // End loop over bodies
    return bodies;                                              // Return bodies
  }

 public:
//! Constructor
  Dataset() : filePosition(0) {}

//! Initialize source values
  void initSource(Bodies &bodies) {
#if MASS
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       //  Loop over bodies
      B->SRC = 1. / bodies.size();                              //   Initialize mass
    }                                                           //  End loop over bodies
#else
    real_t average = 0;                                         //  Initialize average charge
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       //  Loop over bodies
      B->SRC = drand48() - .5;                                  //   Initialize charge
      average += B->SRC;                                        //   Accumulate average
    }                                                           //  End loop over bodies
    average /= bodies.size();                                   //  Normalize average
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       //  Loop over bodies
      B->SRC -= average;                                        //   Subtract average charge
    }                                                           //  End loop over bodies
#endif
  }

//! Initialize target values
  void initTarget(Bodies &bodies) {
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      B->TRG = 0;                                               //  Clear target values
    }                                                           // End loop over bodies
  }

//! Initialize dsitribution, source & target value of bodies
  Bodies initBodies(int numBodies, const char * distribution,
		  int mpirank=0, int mpisize=1, int numSplit=1) {
    Bodies bodies;                                              // Initialize bodies
    switch (distribution[0]) {                                  // Switch between data distribution type
    case 'l':                                                   // Case for lattice
      bodies = lattice(numBodies,mpirank,mpisize);              //  Uniform distribution on [-1,1]^2 lattice
      break;                                                    // End case for lattice
    case 's':                                                   // Case for cube
      bodies = square(numBodies,mpirank,numSplit);              //  Random distribution in [-1,1]^2 square
      break;                                                    // End case for cube
    case 'c':                                                   // Case for circle
      bodies = circle(numBodies,mpirank,numSplit);              //  Random distribution on circumference of r = 1 circle
      break;                                                    // End case for circle
    default:                                                    // If none of the above
      fprintf(stderr, "Unknown data distribution %s\n", distribution);// Print error message
    }                                                           // End switch between data distribution type
    initSource(bodies);                                         // Initialize source values
    initTarget(bodies);                                         // Initialize target values
    return bodies;                                              // Return bodies
  }

//! Downsize target bodies by even sampling 
  void sampleBodies(Bodies &bodies, int numTargets) {
    if (numTargets < int(bodies.size())) {                      // If target size is smaller than current
      int stride = bodies.size() / numTargets;                  //  Stride of sampling
      for (int i=0; i<numTargets; i++) {                        //  Loop over target samples
        bodies[i] = bodies[i*stride];                           //   Sample targets
      }                                                         //  End loop over target samples
      bodies.resize(numTargets);                                //  Resize bodies to target size
    }                                                           // End if for target size
  }

//! Get bodies with positive charges
  Bodies getPositive(Bodies &bodies) {
    Bodies buffer = bodies;                                     // Copy bodies to buffer
    B_iter B2 = buffer.begin();                                 // Initialize iterator for buffer
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      if (B->SRC >= 0) {                                        //  If source is positive
        *B2 = *B;                                               //   Copy data to buffer
        B2++;                                                   //   Increment iterator
      }                                                         //  End if for positive source
    }                                                           // End loop over bodies
    buffer.resize(B2-buffer.begin());                           // Resize buffer
    return buffer;                                              // Return buffer
  }


//! Get bodies with negative charges
  Bodies getNegative(Bodies &bodies) {
    Bodies buffer = bodies;                                     // Copy bodies to buffer
    B_iter B2 = buffer.begin();                                 // Initialize iterator for buffer
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      if (B->SRC < 0) {                                         //  If source is negative
        *B2 = *B;                                               //   Copy data to buffer
        B2++;                                                   //   Increment iterator
      }                                                         //  End if for negative source
    }                                                           // End loop over bodies
    buffer.resize(B2-buffer.begin());                           // Resize buffer
    return buffer;                                              // Return buffer
  }

//! Evaluate relaitve L2 norm error
  void evalError(Bodies &bodies, Bodies &bodies2, double &diff1, double &norm1) {
    B_iter B2 = bodies2.begin();                                // Set iterator for bodies2
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++, B2++) { // Loop over bodies & bodies2
      double dp = (B->TRG - B2->TRG) * (B->TRG - B2->TRG);      // Difference of potential
      double  p = B2->TRG * B2->TRG;                             //  Value of potential
      diff1 += dp;                                              //  Accumulate difference of potential
      norm1 += p;                                               //  Accumulate value of potential
    }                                                           // End loop over bodies & bodies2
  }
};
#endif
