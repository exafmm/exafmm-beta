#ifndef dataset_h
#define dataset_h
#include "types.h"

namespace {
class Dataset {                                                 // Contains all the different datasets
public:
  void initSource(Bodies &bodies) {                             // Initialize source values
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      B->IBODY = B-bodies.begin();                              //  Tag body with initial index
      B->IPROC = MPIRANK;                                       //  Tag body with initial MPI rank
#if Laplace
      B->Q = 1. / bodies.size() / MPISIZE;                      //  Initialize mass/charge
#elif BiotSavart
      B->Q[0] = (rand() / (1. + RAND_MAX) * 2 * M_PI - M_PI)/ bodies.size() / MPISIZE;// Initialize x vortex strength
      B->Q[1] = (rand() / (1. + RAND_MAX) * 2 * M_PI - M_PI)/ bodies.size() / MPISIZE;// Initialize y vortex strength
      B->Q[2] = (rand() / (1. + RAND_MAX) * 2 * M_PI - M_PI)/ bodies.size() / MPISIZE;// Initialize z vortex strength
      B->S    = 2 * powf(bodies.size(),-1./3);                  //  Initialize core radius
#elif Stretching
      B->Q[0] = (rand() / (1. + RAND_MAX) * 2 * M_PI - M_PI)/ bodies.size() / MPISIZE;// Initialize x vortex strength
      B->Q[1] = (rand() / (1. + RAND_MAX) * 2 * M_PI - M_PI)/ bodies.size() / MPISIZE;// Initialize y vortex strength
      B->Q[2] = (rand() / (1. + RAND_MAX) * 2 * M_PI - M_PI)/ bodies.size() / MPISIZE;// Initialize z vortex strength
      B->S    = 2 * powf(bodies.size(),-1./3);                  //  Initialize core radius
#elif Gaussian
      B->Q = 1. / bodies.size() / MPISIZE;                      //  Initialize mass/charge
      B->S = 2 * powf(bodies.size(),-1./3);                     //  Initialize core radius
#endif
    }                                                           // End loop over bodies
  }

  void initTarget(Bodies &bodies, bool IeqJ=true) {             // Initialize target values
    srand(1);                                                   // Set seed for random number generator
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      B->IBODY = B-bodies.begin();                              //  Tag body with initial index
      B->IPROC = MPIRANK;                                       //  Tag body with initial MPI rank
#if Laplace
      B->pot = -B->Q / std::sqrt(EPS2) * IeqJ;                  //  Initialize potential (0 if I != J)
      B->acc = 0;                                               //  Initialize acceleration
#elif BiotSavart
      B->vel = 0 * IeqJ;                                        //  Initialize velocity
#elif Stretching
      if( !IeqJ ) {                                             //  If source and target are different
        B->Q[0] = (rand() / (1. + RAND_MAX) * 2 * M_PI - M_PI)/ bodies.size() / MPISIZE;// Initialize x vortex strength
        B->Q[1] = (rand() / (1. + RAND_MAX) * 2 * M_PI - M_PI)/ bodies.size() / MPISIZE;// Initialize y vortex strength
        B->Q[2] = (rand() / (1. + RAND_MAX) * 2 * M_PI - M_PI)/ bodies.size() / MPISIZE;// Initialize z vortex strength
      }                                                         //  Endif for different source and target
      B->dQdt = 0;                                              //  Initialize change rate of vortex strength
#elif Gaussian
      B->val = 0 * IeqJ;                                        //  Initialize value
#endif
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
        B->X[d] = rand() / (1. + RAND_MAX) * 2 * M_PI - M_PI;   //   Initialize positions
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
        B->X[d] = rand() / (1. + RAND_MAX) * 2 - 1;             //   Initialize positions
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
    std::ifstream file("data",std::ios::in | std::ios::binary); // Open file
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
#if Laplace
      file >> B->pot;                                           //  Read data for potential
      file >> B->acc[0];                                        //  Read data for x acceleration
      file >> B->acc[1];                                        //  Read data for y acceleration
      file >> B->acc[2];                                        //  Read data for z acceleration
#elif BiotSavart
      file >> B->vel[0];                                        //  Read data for x velocity
      file >> B->vel[1];                                        //  Read data for y velocity
      file >> B->vel[2];                                        //  Read data for z velocity
#elif Stretching
      file >> B->dQdt[0];                                       //  Read data for change rate of x vortex strength
      file >> B->dQdt[1];                                       //  Read data for change rate of y vortex strength
      file >> B->dQdt[2];                                       //  Read data for change rate of z vortex strength
#elif Gaussian
      file >> B->val;                                           //  Read data for value
#endif
    }                                                           // End loop over bodies
    file.close();                                               // Close file
  }

  void writeTarget(Bodies &bodies) {                            // Write target values to file
    std::ofstream file("data",std::ios::out | std::ios::app | std::ios::binary);// Open file
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
#if Laplace
      file << B->pot << std::endl;                              //  Write data for potential
      file << B->acc[0] << std::endl;                           //  Write data for x acceleration
      file << B->acc[1] << std::endl;                           //  Write data for y acceleration
      file << B->acc[2] << std::endl;                           //  Write data for z acceleration
#elif BiotSavart
      file << B->vel[0] << std::endl;                           //  Write data for x velocity
      file << B->vel[1] << std::endl;                           //  Write data for y velocity
      file << B->vel[2] << std::endl;                           //  Write data for z velocity
#elif Stretching
      file << B->dQdt[0] << std::endl;                          //  Write data for change rate of x vortex strength
      file << B->dQdt[1] << std::endl;                          //  Write data for change rate of y vortex strength
      file << B->dQdt[2] << std::endl;                          //  Write data for change rate of z vortex strength
#elif Gaussian
      file << B->val << std::endl;                              //  Write data for value
#endif
    }                                                           // End loop over bodies
    file.close();                                               // Close file

  }

  void evalError(Bodies &bodies, Bodies &bodies2,               // Evaluate error
                 real &diff1, real &norm1, real &diff2, real &norm2) {
#if Laplace
    B_iter B2 = bodies2.begin();                                // Set iterator for bodies2
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B, ++B2 ) {// Loop over bodies & bodies2
#ifdef DEBUG
      std::cout << B->ICELL << " " << B->pot << " " << B2->pot << std::endl;// Compare every element
#endif
      diff1 += (B->pot - B2->pot) * (B->pot - B2->pot);         // Difference of potential
      norm1 += B2->pot * B2->pot;                               // Value of potential
      diff2 += norm(B->acc - B2->acc);                          // Difference of acceleration
      norm2 += norm(B2->acc);                                   // Value of acceleration
    }                                                           // End loop over bodies & bodies2
#elif BiotSavart
    diff2 = norm2 = 0;                                          // Set unused values to 0
    B_iter B2 = bodies2.begin();                                // Set iterator for bodies2
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B, ++B2 ) {// Loop over bodies & bodies2
#ifdef DEBUG
      std::cout << B->ICELL << " " << B->vel[0] << " " << B2->vel[0] << std::endl;// Compare every element
#endif
      diff1 += norm(B->vel - B2->vel);                          // Difference of velocity
      norm1 += norm(B2->vel);                                   // Value of velocity
    }                                                           // End loop over bodies & bodies2
#elif Stretching
    diff2 = norm2 = 0;                                          // Set unused values to 0
    B_iter B2 = bodies2.begin();                                // Set iterator for bodies2
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B, ++B2 ) {// Loop over bodies & bodies2
#ifdef DEBUG
      std::cout << B->ICELL << " " << B->dQdt[0] << " " << B2->dQdt[0] << std::endl;// Compare every element
#endif
      diff1 += norm(B->dQdt - B2->dQdt);                        // Difference of change rate of vortex strength
      norm1 += norm(B2->dQdt);                                  // Value of change rate of vortex strength
    }                                                           // End loop over bodies & bodies2
#elif Gaussian
    diff2 = norm2 = 0;                                          // Set unused values to 0
    B_iter B2 = bodies2.begin();                                // Set iterator for bodies2
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B, ++B2 ) {// Loop over bodies & bodies2
#ifdef DEBUG
      std::cout << B->ICELL << " " << B->val << " " << B2->val << std::endl;// Compare every element
#endif
      diff1 += (B->val - B2->val) * (B->val - B2->val);         // Difference of potential
      norm1 += B2->val * B2->val;                               // Value of potential
    }                                                           // End loop over bodies & bodies2
#endif
  }

  void printError(real diff1, real norm1, real diff2, real norm2) {// Print relative L2 norm error
#if Laplace
  std::cout << "Error (pot)   : " << std::sqrt(diff1/norm1) << std::endl;
  std::cout << "Error (acc)   : " << std::sqrt(diff2/norm2) << std::endl;
#elif BiotSavart
  vect dummy = diff2; dummy = norm2;                            // Use the values so compiler does not complain
  std::cout << "Error         : " << std::sqrt(diff1/norm1) << std::endl;
#elif Stretching
  vect dummy = diff2; dummy = norm2;                            // Use the values so compiler does not complain
  std::cout << "Error         : " << std::sqrt(diff1/norm1) << std::endl;
#elif Gaussian
  vect dummy = diff2; dummy = norm2;                            // Use the values so compiler does not complain
  std::cout << "Error         : " << std::sqrt(diff1/norm1) << std::endl;
#endif
  }
};
}

#endif
