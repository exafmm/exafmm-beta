#ifndef dataset_h
#define dataset_h
#include "types.h"

class Dataset {                                                 // Contains all the different datasets
public:
  void random(Bodies &bodies, int seed=1, int numSplit=1) {     // Random distribution in [-1,1]^3 cube
    srand(seed);                                                // Set seed for random number generator
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      if( numSplit != 1 && B-bodies.begin() == int(seed*bodies.size()/numSplit) ) {// Mimic parallel dataset
        seed++;                                                 //   Mimic seed at next rank
        srand(seed);                                            //   Set seed for random number generator
      }                                                         //  Endif for mimicing parallel dataset
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        B->pos[d] = rand() / (1. + RAND_MAX) * 2 - 1;           //   Initialize positions
      }                                                         //  End loop over dimension
      B->scal = 1. / bodies.size() / MPISIZE;                   //  Initialize mass/charge
      B->acc = B->pot = 0;                                      //  Initialize target values
    }                                                           // End loop over bodies
  }

  void sphere(Bodies &bodies, int seed=1, int numSplit=1) {     // Random distribution on r = 1 sphere
    srand(seed);                                                // Set seed for random number generator
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      if( numSplit != 1 && B-bodies.begin() == int(seed*bodies.size()/numSplit) ) {// Mimic parallel dataset
        seed++;                                                 //   Mimic seed at next rank
        srand(seed);                                            //   Set seed for random number generator
      }                                                         //  Endif for mimicing parallel dataset
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        B->pos[d] = rand() / (1. + RAND_MAX) * 2 - 1;           //   Initialize positions
      }                                                         //  End loop over dimension
      real r = std::sqrt(norm(B->pos));                         //  Distance from center
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        B->pos[d] /= r * 1.1;                                   //   Normalize positions
      }                                                         //  End loop over dimension
      B->scal = 1. / bodies.size() / MPISIZE;                   //  Initialize mass/charge
      B->acc = B->pot = 0;                                      //  Initialize target values
    }                                                           // End loop over bodies
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
        B->pos[d] = -1 + (2 * nx[d] + 1.) / (1 << level);       //   Calculate cell center from 3-D cell index
      }                                                         //  End loop over dimensions
      B->scal = 1. / bodies.size() / MPISIZE;                   //  Initialize mass/charge
      B->acc = B->pot = 0;                                      //  Initialize target values
    }                                                           // End loop over bodies
  }

};

#endif
