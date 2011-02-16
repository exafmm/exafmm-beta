#ifndef dataset_h
#define dataset_h
#include "types.h"

class Dataset {
public:
  void random(Bodies &bodies) {                                 // Random distribution in [-1,1]^3 cube
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over all bodies
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        B->pos[d] = rand() / (1. + RAND_MAX) * 2 - 1;           //   Initialize positions
      }                                                         //  End loop over dimension
      B->scal = 1. / bodies.size() / MPISIZE;                   //  Initialize source value
      B->acc = B->pot = 0;                                      //  Initialize target value
    }                                                           // End loop over all bodies
  }

  void sphere(Bodies &bodies) {                                 // Random distribution on r = 1 sphere
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over all bodies
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        B->pos[d] = rand() / (1. + RAND_MAX) * 2 - 1;           //   Initialize positions
      }                                                         //  End loop over dimension
      real r = std::sqrt(norm(B->pos));                         //  Distance from center
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        B->pos[d] /= r * 1.1;                                   //   Normalize positions
      }                                                         //  End loop over dimension
      B->scal = 1. / bodies.size() / MPISIZE;                   //  Initialize source value
      B->acc = B->pot = 0;                                      //  Initialize target value
    }                                                           // End loop over all bodies
  }

  void lattice(Bodies &bodies) {                                // Uniform distribution on [-1,1]^3 lattice (for debug)
    int level = int(log(bodies.size()*MPISIZE+1.)/M_LN2/3);     // Level of tree
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over all bodies
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
      B->scal = 1. / bodies.size() / MPISIZE;                   //  Initialize source value
      B->acc = B->pot = 0;                                      //  Initialize target value
    }                                                           // End loop over bodies
  }
};

#endif
