#ifndef dataset_h
#define dataset_h
#include "types.h"

class Dataset {
public:
  void random(Bodies &bodies) {
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over all bodies
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        B->pos[d] = rand() / (1. + RAND_MAX) * 2 - 1;           //   Initialize positions
      }                                                         //  End loop over dimension
      B->scal = 1. / bodies.size();                             //  Initialize source value
      B->acc = B->pot = 0;                                      //  Initialize target value
    }                                                           // End loop over all bodies
  }

  void sphere(Bodies &bodies) {
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over all bodies
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        B->pos[d] = rand() / (1. + RAND_MAX) * 2 - 1;           //   Initialize positions
      }                                                         //  End loop over dimension
      real r = std::sqrt(norm(B->pos));                         //  Distance from center
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        B->pos[d] /= r * 1.1;                                   //   Normalize positions
      }                                                         //  End loop over dimension
      B->scal = 1. / bodies.size();                             //  Initialize source value
      B->acc = B->pot = 0;                                      //  Initialize target value
    }                                                           // End loop over all bodies
  }

};

#endif
