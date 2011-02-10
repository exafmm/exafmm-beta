#ifndef dataset_h
#define dataset_h
#include "types.h"

class Dataset {
private:
  Bodies &bodies;
public:
  Dataset(Bodies &b) : bodies(b) {}                             // Constructor
  ~Dataset() {}                                                 // Destructor

  void random() {
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over all bodies
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        B->pos[d] = rand() / (1. + RAND_MAX) * 2 - 1;           //   Initialize positions
      }                                                         //  End loop over dimension
      B->scal = 1. / bodies.size();                             //  Initialize source value
    }                                                           // End loop over all bodies
  }

  void sphere() {
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over all bodies
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        B->pos[d] = rand() / (1. + RAND_MAX) * 2 - 1;           //   Initialize positions
      }                                                         //  End loop over dimension
      real r = std::sqrt(norm(B->pos));                         //  Distance from center
      for( int d=0; d!=3; ++d ) {                               //  Loop over dimension
        B->pos[d] /= r * 1.1;                                   //   Normalize positions
      }                                                         //  End loop over dimension
      B->scal = 1. / bodies.size();                             //  Initialize source value
    }                                                           // End loop over all bodies
  }

};

#endif
