#ifndef periodic_h
#define periodic_h
#define EVALUATOR
#include "kernel.h"
#undef EVALUATOR

class Periodic : public Kernel {                                // Periodic boundary condition routines
protected:
  int       Iperiodic;                                          // Periodic image flag
  const int Icenter;                                            // Periodic image flag at center
  Maps      flagM2L;                                            // Flag indicating existance of periodic image for M2L
  Maps      flagM2P;                                            // Flag indicating existance of periodic image for M2P
  Maps      flagP2P;                                            // Flag indicating existance of periodic image for P2P

public:
  Periodic() : Icenter(1 << 13) {}                              // Constructor
  ~Periodic() {}                                                // Destructor

  int getPeriodicRange() {                                      // Get range of periodic images
    int prange = 0;                                             //  Range of periodic images
    for( int i=0; i!=IMAGES; ++i ) {                            //  Loop over periodic image sublevels
      prange += std::pow(3,i);                                  //   Accumulate range of periodic images
    }                                                           //  End loop over perioidc image sublevels
    return prange;                                              // Return range of periodic images
  }

  Bodies periodicBodies(Bodies &bodies) {                       // Create periodic images of bodies
    Bodies jbodies;                                             // Vector for periodic images of bodies
    int prange = getPeriodicRange();                            // Get range of periodic images
    for( int ix=-prange; ix<=prange; ++ix ) {                   // Loop over x periodic direction
      for( int iy=-prange; iy<=prange; ++iy ) {                 //  Loop over y periodic direction
        for( int iz=-prange; iz<=prange; ++iz ) {               //   Loop over z periodic direction
          for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {//    Loop over bodies
            Body body = *B;                                     //     Copy current body
            body.pos[0] += ix * 2 * R0;                         //     Shift x position
            body.pos[1] += iy * 2 * R0;                         //     Shift y position
            body.pos[2] += iz * 2 * R0;                         //     Shift z position
            jbodies.push_back(body);                            //     Push shifted body into jbodies
          }                                                     //    End loop over bodies
        }                                                       //   End loop over z periodic direction
      }                                                         //  End loop over y periodic direction
    }                                                           // End loop over x periodic direction
    return jbodies;                                             // Return vector for periodic images of bodies
  }

};

#endif
