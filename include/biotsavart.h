#include "dataset.h"

void Dataset::initSource(Bodies &bodies) {                      // Initialize source values
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {        // Loop over bodies
    B->Q[0] = (rand() / (1. + RAND_MAX) * 2 - 1)/ bodies.size() / MPISIZE;//  Initialize x vortex strength
    B->Q[1] = (rand() / (1. + RAND_MAX) * 2 - 1)/ bodies.size() / MPISIZE;//  Initialize y vortex strength
    B->Q[2] = (rand() / (1. + RAND_MAX) * 2 - 1)/ bodies.size() / MPISIZE;//  Initialize z vortex strength
    B->S    = 2 * powf(bodies.size(),-1.0/3);                   // Initialize core radius
  }                                                             // End loop over bodies
}

void Dataset::initTarget(Bodies &bodies, bool IeqJ) {           // Initialize target values
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {        // Loop over bodies
    B->vel = 0 * IeqJ;                                          //  Initialize velocity
  }                                                             // End loop over bodies
}

void Dataset::readTarget(Bodies &bodies) {                      // Read target data from file
  std::ifstream file("data",std::ios::in | std::ios::ate | std::ios::binary);// Open file
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {        // Loop over bodies
    file >> B->vel[0];                                          //  Read data for x velocity
    file >> B->vel[1];                                          //  Read data for y velocity
    file >> B->vel[2];                                          //  Read data for z velocity
  }                                                             // End loop over bodies
  file.close();                                                 // Close file
}

void Dataset::writeTarget(Bodies &bodies) {                     // Write target data to file
  std::ofstream file("data",std::ios::out | std::ios::app | std::ios::binary);// Open file
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {        // Loop over bodies
    file << B->vel[0] << std::endl;                             //  Write data for x velocity
    file << B->vel[1] << std::endl;                             //  Write data for y velocity
    file << B->vel[2] << std::endl;                             //  Write data for z velocity
  }                                                             // End loop over bodies
  file.close();                                                 // Close file
}

void Dataset::evalError(Bodies &bodies, Bodies &bodies2,        // Evaluate error
                        real &velDiff, real &velNorm, real &diff2, real &norm2) {
  diff2 = norm2 = 0;                                            // Set unused values to 0
  B_iter B2 = bodies2.begin();                                  // Set iterator for bodies2
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B, ++B2 ) {  // Loop over bodies & bodies2
#ifdef DEBUG
    std::cout << B->I << " " << B->vel[0] << " " << B2->vel[0] << std::endl;// Compare every element
#endif
    velDiff += norm(B->vel - B2->vel);                          // Difference of velocity
    velNorm += norm(B2->vel);                                   // Value of velocity
  }                                                             // End loop over bodies & bodies2
}

void Dataset::printError(real velDiff, real velNorm, real diff2, real norm2) {// Print relative L2 norm error
  vect dummy = diff2; dummy = norm2;                            // Use the values so compiler does not complain
  std::cout << "Error         : " << std::sqrt(velDiff/velNorm) << std::endl;
}

void Kernel::P2P_CPU() {
  for( B_iter BI=BI0; BI!=BIN; ++BI ) {
    for( B_iter BJ=BJ0; BJ!=BJN; ++BJ ) {
      vect dist = BI->X - BJ->X - Xperiodic;
      real S2 = 2 * BJ->S * BJ->S;
      real R2  = norm(dist) + EPS2;
      real RS = R2 / S2;
      real cutoff = 0.25 / M_PI / R2 / std::sqrt(R2) * (erf( std::sqrt(RS) )
                  - std::sqrt(4 / M_PI * RS) * exp(-RS));
      BI->vel[0] += (dist[1] * BJ->Q[2] - dist[2] * BJ->Q[1]) * cutoff;
      BI->vel[1] += (dist[2] * BJ->Q[0] - dist[0] * BJ->Q[2]) * cutoff;
      BI->vel[2] += (dist[0] * BJ->Q[1] - dist[1] * BJ->Q[0]) * cutoff;
    }
  }
}
