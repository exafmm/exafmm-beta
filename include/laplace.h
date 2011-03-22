#include "dataset.h"

void Dataset::initSource(Bodies &bodies) {                      // Initialize source values
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {        // Loop over bodies
    B->Q = 1. / bodies.size() / MPISIZE;                        //  Initialize mass/charge
  }                                                             // End loop over bodies
}

void Dataset::initTarget(Bodies &bodies, bool IeqJ) {           // Initialize target values
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {        // Loop over bodies
    B->pot = -B->Q / std::sqrt(EPS2) * IeqJ;                    //  Initialize potential
    B->acc = 0;                                                 //  Initialize acceleration
  }                                                             // End loop over bodies
}

void Dataset::readTarget(Bodies &bodies) {
  std::ifstream file("data",std::ios::in | std::ios::ate | std::ios::binary);
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    file >> B->pot;
    file >> B->acc[0];
    file >> B->acc[1];
    file >> B->acc[2];
  }
  file.close();
}

void Dataset::writeTarget(Bodies &bodies) {
  std::ofstream file("data",std::ios::out | std::ios::app | std::ios::binary);
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
    file << B->pot << std::endl;
    file << B->acc[0] << std::endl;
    file << B->acc[1] << std::endl;
    file << B->acc[2] << std::endl;
  }
  file.close();
}

void Dataset::evalError(Bodies &bodies, Bodies &bodies2,        // Evaluate error
                        real &potDiff, real &potNorm, real &accDiff, real &accNorm) {
  B_iter B2 = bodies2.begin();                                  // Set iterator for bodies2
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B, ++B2 ) {  // Loop over bodies & bodies2
#ifdef DEBUG
    std::cout << B->I << " " << B->pot << " " << B2->pot << std::endl;// Compare every element
#endif
    potDiff += (B->pot - B2->pot) * (B->pot - B2->pot);         // Difference of potential
    potNorm += B2->pot * B2->pot;                               // Value of potential
    accDiff += norm(B->acc - B2->acc);                          // Difference of acceleration
    accNorm += norm(B2->acc);                                   // Value of acceleration
  }                                                             // End loop over bodies & bodies2
}

void Dataset::printError(real potDiff, real potNorm, real accDiff, real accNorm) {// Print relative L2 norm error
  std::cout << "Error (pot)   : " << std::sqrt(potDiff/potNorm) << std::endl;
  std::cout << "Error (acc)   : " << std::sqrt(accDiff/accNorm) << std::endl;
}

void Kernel::P2P_CPU() {
  for( B_iter BI=BI0; BI!=BIN; ++BI ) {
    for( B_iter BJ=BJ0; BJ!=BJN; ++BJ ) {
      vect dist = BI->X - BJ->X - Xperiodic;
      real invR = 1 / std::sqrt(norm(dist) + EPS2);
      real invR3 = BJ->Q * invR * invR * invR;
      BI->pot += BJ->Q * invR;
      BI->acc -= dist * invR3;
    }
  }
}
