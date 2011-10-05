#include "dataset.h"
#include "tree.h"

int main() {
#ifdef MANY
  for ( int it=0; it<25; it++ ) {
#else
#if BUILD
  for ( int it=32; it<33; it++ ) {
#else
  for ( int it=8; it<9; it++ ) {
#endif
#endif
  int numBodies = int(pow(10,(it+24)/8.0));
  std::cout << "N             : " << numBodies << std::endl;
  IMAGES = 0;
  THETA = 0.6;
  Bodies bodies(numBodies);
  Bodies bodies2;
  Dataset D;
  D.kernelName = "Laplace";
  TreeConstructor T;
  D.random(bodies);
  T.startTimer("FMM          ");
  T.setDomain(bodies);
  T.build();
  T.link();
#if BUILD
  approx = 0;
#else
  T.approximate(bodies);
  T.stopTimer("FMM          ",true);
  T.eraseTimer("FMM          ");

#ifdef DIRECT
  bodies2 = bodies;
  T.startTimer("Direct sum   ");
  T.exact(bodies2);
  T.stopTimer("Direct sum   ",true);
  T.eraseTimer("Direct sum   ");
  T.writeTime();

  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
  D.evalError(bodies,bodies2,diff1,norm1,diff2,norm2);
  D.printError(diff1,norm1,diff2,norm2);

/*
  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
  for( B_iter B=bodies.begin(), B2=bodies2.begin(); B!=bodies.end(); ++B, ++B2 ) {// Loop over bodies
    diff1 += (B->TRG[0] - B2->TRG[0]) * (B->TRG[0] - B2->TRG[0]);// Difference of potential
    norm1 += B2->TRG[0] * B2->TRG[0];                           //  Value of potential
    diff2 += (B->TRG[1] - B2->TRG[1]) * (B->TRG[1] - B2->TRG[1]);// Difference of x acceleration
    diff2 += (B->TRG[2] - B2->TRG[2]) * (B->TRG[2] - B2->TRG[2]);// Difference of y acceleration
    diff2 += (B->TRG[3] - B2->TRG[3]) * (B->TRG[3] - B2->TRG[3]);// Difference of z acceleration
    norm2 += B2->TRG[1] * B2->TRG[1];                           //  Value of x acceleration
    norm2 += B2->TRG[2] * B2->TRG[2];                           //  Value of y acceleration
    norm2 += B2->TRG[3] * B2->TRG[3];                           //  Value of z acceleration
  }
  std::cout << "Error (pot)   : " << std::sqrt(diff1/norm1) << std::endl;
  std::cout << "Error (acc)   : " << std::sqrt(diff2/norm2) << std::endl;
*/
#endif
#endif
  }
}
