#include <fstream>
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
  double tic,toc,tree,approx;
  Bodies bodies(numBodies);
  Bodies bodies2;
  srand(1);
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {        // Loop over bodies
    for( int d=0; d!=3; ++d ) {                                 //  Loop over dimension
      B->X[d] = rand() / (1. + RAND_MAX);                       //   Initialize positions
    }                                                           //  End loop over dimension
    B->SRC[0] = 1. / bodies.size();                             //  Initialize mass/charge
  }                                                             // End loop over bodies

  tic = get_time();
  TreeConstructor T;
  T.setDomain(bodies);
  T.build();
  toc = get_time();
  std::cout << "build  : " << toc-tic << std::endl;
  tree = toc-tic;
  tic = get_time();
  T.link();
  toc = get_time();
  std::cout << "link   : " << toc-tic << std::endl;
  tree += toc-tic;
#if BUILD
  approx = 0;
#else
  tic = get_time();
  T.approximate(bodies);
  toc = get_time();
  approx = toc-tic;

#ifdef DIRECT
  bodies2 = bodies;
  tic = get_time();
  T.exact(bodies2);
  toc = get_time();
  std::cout << "direct time   : " << toc-tic << std::endl;
#endif
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
#ifdef DIRECT
  std::cout << "FMM time      : " << tree+approx << std::endl;
  std::cout << "Error (pot)   : " << std::sqrt(diff1/norm1) << std::endl;
  std::cout << "Error (acc)   : " << std::sqrt(diff2/norm2) << std::endl;
#else
  std::cout << tree+approx << std::endl;
#endif
  std::ofstream file("time.dat",std::ios::out | std::ios::app);
  file << tree+approx << std::endl;
  file.close();
#endif
  }
}
