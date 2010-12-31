#include "kernel.h"

int main()
{
  double tic,toc;
  int const Nbody=10000;
  Bodies bodies(Nbody);
  B_iter B;
  Kernels K(bodies);

  tic = get_time();
  for( B=bodies.begin(); B!=bodies.end(); ++B ) {               // Loop over all bodies
    for( int d=0; d!=3; ++d ) {                                 //  Loop over each dimension
      B->pos[d] = rand()/(1.+RAND_MAX)*2-1;                     //   Initialize positions
      B->acc[d] = 0;
    }
    real r = sqrt(B->pos[0]*B->pos[0]+B->pos[1]*B->pos[1]+B->pos[2]*B->pos[2]);
    for( int d=0; d!=3; ++d )                                   //  Loop over each dimension
      B->pos[d] /= r;                                           //   Normalize positions
    B->scal = 1./bodies.size();                                 //  Initialize source value
    B->pot = 0;
  }                                                             // End loop over all bodies
  toc = get_time();
  std::cout << "Initialize    : " << toc-tic << std::endl;

  tic = get_time();
  K.direct();
  toc = get_time();
  std::cout << "Direct sum    : " << toc-tic << std::endl;

  for( B=bodies.begin(); B!=bodies.end(); ++B ) {
//    std::cout << B-bodies.begin() << " " << B->acc << std::endl;
  }

}
