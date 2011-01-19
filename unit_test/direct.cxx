#include "bottomup.h"
#include "kernel.h"

int main()
{
  double tic,toc;
  int const Nbody=10000;
  tic = get_time();
  Bodies bodies(Nbody);
  Kernels K;
  toc = get_time();
  std::cout << "Allocate      : " << toc-tic << std::endl;

  tic = get_time();
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {        // Loop over all bodies
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
  K.direct(bodies.begin(),bodies.end());
  toc = get_time();
  std::cout << "Direct sum    : " << toc-tic << std::endl;

  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
//    std::cout << B-bodies.begin() << " " << B->acc << std::endl;
  }

}
