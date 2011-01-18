#include "topdown.h"
#ifdef VTK
#include "vtk.h"
#endif

int main()
{
  double tic,toc;
  int Nbody=50;
  Bodies bodies(Nbody);
  B_iter B;
  TopDownTreeConstructor T(bodies);

  tic = get_time();
  for( B=bodies.begin(); B!=bodies.end(); ++B ) {               // Loop over all bodies
    for( int d=0; d!=3; ++d )                                   //  Loop over each dimension
      B->pos[d] = rand()/(1.+RAND_MAX)*2-1;                     //   Initialize positions
    real r = sqrt(B->pos[0]*B->pos[0]+B->pos[1]*B->pos[1]+B->pos[2]*B->pos[2]);
    for( int d=0; d!=3; ++d )                                   //  Loop over each dimension
      B->pos[d] /= r*1.1;                                       //   Normalize positions
    B->scal = 1./bodies.size();                                 //  Initialize source value
  }                                                             // End loop over all bodies
  toc = get_time();
  std::cout << "Initialize    : " << toc-tic << std::endl;

  tic = get_time();
  T.setDomain();
  toc = get_time();
  std::cout << "Set domain    : " << toc-tic << std::endl;

  tic = get_time();
  T.grow();
  toc = get_time();
  std::cout << "Grow tree     : " << toc-tic << std::endl;
}
