#include "bottomup.h"
#ifdef VTK
#include "vtk.h"
#endif

int main()
{
  double tic,toc;
  int const Nbody=1000;
  Bodies bodies(Nbody);
  Bodies bodies2(Nbody);
  B_iter B;
  BottomUpTreeConstructor T(bodies);

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
  T.setMorton();
  toc = get_time();
  std::cout << "Get Morton    : " << toc-tic << std::endl;

  tic = get_time();
  T.sort(T.Ibody,bodies,bodies2);
  toc = get_time();
  std::cout << "Sort Morton   : " << toc-tic << std::endl;

  tic = get_time();
  T.prune();
  toc = get_time();
  std::cout << "Prune tree    : " << toc-tic << std::endl;

  tic = get_time();
  T.grow(bodies2);
  toc = get_time();
  std::cout << "Grow tree     : " << toc-tic << std::endl;

  tic = get_time();
  T.sort(T.Ibody,bodies,bodies2,false);
  toc = get_time();
  std::cout << "Sort descend  : " << toc-tic << std::endl;

  tic = get_time();
  T.link();
  toc = get_time();
  std::cout << "Link cells    : " << toc-tic << std::endl;

#ifdef VTK
  T.sort(T.Ibody,bodies,bodies2);
  int Ncell(0);
  vtkPlot vtk;
  vtk.setDomain(T.getR0(),T.getX0());
  vtk.setGroupOfPoints(T.Ibody,bodies,Ncell);
  vtk.plot(Ncell);
#endif
}
