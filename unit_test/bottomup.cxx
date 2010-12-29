#include "bottomup.h"
#ifdef VTK
#include "vtk.h"
#endif

int main()
{
  double tic,toc;
  int const Nbody=100;
  bodies B(Nbody,Nbody);
  bodies B2(Nbody,Nbody);
  BottomUpTreeConstructor T(B);

  tic = get_time();
  for( B=B.begin(); B!=B.end(); ++B ) {                         // Loop over all bodies
    for( int d=0; d!=3; ++d )                                   //  Loop over each dimension
      B.pos()[d] = rand()/(1.+RAND_MAX)*2-1;                    //   Initialize positions
    real r = sqrt(B.pos()[0]*B.pos()[0]+B.pos()[1]*B.pos()[1]+B.pos()[2]*B.pos()[2]);
    for( int d=0; d!=3; ++d )                                   //  Loop over each dimension
      B.pos()[d] /= r*1.1;                                      //   Normalize positions
    B.scal() = 1./B.size();                                     //  Initialize source value
  }                                                             // End loop over all bodies
  toc = get_time();
  std::cout << "Initialize    : " << toc-tic << std::endl;

  tic = get_time();
  T.setDomain();
  toc = get_time();
  std::cout << "Set domain    : " << toc-tic << std::endl;

  tic = get_time();
  bigint *index = new bigint [Nbody];
  T.setMorton(index);
  toc = get_time();
  std::cout << "Get Morton    : " << toc-tic << std::endl;

  tic = get_time();
  T.sort(index,B,B2,true);
  toc = get_time();
  std::cout << "Sort Morton   : " << toc-tic << std::endl;

  tic = get_time();
  T.prune(index);
  toc = get_time();
  std::cout << "Prune tree    : " << toc-tic << std::endl;

  tic = get_time();
  T.grow(index,B2);
  toc = get_time();
  std::cout << "Grow tree     : " << toc-tic << std::endl;

  tic = get_time();
  T.sort(index,B,B2,false);
  toc = get_time();
  std::cout << "Sort descend  : " << toc-tic << std::endl;

  T.allocateCells();
  T.link(index);
  T.deallocateCells();

#ifdef VTK
  T.sort(index,B,B2,true);
  int Ncell(0);
  vtkPlot vtk;
  vtk.setDomain(T.getR0(),T.getX0());
  vtk.setGroupOfPoints(index,B,Ncell);
  vtk.plot(Ncell);
#endif
  T.sortDealloc();
  delete[] index;
}
