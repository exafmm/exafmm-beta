#include "let.h"
#include "dataset.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  double tic,toc;
  int const numBodies(1000);
  tic = get_time();
  Bodies bodies(numBodies);
  Dataset D(bodies);
  LocalEssentialTree P(bodies);
  Kernel K;
  bool print(true);
  if( P.commRank() != 0 ) print = false;
  toc = get_time();
  if(print) std::cout << "Allocate      : " << toc-tic << std::endl;

  tic = get_time();
  srand(P.commRank()+1);
  D.random();
  toc = get_time();
  if(print) std::cout << "Set bodies    : " << toc-tic << std::endl;

  tic = get_time();
  P.setGlobDomain();
  toc = get_time();
  if(print) std::cout << "Set domain    : " << toc-tic << std::endl;

  tic = get_time();
  P.bisection();
  toc = get_time();
  if(print) std::cout << "Partition     : " << toc-tic << std::endl;

#ifdef TOPDOWN
  P.topdown(print);
#else
  P.bottomup(print);
#endif

  tic = get_time();
  P.P2M();
  toc = get_time();
  if(print) std::cout << "P2M           : " << toc-tic << std::endl;

  tic = get_time();
  P.M2M();
  toc = get_time();
  if(print) std::cout << "M2M           : " << toc-tic << std::endl;

  tic = get_time();
  P.commBodies();
  toc = get_time();
  if(print) std::cout << "Comm bodies   : " << toc-tic << std::endl;

  tic = get_time();
  P.commCells();
  toc = get_time();
  if(print) std::cout << "Comm cells    : " << toc-tic << std::endl;

  tic = get_time();
  P.evaluate(1);
  toc = get_time();
  if(print) std::cout << "Evaluate      : " << toc-tic << std::endl;

  tic = get_time();
  srand(P.commRank()+1);
  bodies.resize(numBodies);
  D.random();
  P.buffer = bodies;
  for( int i=0; i!=P.commSize(); ++i ) {
    P.shiftBodies();
    K.P2P(P.buffer.begin(),P.buffer.end(),bodies.begin(),bodies.end());
  }
  toc = get_time();
  if(print) std::cout << "Direct sum    : " << toc-tic << std::endl;

  B_iter B  = bodies.begin();
  B_iter B2 = P.buffer.begin();
  real err(0),rel(0);
  for( int i=0; i!=numBodies; ++i,++B,++B2 ) {
    B->pot  -= B->scal  / std::sqrt(EPS2);                      //  Initialize body values
    B2->pot -= B2->scal / std::sqrt(EPS2);                      //  Initialize body values
    err += (B->pot - B2->pot) * (B->pot - B2->pot);
    rel += B2->pot * B2->pot;
  }
  if(print) std::cout << "Error         : " << std::sqrt(err/rel) << std::endl;

#ifdef VTK
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) B->I = 0;
  for( B_iter B=P.buffer.begin(); B!=P.buffer.end(); ++B ) B->I = 1;
  bodies.insert(bodies.end(),P.buffer.begin(),P.buffer.end());

  int Ncell(0);
  vtkPlot vtk;
  if( P.commRank() == 0 ) {
    vtk.setDomain(P.getR0(),P.getX0());
    vtk.setGroupOfPoints(bodies,Ncell);
  }
  tic = get_time();
  for( int i=1; i!=P.commSize(); ++i ) {
    P.shiftBodies();
    if( P.commRank() == 0 ) {
      vtk.setGroupOfPoints(bodies,Ncell);
    }
  }
  toc = get_time();
  if(print) std::cout << "Shift bodies  : " << toc-tic << std::endl;
  if( P.commRank() == 0 ) {
    vtk.plot(Ncell);
  }
#endif
}
