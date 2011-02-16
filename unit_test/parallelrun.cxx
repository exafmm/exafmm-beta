#include "let.h"
#include "dataset.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  double tic,toc;
  int const numBodies(4096);
  tic = get_time();
  Bodies bodies(numBodies);
  Cells cells;
  Dataset D;
  LocalEssentialTree P;
  Kernel K;
  bool print(true);
  if( P.commRank() != 0 ) print = false;
  toc = get_time();
  if(print) std::cout << "Allocate      : " << toc-tic << std::endl;

  tic = get_time();
  srand(P.commRank()+1);
  D.random(bodies);
  toc = get_time();
  if(print) std::cout << "Set bodies    : " << toc-tic << std::endl;

  tic = get_time();
  P.setGlobDomain(bodies);
  toc = get_time();
  if(print) std::cout << "Set domain    : " << toc-tic << std::endl;

  tic = get_time();
  P.bisection(bodies);
  toc = get_time();
  if(print) std::cout << "Partition     : " << toc-tic << std::endl;

#ifdef TOPDOWN
  P.topdown(bodies,cells,print);
#else
  P.bottomup(bodies,cells,print);
#endif

  tic = get_time();
  P.setCI0(cells.begin());
  P.setCJ0(cells.begin());
  P.P2M(cells);
  toc = get_time();
  if(print) std::cout << "P2M           : " << toc-tic << std::endl;

  tic = get_time();
  P.M2M(cells);
  toc = get_time();
  if(print) std::cout << "M2M           : " << toc-tic << std::endl;

  tic = get_time();
  P.commBodies(cells);
  toc = get_time();
  if(print) std::cout << "Comm bodies   : " << toc-tic << std::endl;

  tic = get_time();
  Bodies bodies2 = bodies;
  Cells jcells = cells;
  P.commCells(bodies2,jcells);
  P.setCI0(cells.begin());
  P.setCJ0(jcells.begin());
  toc = get_time();
  if(print) std::cout << "Comm cells    : " << toc-tic << std::endl;

  tic = get_time();
  P.evaluate(cells,jcells,1);
  toc = get_time();
  if(print) std::cout << "Evaluate      : " << toc-tic << std::endl;

  tic = get_time();
  bodies2 = bodies;
  for( B_iter B=bodies2.begin(); B!=bodies2.end(); ++B ) {
    B->acc = B->pot = 0;
  }
  for( int i=0; i!=P.commSize(); ++i ) {
    P.shiftBodies(bodies);
    K.P2P(bodies2.begin(),bodies2.end(),bodies.begin(),bodies.end());
    if(print) std::cout << "Direct loop   : " << i+1 << "/" << P.commSize() << std::endl;
  }
  toc = get_time();
  if(print) std::cout << "Direct sum    : " << toc-tic << std::endl;

  B_iter B  = bodies.begin();
  B_iter B2 = bodies2.begin();
  real err(0),rel(0),err2,rel2;
  for( int i=0; i!=int(bodies.size()); ++i,++B,++B2 ) {
    B->pot  -= B->scal  / std::sqrt(EPS2);                      //  Initialize body values
    B2->pot -= B2->scal / std::sqrt(EPS2);                      //  Initialize body values
//    if(MPIRANK==0) std::cout << i << " " << B->pot << " " << B2->pot << std::endl;
    err += (B->pot - B2->pot) * (B->pot - B2->pot);
    rel += B2->pot * B2->pot;
  }
  int MPI_TYPE = P.getType(err);
  MPI_Reduce(&err,&err2,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&rel,&rel2,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  if(print) std::cout << "Error         : " << std::sqrt(err2/rel2) << std::endl;
  P.print(std::sqrt(err/rel));

#ifdef VTK
  for( B=bodies.begin(); B!=bodies.end(); ++B ) B->I = 0;
  for( C_iter C=jcells.begin(); C!=jcells.end(); ++C ) {
    Body body;
    body.I = 1;
    body.pos  = C->X;
    body.scal = 0;
    bodies.push_back(body);
  }

  int Ncell(0);
  vtkPlot vtk;
  if( P.commRank() == 0 ) {
    vtk.setDomain(P.getR0(),P.getX0());
    vtk.setGroupOfPoints(bodies,Ncell);
  }
  for( int i=1; i!=P.commSize(); ++i ) {
    P.shiftBodies(bodies);
    if( P.commRank() == 0 ) {
      vtk.setGroupOfPoints(bodies,Ncell);
    }
  }
  if( P.commRank() == 0 ) {
    vtk.plot(Ncell);
  }
#endif
}
