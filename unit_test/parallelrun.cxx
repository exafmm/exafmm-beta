#include "let.h"
#include "dataset.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  double tic,toc;
  int const numBodies(10000);
  tic = get_time();
  Bodies bodies(numBodies);
  Cells cells;
  Dataset D;
  LocalEssentialTree T;
  Evaluator E;
  bool print(true);
  if( T.commRank() != 0 ) print = false;
  toc = get_time();
  if(print) std::cout << "Allocate      : " << toc-tic << std::endl;

  tic = get_time();
  D.sphere(bodies,T.commRank()+1);
  toc = get_time();
  if(print) std::cout << "Set bodies    : " << toc-tic << std::endl;

  tic = get_time();
  T.setGlobDomain(bodies);
  toc = get_time();
  if(print) std::cout << "Set domain    : " << toc-tic << std::endl;

  tic = get_time();
  T.bisection(bodies);
  toc = get_time();
  if(print) std::cout << "Partition     : " << toc-tic << std::endl;

#ifdef TOPDOWN
  T.topdown(bodies,cells,print);
#else
  T.bottomup(bodies,cells,print);
#endif

  tic = get_time();
  T.commBodies(cells);
  toc = get_time();
  if(print) std::cout << "Comm bodies   : " << toc-tic << std::endl;

  tic = get_time();
  Bodies bodies2 = bodies;
  Cells jcells = cells;
  T.commCells(bodies2,jcells);
  toc = get_time();
  if(print) std::cout << "Comm cells    : " << toc-tic << std::endl;

  tic = get_time();
  T.downward(cells,jcells,1);
  toc = get_time();
  if(print) std::cout << "Downward      : " << toc-tic << std::endl;

  tic = get_time();
  bodies2 = bodies;
  for( B_iter B=bodies2.begin(); B!=bodies2.end(); ++B ) {
    B->pot = -B->scal / std::sqrt(EPS2);
  }
  for( int i=0; i!=T.commSize(); ++i ) {
    T.shiftBodies(bodies);
    E.evalP2P(bodies2,bodies);
    if(print) std::cout << "Direct loop   : " << i+1 << "/" << T.commSize() << std::endl;
  }
  toc = get_time();
  if(print) std::cout << "Direct sum    : " << toc-tic << std::endl;

  B_iter B  = bodies.begin();
  B_iter B2 = bodies2.begin();
  real err(0),rel(0),err2,rel2;
  for( int i=0; i!=int(bodies.size()); ++i,++B,++B2 ) {
    B->pot  -= B->scal  / std::sqrt(EPS2);
#ifdef DEBUG
    if(MPIRANK==0) std::cout << B->I << " " << B->pot << " " << B2->pot << std::endl;
#endif
    err += (B->pot - B2->pot) * (B->pot - B2->pot);
    rel += B2->pot * B2->pot;
  }
  int MPI_TYPE = T.getType(err);
  MPI_Reduce(&err,&err2,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  MPI_Reduce(&rel,&rel2,1,MPI_TYPE,MPI_SUM,0,MPI_COMM_WORLD);
  if(print) std::cout << "Error         : " << std::sqrt(err2/rel2) << std::endl;
#ifdef DEBUG
  T.print(std::sqrt(err/rel));
#endif

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
  if( T.commRank() == 0 ) {
    vtk.setDomain(T.getR0(),T.getX0());
    vtk.setGroupOfPoints(bodies,Ncell);
  }
  for( int i=1; i!=T.commSize(); ++i ) {
    T.shiftBodies(bodies);
    if( T.commRank() == 0 ) {
      vtk.setGroupOfPoints(bodies,Ncell);
    }
  }
  if( T.commRank() == 0 ) {
    vtk.plot(Ncell);
  }
#endif
}
