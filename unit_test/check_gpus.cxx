#include "mympi.h"
#include "dataset.h"
#include "evaluator.h"

int main() {
  char hostname[256];
  double tic,toc;
  const int numBodies = 10000;
  tic = get_time();
  Bodies bodies(numBodies);
  Dataset D;
  Evaluator E;
  MyMPI M;
  E.initialize();
  gethostname(hostname,sizeof(hostname));
  bool print = true;
  if( M.commRank() != 0 ) print = false;
  toc = get_time();
  if(print) std::cout << "Allocate      : " << toc-tic << std::endl;

  tic = get_time();
  D.sphere(bodies);
  toc = get_time();
  if(print) std::cout << "Set bodies    : " << toc-tic << std::endl;

  tic = get_time();
  E.evalP2P(bodies,bodies);
  toc = get_time();
  if(print) std::cout << "Direct GPU    : " << toc-tic << std::endl;

  tic = get_time();
  real err = 0, rel = 0;
  for( B_iter BI=bodies.begin(); BI!=bodies.end(); ++BI ) {
    real pot = -BI->scal / std::sqrt(EPS2);
    BI->pot -= BI->scal / std::sqrt(EPS2);
    for( B_iter BJ=bodies.begin(); BJ!=bodies.end(); ++BJ ) {
      vect dist = BI->pos - BJ->pos;
      real r = std::sqrt(norm(dist) + EPS2);
      pot += BJ->scal / r;
    }
//    std::cout << BI-bodies.begin() << " " << BI->pot << " " << pot << std::endl;
    err += (BI->pot - pot) * (BI->pot - pot);
    rel += pot * pot;
  }
  toc = get_time();
  if(print) std::cout << "Direct CPU    : " << toc-tic << std::endl;

  for( int irank=0; irank!=M.commSize(); ++irank ) {
    MPI_Barrier(MPI_COMM_WORLD);
    if( M.commRank() == irank ) {
      std::cout << hostname << " @ rank : " << M.commRank() << " / " << M.commSize();
      std::cout << " @ device : " << M.commRank()%GPUS << " / " << GPUS << std::endl;
      std::cout << "Error         : " << std::sqrt(err/rel) << std::endl;
    }
    usleep(100);
  }
  E.finalize();
}
