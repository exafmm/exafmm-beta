#include "mympi.h"
#include "dataset.h"
#include "evaluator.h"

int main() {
  char hostname[256];
  const int numBodies = 10000;
  Bodies bodies(numBodies);
  Dataset D;
  Evaluator E;
  MyMPI M;
  E.initialize();
  gethostname(hostname,sizeof(hostname));
  if( M.commRank() == 0 ) E.printNow = true;

  E.startTimer("Set bodies   ");
  D.sphere(bodies);
  E.stopTimer("Set bodies   ",E.printNow);

  E.startTimer("Direct GPU   ");
  E.evalP2P(bodies,bodies);
  E.stopTimer("Direct GPU   ",E.printNow);

  E.startTimer("Direct CPU   ");
  real potDiff = 0, potNorm = 0, accDiff = 0, accNorm = 0;
  for( B_iter BI=bodies.begin(); BI!=bodies.end(); ++BI ) {
    real pot = -BI->scal / std::sqrt(EPS2);
    vect acc = 0;
    BI->pot -= BI->scal / std::sqrt(EPS2);
    for( B_iter BJ=bodies.begin(); BJ!=bodies.end(); ++BJ ) {
      vect dist = BI->pos - BJ->pos;
      real invR = 1 / std::sqrt(norm(dist) + EPS2);
      real invR3 = BJ->scal * invR * invR * invR;
      pot += BJ->scal * invR;
      acc -= dist * invR3;
      
    }
//    std::cout << BI-bodies.begin() << " " << BI->pot << " " << pot << std::endl;
    potDiff += (BI->pot - pot) * (BI->pot - pot);
    potNorm += pot * pot;
    accDiff += norm(BI->acc - acc);
    accNorm += norm(acc);
  }
  E.stopTimer("Direct CPU   ",E.printNow);

  for( int irank=0; irank!=M.commSize(); ++irank ) {
    MPI_Barrier(MPI_COMM_WORLD);
    if( M.commRank() == irank ) {
      std::cout << hostname << " @ rank : " << M.commRank() << " / " << M.commSize();
      std::cout << " @ device : " << M.commRank()%GPUS << " / " << GPUS << std::endl;
      std::cout << "Error (pot)   : " << std::sqrt(potDiff/potNorm) << std::endl;
      std::cout << "Error (acc)   : " << std::sqrt(accDiff/accNorm) << std::endl;
    }
    usleep(100);
  }
  E.finalize();
}
