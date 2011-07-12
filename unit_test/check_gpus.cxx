#include "mympi.h"
#include "dataset.h"
#include "evaluator.h"

int main() {
  char hostname[256];
  const int numBodies = 1000;
  std::string kernelName = "Laplace";
  IMAGES = 0;
  THETA = 1/sqrtf(3);
  Bodies bodies(numBodies);
  Dataset D;
  Evaluator E;
  MyMPI M;
  E.setKernel(kernelName);
  E.initialize();
  D.kernelName = kernelName;
  E.preCalculation();
  gethostname(hostname,sizeof(hostname));
  if( MPIRANK == 0 ) E.printNow = true;

  E.startTimer("Set bodies   ");
  D.sphere(bodies);
  E.stopTimer("Set bodies   ",E.printNow);

  E.startTimer("Set domain   ");
  E.setDomain(bodies);
  E.stopTimer("Set domain   ",E.printNow);

  E.startTimer("Direct GPU   ");
  E.evalP2P(bodies,bodies);
  E.stopTimer("Direct GPU   ",E.printNow);

  E.startTimer("Direct CPU   ");
  bool onCPU = true;
  Bodies bodies2 = bodies;
  D.initTarget(bodies2);
  E.evalP2P(bodies2,bodies2,onCPU);
  E.stopTimer("Direct CPU   ",E.printNow);

  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
  D.evalError(bodies,bodies2,diff1,norm1,diff2,norm2);

  for( int irank=0; irank!=MPISIZE; ++irank ) {
    MPI_Barrier(MPI_COMM_WORLD);
    if( MPIRANK == irank ) {
      std::cout << hostname << " @ rank : " << MPIRANK << " / " << MPISIZE;
      std::cout << " @ device : " << DEVICE << " / " << GPUS << std::endl;
      std::cout << "Error         : " << std::sqrt(diff1/norm1) << std::endl;
    }
    usleep(100);
  }
  E.postCalculation();
  E.finalize();
}
