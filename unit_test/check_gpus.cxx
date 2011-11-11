#include "mympi.h"
#include "dataset.h"
#include "evaluator.h"

int main() {
  char hostname[256];
  const int numBodies = 1000;
  IMAGES = 0;
  THETA = 1/sqrtf(3);
  Bodies bodies(numBodies);
  Dataset dataset;
  dataset.kernelName = "Laplace";
  Evaluator FMM;
  MyMPI M;
  FMM.setKernel(dataset.kernelName);
  FMM.initialize();
  FMM.preCalculation();
  gethostname(hostname,sizeof(hostname));
  if( MPIRANK == 0 ) FMM.printNow = true;

  FMM.startTimer("Set bodies   ");
  dataset.sphere(bodies);
  FMM.stopTimer("Set bodies   ",FMM.printNow);

  FMM.startTimer("Set domain   ");
  FMM.setDomain(bodies);
  FMM.stopTimer("Set domain   ",FMM.printNow);

  FMM.startTimer("Direct GPU   ");
  FMM.evalP2P(bodies,bodies);
  FMM.stopTimer("Direct GPU   ",FMM.printNow);

  FMM.startTimer("Direct CPU   ");
  bool onCPU = true;
  Bodies bodies2 = bodies;
  dataset.initTarget(bodies2);
  FMM.evalP2P(bodies2,bodies2,onCPU);
  FMM.stopTimer("Direct CPU   ",FMM.printNow);

  real diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
  dataset.evalError(bodies,bodies2,diff1,norm1,diff2,norm2);

  for( int irank=0; irank!=MPISIZE; ++irank ) {
    MPI_Barrier(MPI_COMM_WORLD);
    if( MPIRANK == irank ) {
      std::cout << hostname << " @ rank : " << MPIRANK << " / " << MPISIZE;
      std::cout << " @ device : " << DEVICE << " / " << GPUS << std::endl;
      std::cout << "Error         : " << std::sqrt(diff1/norm1) << std::endl;
    }
    usleep(100);
  }
  FMM.postCalculation();
  FMM.finalize();
}
