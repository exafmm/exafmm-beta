#include "dataset.h"
#include "evaluator.h"

int main() {
  const int numBodies = 1000;
  Bodies bodies(numBodies);
  Bodies jbodies;
  Dataset D;
  Evaluator E;
  E.initialize();
  E.printNow = true;

  E.startTimer("Set bodies   ");
  D.sphere(bodies);
  E.stopTimer("Set bodies   ",E.printNow);

  E.startTimer("Set domain   ");
  E.setDomain(bodies);
  E.stopTimer("Set domain   ",E.printNow);

  if( IMAGES != 0 ) {
    E.startTimer("Set periodic ");
    jbodies = E.periodicBodies(bodies);
    E.stopTimer("Set periodic ",E.printNow);
  } else {
    jbodies = bodies;
  }

  E.startTimer("Direct GPU   ");
  E.evalP2P(bodies,jbodies);
  E.stopTimer("Direct GPU   ",E.printNow);

  E.startTimer("Direct CPU   ");
  int prange = E.getPeriodicRange();
  real potDiff = 0, potNorm = 0, accDiff = 0, accNorm = 0;
  for( B_iter BI=bodies.begin(); BI!=bodies.end(); ++BI ) {
    real pot = -BI->scal / std::sqrt(EPS2);
    vect acc = 0;
    BI->pot -= BI->scal / std::sqrt(EPS2);
    for( int ix=-prange; ix<=prange; ++ix ) {
      for( int iy=-prange; iy<=prange; ++iy ) {
        for( int iz=-prange; iz<=prange; ++iz ) {
          for( B_iter BJ=bodies.begin(); BJ!=bodies.end(); ++BJ ) {
            vect dist = BI->X - BJ->X;
            dist[0] -= ix * 2 * E.getR0();
            dist[1] -= iy * 2 * E.getR0();
            dist[2] -= iz * 2 * E.getR0();
            real invR = 1 / std::sqrt(norm(dist) + EPS2);
            real invR3 = BJ->scal * invR * invR * invR;
            pot += BJ->scal * invR;
            acc -= dist * invR3;
          }
        }
      }
    }
//    std::cout << BI-bodies.begin() << " " << BI->pot << " " << pot << std::endl;
    potDiff += (BI->pot - pot) * (BI->pot - pot);
    potNorm += pot * pot;
    accDiff += norm(BI->acc - acc);
    accNorm += norm(acc);
  }
  E.stopTimer("Direct CPU   ",E.printNow);
  std::cout << "Error (pot)   : " << std::sqrt(potDiff/potNorm) << std::endl;
  std::cout << "Error (acc)   : " << std::sqrt(accDiff/accNorm) << std::endl;
  E.finalize();
}
