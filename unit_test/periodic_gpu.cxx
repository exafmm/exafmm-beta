#include "dataset.h"
#include "evaluator.h"

int main() {
  const int numBodies = 1000;
  Bodies bodies(numBodies);
  Dataset D;
  Evaluator E;
  E.initialize();
  E.printNow = true;

  E.startTimer("Set bodies   ");
  D.sphere(bodies);
  E.setX0(0);
  E.setR0(1);
  E.stopTimer("Set bodies   ",E.printNow);

  E.startTimer("Direct GPU   ");
  E.evalP2P(bodies,bodies,true);
  E.stopTimer("Direct GPU   ",E.printNow);

  E.startTimer("Direct CPU   ");
  real potDiff = 0, potNorm = 0, accDiff = 0, accNorm = 0;
  for( B_iter BI=bodies.begin(); BI!=bodies.end(); ++BI ) {
    real pot = -BI->scal / std::sqrt(EPS2);
    vect acc = 0;
    BI->pot -= BI->scal / std::sqrt(EPS2);
    int I=0;
    for( int ix=-1; ix<=1; ++ix ) {
      for( int iy=-1; iy<=1; ++iy ) {
        for( int iz=-1; iz<=1; ++iz, ++I ) {
          for( B_iter BJ=bodies.begin(); BJ!=bodies.end(); ++BJ ) {
            vect dist = BI->pos - BJ->pos;
            dist[0] -= ix * 2;
            dist[1] -= iy * 2;
            dist[2] -= iz * 2;
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
