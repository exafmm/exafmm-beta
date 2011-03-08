#include "dataset.h"
#include "evaluator.h"

int main() {
  const int numBodies = 10000;
  Bodies bodies(numBodies);
  Dataset D;
  Evaluator E;
  E.initialize();

  E.startTimer("Set bodies   ");
  D.sphere(bodies);
  E.stopTimer("Set bodies   ",true);

  E.startTimer("Direct GPU   ");
  E.evalP2P(bodies,bodies);
  E.stopTimer("Direct GPU   ",true);

  E.startTimer("Direct CPU   ");
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
  E.stopTimer("Direct CPU   ",true);
  std::cout << "Error         : " << std::sqrt(err/rel) << std::endl;
  E.finalize();
}
