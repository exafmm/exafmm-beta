#include "dataset.h"
#include "kernel.h"

int main() {
  double tic,toc;
  int const numBodies(10000);
  tic = get_time();
  Bodies bodies(numBodies);
  Dataset D;
  Kernel K;
  K.setup();
  toc = get_time();
  std::cout << "Allocate      : " << toc-tic << std::endl;

  tic = get_time();
  D.sphere(bodies);
  toc = get_time();
  std::cout << "Set bodies    : " << toc-tic << std::endl;

  tic = get_time();
  K.P2P(bodies.begin(),bodies.end());
  toc = get_time();
  std::cout << "Direct GPU    : " << toc-tic << std::endl;

  tic = get_time();
  real err(0),rel(0);
  for( B_iter Bi=bodies.begin(); Bi!=bodies.end(); ++Bi ) {
    real pot = -Bi->scal / std::sqrt(EPS2);
    for( B_iter Bj=bodies.begin(); Bj!=bodies.end(); ++Bj ) {
      vect dist = Bi->pos - Bj->pos;
      real r = std::sqrt(norm(dist) + EPS2);
      pot += Bj->scal / r;
    }
    err += (Bi->pot - pot) * (Bi->pot - pot);
    rel += pot * pot;
  }
  toc = get_time();
  std::cout << "Direct CPU    : " << toc-tic << std::endl;
  std::cout << "Error         : " << std::sqrt(err/rel) << std::endl;
}
