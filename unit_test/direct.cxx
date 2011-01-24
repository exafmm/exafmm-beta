#include "dataset.h"
#include "kernel.h"

int main() {
  double tic,toc;
  int const numBodies(10000);
  tic = get_time();
  Bodies bodies(numBodies);
  Dataset D(bodies);
  Kernel K;
  toc = get_time();
  std::cout << "Allocate      : " << toc-tic << std::endl;

  tic = get_time();
  D.sphere();
  toc = get_time();
  std::cout << "Set bodies    : " << toc-tic << std::endl;

  tic = get_time();
  K.P2P(bodies.begin(),bodies.end());
  toc = get_time();
  std::cout << "Direct sum    : " << toc-tic << std::endl;

  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
//    std::cout << B-bodies.begin() << " " << B->acc << std::endl;
  }

}
