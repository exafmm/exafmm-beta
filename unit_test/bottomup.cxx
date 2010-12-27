#include "body.h"
#include "tree.h"

class cell {
  int NLEAF;
  int NCHILD;
  bigint I;
  bodies *LEAF;
  cell *PARENT;
  cell *CHILD[8];
  coef *M,*L;
};

int main()
{
  double tic,toc;
  int const Nbody=50;
  bodies B(Nbody,Nbody);
  Tree tree(B);

  tic = get_time();
  for( B=B.begin(); B!=B.end(); ++B ) {                         // Loop over all bodies
    for( int d=0; d!=3; ++d )                                   //  Loop over each dimension
      B.pos()[d] = rand()/(1.+RAND_MAX)*2-1;                    //   Initialize positions
    B.scal() = 1./B.size();                                     //  Initialize source value
  }                                                             // End loop over all bodies
  toc = get_time();
  std::cout << "Initialize    : " << toc-tic << std::endl;

  tic = get_time();
  tree.setDomain();
  toc = get_time();
  std::cout << "Set domain    : " << toc-tic << std::endl;

  tic = get_time();
  int level,digit,Nbucket;
  bigint *index;
  index = new bigint [Nbody];
  tree.get_morton(index);
  level = tree.get_max_level();
  Nbucket = 1 << 3*level;
  tree.sort(index,B,Nbucket);
  digit = floor(log10(1 << 3*level)/2)+1;
  Nbucket = pow(10,digit);
//  tree.sortLarge(index,B,Nbucket);
  toc = get_time();
  std::cout << "Construct tree: " << toc-tic << std::endl;
  for( B=B.begin(); B!=B.end(); ++B ) {
    std::cout << B << " " << index[B] << std::endl;
  }
  delete[] index;
}
