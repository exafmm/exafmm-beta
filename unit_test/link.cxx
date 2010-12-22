#include <iostream>
#include "gettime.h"

struct node {
  int I;
  node *NEXT;
  double shit[50];
};

int main() {
  int const size = 10000000;
  double tic,toc;
  node *N0,*NN,*N;
  NN = NULL;
  N = N0 = new node [size];
  for( int i=0; i!=size; ++i,++N ) {
    N->I = i;
    N->NEXT = NN;
    NN = N;
  }
  tic = get_time();
  for( int i=0; i!=size; ++i ) {
    for( int j=0; j!=10; ++j )
      N0[i].I = 0;
//    std::cout << i << std::endl;
  }
  toc = get_time();
  std::cout << "iloop: " << toc-tic << std::endl;
  tic = get_time();
  for( N=N0; N!=NN; ++N ) {
    for( int j=0; j!=10; ++j )
      N->I = 1;
//    std::cout << N->I << std::endl;
  }
  toc = get_time();
  std::cout << "ploop: " << toc-tic << std::endl;
  tic = get_time();
  for ( N=NN; N!=NULL; N=N->NEXT ) {
    for( int j=0; j!=10; ++j )
      N->I = 2;
//    std::cout << N->I << std::endl;
  }
  toc = get_time();
  delete[] N0;
  std::cout << "plink: " << toc-tic << std::endl;
}
