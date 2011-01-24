#include "mympi.h"

int main() {
  int const numData(16);
  int *data,*scnt;
  MyMPI mpi;

  data = new int [numData];
  scnt = new int [mpi.size()];
  for( int i=0; i!=numData; ++i ) data[i] = i;
  for( int irank=0; irank!=mpi.size(); ++irank ) scnt[irank] = irank + 1;
  mpi.print("(before)\n",0);
  mpi.print(data,0,numData);

  mpi.alltoallv(data,numData,scnt);

  mpi.print("(after)\n",0);
  mpi.print(data,0,(mpi.rank() + 1) * mpi.size());
  delete[] data;
  delete[] scnt;
}
