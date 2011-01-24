#include "mympi.h"

int main() {
  int const numData(16);
  int count;
  int *data;
  MyMPI mpi;

  data = new int [numData];
  count = numData / mpi.size();
  for( int i=0; i!=numData; ++i ) data[i] = i / count;
  mpi.print("(before)\n",0);
  mpi.print(data,0,numData);

  mpi.alltoall(data,numData,count);

  mpi.print("(after)\n",0);
  mpi.print(data,0,numData);
  delete[] data;
}
