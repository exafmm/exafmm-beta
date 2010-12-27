#include "mympi.h"

int main(int argc, char **argv) {
  int const N(16);
  int nprocs,myrank,size;
  int *data;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  MyMPI mpi;

  data = new int [N];
  size = N/nprocs;
  for( int i=0; i!=N; ++i ) data[i] = i/size;
  mpi.print("(before)\n",0);
  mpi.print(data,0,N);

  mpi.alltoall(data,N,size);

  mpi.print("(after)\n",0);
  mpi.print(data,0,N);
  delete[] data;
  MPI_Finalize();
}
