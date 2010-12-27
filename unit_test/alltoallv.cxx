#include "mympi.h"

int main(int argc, char **argv) {
  int const N(16);
  int nprocs,myrank;
  int *data,*scnt;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  MyMPI mpi;

  data = new int [N];
  scnt = new int [nprocs];
  for( int i=0; i!=N; ++i ) data[i] = i;
  for( int irank=0; irank!=nprocs; ++irank ) scnt[irank] = irank+1;
  mpi.print("(before)\n",0);
  mpi.print(data,0,N);

  mpi.alltoallv(data,N,scnt);

  mpi.print("(after)\n",0);
  mpi.print(data,0,(myrank+1)*nprocs);
  delete[] data;
  delete[] scnt;
  MPI_Finalize();
}
