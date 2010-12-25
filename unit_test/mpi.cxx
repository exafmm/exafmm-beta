#include <iostream>
#include <mpi.h>

int main(int argc, char **argv) {
  char hostname[256];
  int nprocs,myrank;
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  gethostname(hostname,sizeof(hostname));
  for( int irank=0; irank!=nprocs; ++irank ) {
    MPI_Barrier(MPI_COMM_WORLD);
    if( myrank == irank )
      std::cout << hostname << " " << myrank << " / " << nprocs << std::endl;
    usleep(100);
  }
  MPI_Finalize();
}
