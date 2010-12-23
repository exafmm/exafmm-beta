#include <iostream>
#include <cmath>
#include <mpi.h>

int main(int argc, char **argv) {
  int const N(32);
  int nprocs,myrank,size;
  int idata[N],isend[N],irecv[N];
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  size = N/nprocs;
  for( int i=0; i!=N; ++i ) idata[i] = i/size;
  if( myrank == 0 ) std::cout << "before" << std::endl;
  for( int irank=0; irank!=nprocs; ++irank ) {
    MPI_Barrier(MPI_COMM_WORLD);
    sleep(.1);
    if( myrank == irank ) {
      std::cout << myrank << " : ";
      for( int i=0; i!=N; ++i ) std::cout << idata[i] << " ";
      std::cout << std::endl;
    }
  }

  MPI_Comm MPI_COMM;
  int color[nprocs],key[nprocs];
  int level,npart,proc,rank,ic;
  level = log(1.+nprocs)/M_LN2;
  for( int l=0; l!=level; ++l ) {
    npart = 1 << l;
    for( int irank=0; irank!=nprocs; ++irank ) {
      color[irank] = (irank/npart/2)*npart+irank%npart;
      key[irank] = irank/npart%2;
    }
    MPI_Comm_split(MPI_COMM_WORLD,color[myrank],key[myrank],&MPI_COMM);
    MPI_Comm_size(MPI_COMM,&proc);
    MPI_Comm_rank(MPI_COMM,&rank);
    if( myrank == 0 ) {
      std::cout << "level : " << level-l << std::endl;
      std::cout << "color : " ;
      for( int irank=0; irank!=nprocs; ++irank )
        std::cout << color[irank] << " ";
      std::cout << std::endl << "key   : " ;
      for( int irank=0; irank!=nprocs; ++irank )
        std::cout << key[irank] << " ";
      std::cout << std::endl;
    }
    for( int ii=0; ii!=nprocs/2; ++ii ) {
      ic=0;
      for( int i=0; i!=nprocs; ++i ) {
        if( color[i] == ii ) {
          for( int isize=0; isize!=size; ++isize ) {
            isend[ic] = idata[i*size+isize];
            ic++;
          }
        }
      }
      MPI_Alltoall(isend,size,MPI_INT,irecv,size,MPI_INT,MPI_COMM);
      ic=0;
      for( int i=0; i!=nprocs; ++i ) {
        if( color[i] == ii ) {
          for( int isize=0; isize!=size; ++isize ) {
            idata[i*size+isize] = irecv[ic];
            ic++;
          }
        }
      }
    }
  }
  if( myrank == 0 ) std::cout << "after" << std::endl;
  for( int irank=0; irank!=nprocs; ++irank ) {
    MPI_Barrier(MPI_COMM_WORLD);
    sleep(.1);
    if( myrank == irank ) {
      std::cout << myrank << " : ";
      for( int i=0; i!=N; ++i ) std::cout << idata[i] << " ";
      std::cout << std::endl;
    }
  }
  MPI_Finalize();
}
