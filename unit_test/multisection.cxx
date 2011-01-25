#include "mympi.h"

void segScan(int *send, int *recv) {
  if( send[0] == recv[0] )
    recv[1] += send[1];
  recv[0] = send[0];
}

int main() {
  MyMPI mpi;
  int level = int(log(1. + mpi.size()) / M_LN2);
  for( int l=0; l!=level; ++l ) {
    int npart = 1 << l;
    int color = (mpi.rank() / npart / 2) * npart + mpi.rank() % npart;
    int key = mpi.rank() / npart % 2;
    int scan;

    MPI_Datatype MPI_INT2;
    MPI_Op       MPI_SEG_SCAN;
    int send[2],recv[2];
    send[0] = key;
    send[1] = 1;
    MPI_Type_vector(1,2,2,MPI_INT,&MPI_INT2);
    MPI_Type_commit(&MPI_INT2);
    MPI_Op_create((MPI_User_function *)segScan,0,&MPI_SEG_SCAN);
    MPI_Scan(&send,&recv,1,MPI_INT2,MPI_SEG_SCAN,MPI_COMM_WORLD);

    MPI_Scan(&key,&scan,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    mpi.print("level : ",0);
    mpi.print(level-l,0);
    mpi.print("\n",0);
    mpi.print("key   : ",0);
    mpi.print(key);
    mpi.print("color : ",0);
    mpi.print(color);
    mpi.print("scan  : ",0);
    mpi.print(recv[1]-1);
  }
}
