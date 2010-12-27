#ifndef mympi_h
#define mympi_h
#include <mpi.h>
#include <iostream>
#include <cmath>

class MyMPI {                                                   // My own MPI utilities
  int const WAIT;                                               // Waiting time between output of different ranks
public:
  MyMPI() : WAIT(100) {}                                        // Constructor, initialize WAIT time
  ~MyMPI() {}                                                   // Destructor

  bool isPowerOfTwo(int const n) {                              // If nprocs is power of two
    return ((n != 0) && !(n & (n - 1)));                        // Decrement and compare bits
  }

  template<typename T>
  int setMPItype(T sample) {
    int type;
    switch ( sizeof(sample) ) {
      case  1 : type = MPI_CHAR; break;
      case  2 : type = MPI_SHORT; break;
      case  4 : type = MPI_INT; break;
      case  8 : type = MPI_DOUBLE; break;
      case 16 : type = MPI_LONG_DOUBLE; break;
    }
    return type;
  }

  template<typename T>                                          // Auto-detect data type to print
  void print(T data) {                                          // Print a scalar value on all ranks
    int nprocs,myrank;                                          // Define number of processes and process ID
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);                      // Get number of processes
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);                      // Get process ID
    for( int irank=0; irank!=nprocs; ++irank ) {                // Loop over ranks
      MPI_Barrier(MPI_COMM_WORLD);                              // Sync processes
      usleep(WAIT);                                             // Wait "WAIT" milliseconds
      if( myrank == irank ) std::cout << data << " ";           // Print "data"
    }                                                           // End loop over ranks
    MPI_Barrier(MPI_COMM_WORLD);                                // Sync processes
    usleep(WAIT);                                               // Wait "WAIT" milliseconds
    if( myrank == 0 ) std::cout << std::endl;                   // New lint
  }

  template<typename T>
  void print(T data, int const irank) {
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(WAIT);
    if( myrank == irank ) std::cout << data;
  }

  template<typename T>
  void print(T *data, int const ista, int const iend) {
    int nprocs,myrank;
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    for( int irank=0; irank!=nprocs; ++irank ) {
      MPI_Barrier(MPI_COMM_WORLD);
      usleep(WAIT);
      if( myrank == irank ) {
        std::cout << myrank << " : ";
        for( int i=ista; i!=iend; ++i ) std::cout << data[i] << " ";
        std::cout << std::endl;
      }
    }
  }

  template<typename T>
  void print(T *data, int const ista, int const iend, int const irank) {
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(WAIT);
    if( myrank == irank ) {
      std::cout << myrank << " : ";
      for( int i=ista; i!=iend; ++i ) std::cout << data[i] << " ";
      std::cout << std::endl;
    }
  }

  template<typename T>
  void alltoall(T *data, int const N, int const size) {
    int nprocs,myrank,type;
    T *send,*recv,sample(0);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    type = setMPItype(sample);
    send = new T [N];
    recv = new T [N];

    if( isPowerOfTwo(nprocs) ) {
      MPI_Comm MPI_COMM;
      int level,npart,color,key,proc,rank,idata,isend,irecv,chunk;
      level = log(1.+nprocs)/M_LN2;
      chunk = size*nprocs/2;
      for( int l=0; l!=level; ++l ) {
        npart = 1 << l;
        color = (myrank/npart/2)*npart+myrank%npart;
        key = myrank/npart%2;
        MPI_Comm_split(MPI_COMM_WORLD,color,key,&MPI_COMM);
        MPI_Comm_size(MPI_COMM,&proc);
        MPI_Comm_rank(MPI_COMM,&rank);
#ifdef DEBUG
        print("level : ",0);
        print(level-l,0);
        print("\n",0);
        print("color : ",0);
        print(color);
        print("key   : ",0);
        print(key);
#endif
        for( int irank=0; irank!=nprocs/2; ++irank ) {
          for( int i=0; i!=2; ++i ) {
            idata = (irank/npart)*2*npart+irank%npart+i*npart;
            isend = i*nprocs/2+irank;
            for( int isize=0; isize!=size; ++isize ) {
              send[isend*size+isize] = data[idata*size+isize];
            }
          }
        }
        MPI_Alltoall(send,chunk,type,recv,chunk,type,MPI_COMM);
        for( int irank=0; irank!=nprocs/2; ++irank ) {
          for( int i=0; i!=2; ++i ) {
            idata = (irank/npart)*2*npart+irank%npart+i*npart;
            irecv = i*nprocs/2+irank;
            for( int isize=0; isize!=size; ++isize ) {
              data[idata*size+isize] = recv[irecv*size+isize];
            }
          }
        }
      }
    } else {
      MPI_Alltoall(data,size,type,recv,size,type,MPI_COMM_WORLD);
      for( int i=0; i!=N; ++i ) data[i] = recv[i];
    }
    delete[] send;
    delete[] recv;
  }

  template<typename T>
  void alltoallv(T *data, int const N, int *scnt) {
    int nprocs,myrank,type;
    int *sdsp,*rcnt,*rdsp;
    T *send,*recv,sample(0);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
    type = setMPItype(sample);
    send = new T [N];
    recv = new T [N];
    sdsp = new int [nprocs];
    rcnt = new int [nprocs];
    rdsp = new int [nprocs];

    sdsp[0] = 0;
    for( int irank=0; irank!=nprocs-1; ++irank ) {
      sdsp[irank+1] = sdsp[irank]+scnt[irank];
    }

    if( isPowerOfTwo(nprocs) ) {
      MPI_Comm MPI_COMM;
      int level,npart,color,key,proc,rank,idata,isend,irecv;
      int scnt2[2],sdsp2[2],rcnt2[2],rdsp2[2];
      int scntd[nprocs],rcntd[nprocs],irev[nprocs];

      level = log(1.+nprocs)/M_LN2;
      for( int l=0; l!=level; ++l ) {
        npart = 1 << l;
        color = (myrank/npart/2)*npart+myrank%npart;
        key = myrank/npart%2;
        MPI_Comm_split(MPI_COMM_WORLD,color,key,&MPI_COMM);
        MPI_Comm_size(MPI_COMM,&proc);
        MPI_Comm_rank(MPI_COMM,&rank);
#ifdef DEBUG
        print("level : ",0);
        print(level-l,0);
        print("\n",0);
        print("color : ",0);
        print(color);
        print("key   : ",0);
        print(key);
#endif
        int ic = 0;
        for( int i=0; i!=2; ++i ) {
          scnt2[i] = 0;
          for( int irank=0; irank!=nprocs/2; ++irank ) {
            idata = (irank/npart)*2*npart+irank%npart+i*npart;
            isend = i*nprocs/2+irank;
            irev[idata] = isend;
            scntd[isend] = scnt[idata];
            scnt2[i] += scnt[idata];
            for( int id=sdsp[idata]; id!=sdsp[idata]+scnt[idata]; ++id ) {
              send[ic] = data[id];
              ++ic;
            }
          }
        }
        MPI_Alltoall(scntd,nprocs/2,MPI_INT,rcnt,nprocs/2,MPI_INT,MPI_COMM);
        MPI_Alltoall(scnt2,1,MPI_INT,rcnt2,1,MPI_INT,MPI_COMM);
        sdsp2[0] = 0; sdsp2[1] = scnt2[0];
        rdsp2[0] = 0; rdsp2[1] = rcnt2[0];
        MPI_Alltoallv(send,scnt2,sdsp2,type,recv,rcnt2,rdsp2,type,MPI_COMM);
        rdsp[0] = 0;
        for( int irank=0; irank!=nprocs-1; ++irank ) {
          rdsp[irank+1] = rdsp[irank]+rcnt[irank];
        }
        ic = 0;
        for( int i=0; i!=2; ++i ) {
          for( int irank=0; irank!=nprocs/2; ++irank ) {
            idata = (irank/npart)*2*npart+irank%npart+i*npart;
            irecv = i*nprocs/2+irank;
            rcntd[idata] = rcnt[irecv];
            idata = irev[irecv];
            for( int id=rdsp[idata]; id!=rdsp[idata]+rcnt[idata]; ++id ) {
              data[ic] = recv[id];
              ++ic;
            }
          }
        }
        rdsp[0] = 0;
        for( int irank=0; irank!=nprocs-1; ++irank ) {
          rdsp[irank+1] = rdsp[irank]+rcntd[irank];
        }
        for( int irank=0; irank!=nprocs; ++irank ) {
          scnt[irank] = rcntd[irank];
          sdsp[irank] = rdsp[irank];
        }
#ifdef DEBUG
        print("rcnt\n",0);
        print(rcnt,0,nprocs);
        print("rcnt\n",0);
        print(rcntd,0,nprocs);
        print("send\n",0);
        print(send,0,N);
        print("recv\n",0);
        print(recv,0,N);
        print("data\n",0);
        print(data,0,N);
#endif
      }
    } else {
      MPI_Alltoall(scnt,1,MPI_INT,rcnt,1,MPI_INT,MPI_COMM_WORLD);
      rdsp[0] = 0;
      for( int irank=0; irank!=nprocs-1; ++irank ) {
        rdsp[irank+1] = rdsp[irank]+rcnt[irank];
      }

      MPI_Alltoallv(data,scnt,sdsp,type,recv,rcnt,rdsp,type,MPI_COMM_WORLD);
      for( int i=0; i!=N; ++i ) data[i] = recv[i];
    }
    delete[] send;
    delete[] recv;
    delete[] sdsp;
    delete[] rcnt;
    delete[] rdsp;
  }
};

#endif
