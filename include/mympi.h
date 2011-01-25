#ifndef mympi_h
#define mympi_h
#include <mpi.h>
#include <iostream>
#include <cmath>

class MyMPI {                                                   // My own MPI utilities
private:
  int const WAIT;                                               // Waiting time between output of different ranks
protected:
  int       SIZE;                                               // Number of MPI processes
  int       RANK;                                               // Index of current MPI process
public:
  MyMPI() : WAIT(100) {                                         // Constructor, initialize WAIT time
    int argc(0);                                                // Dummy argument count
    char **argv;                                                // Dummy argument value
    MPI_Init(&argc,&argv);                                      // Initialize MPI communicator
    MPI_Comm_size(MPI_COMM_WORLD,&SIZE);                        // Get number of MPI processes
    MPI_Comm_rank(MPI_COMM_WORLD,&RANK);                        // Get index of current MPI process
  }

  ~MyMPI() {                                                    // Destructor
    MPI_Finalize();                                             // Finalize MPI communicator
  }

  int size() { return SIZE; }                                   // Number of MPI processes
  int rank() { return RANK; }                                   // Index of current MPI process

  bool isPowerOfTwo(int const n) {                              // If n is power of two return true
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
    for( int irank=0; irank!=SIZE; ++irank ) {                  // Loop over ranks
      MPI_Barrier(MPI_COMM_WORLD);                              // Sync processes
      usleep(WAIT);                                             // Wait "WAIT" milliseconds
      if( RANK == irank ) std::cout << data << " ";             // Print "data"
    }                                                           // End loop over ranks
    MPI_Barrier(MPI_COMM_WORLD);                                // Sync processes
    usleep(WAIT);                                               // Wait "WAIT" milliseconds
    if( RANK == 0 ) std::cout << std::endl;                     // New lint
  }

  template<typename T>
  void print(T data, int const irank) {
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(WAIT);
    if( RANK == irank ) std::cout << data;
  }

  template<typename T>
  void print(T *data, int const begin, int const end) {
    for( int irank=0; irank!=SIZE; ++irank ) {
      MPI_Barrier(MPI_COMM_WORLD);
      usleep(WAIT);
      if( RANK == irank ) {
        std::cout << RANK << " : ";
        for( int i=begin; i!=end; ++i ) std::cout << data[i] << " ";
        std::cout << std::endl;
      }
    }
  }

  template<typename T>
  void print(T *data, int const begin, int const end, int const irank) {
    MPI_Barrier(MPI_COMM_WORLD);
    usleep(WAIT);
    if( RANK == irank ) {
      std::cout << RANK << " : ";
      for( int i=begin; i!=end; ++i ) std::cout << data[i] << " ";
      std::cout << std::endl;
    }
  }

  template<typename T>
  void alltoall(T *data, int const numData, int const count) {
    int type;
    T *send,*recv,sample(0);
    type = setMPItype(sample);
    send = new T [numData];
    recv = new T [numData];

    if( isPowerOfTwo(SIZE) ) {
      MPI_Comm MPI_COMM;
      int level,npart,color,key,idata,isend,irecv,chunk;
      level = int(log(1. + SIZE) / M_LN2);
      chunk = count * SIZE / 2;
      for( int l=0; l!=level; ++l ) {
        npart = 1 << l;
        color = (RANK / npart / 2) * npart + RANK % npart;
        key = RANK / npart % 2;
        MPI_Comm_split(MPI_COMM_WORLD,color,key,&MPI_COMM);
#ifdef DEBUG
        print("level : ",0);
        print(level-l,0);
        print("\n",0);
        print("color : ",0);
        print(color);
        print("key   : ",0);
        print(key);
#endif
        for( int irank=0; irank!=SIZE/2; ++irank ) {
          for( int i=0; i!=2; ++i ) {
            idata = (irank / npart) * 2 * npart + irank % npart + i * npart;
            isend = i * SIZE / 2 + irank;
            for( int icount=0; icount!=count; ++icount )
              send[isend*count+icount] = data[idata*count+icount];
          }
        }
        MPI_Alltoall(send,chunk,type,recv,chunk,type,MPI_COMM);
        for( int irank=0; irank!=SIZE/2; ++irank ) {
          for( int i=0; i!=2; ++i ) {
            idata = (irank / npart) * 2 * npart + irank % npart + i * npart;
            irecv = i * SIZE / 2 + irank;
            for( int icount=0; icount!=count; ++icount )
              data[idata*count+icount] = recv[irecv*count+icount];
          }
        }
      }
    } else {
      MPI_Alltoall(data,count,type,recv,count,type,MPI_COMM_WORLD);
      for( int i=0; i!=numData; ++i ) data[i] = recv[i];
    }
    delete[] send;
    delete[] recv;
  }

  template<typename T>
  void alltoallv(T *data, int const numData, int *scnt) {
    int type;
    int *sdsp,*rcnt,*rdsp;
    T *send,*recv,sample(0);
    type = setMPItype(sample);
    send = new T [numData];
    recv = new T [numData];
    sdsp = new int [SIZE];
    rcnt = new int [SIZE];
    rdsp = new int [SIZE];

    sdsp[0] = 0;
    for( int irank=0; irank!=SIZE-1; ++irank ) {
      sdsp[irank+1] = sdsp[irank] + scnt[irank];
    }

    if( isPowerOfTwo(SIZE) ) {
      MPI_Comm MPI_COMM;
      int level,npart,color,key,idata,isend,irecv;
      int scnt2[2],sdsp2[2],rcnt2[2],rdsp2[2];
      int scntd[SIZE],rcntd[SIZE],irev[SIZE];

      level = int(log(1. + SIZE) / M_LN2);
      for( int l=0; l!=level; ++l ) {
        npart = 1 << l;
        color = (RANK / npart / 2) * npart + RANK % npart;
        key = RANK / npart % 2;
        MPI_Comm_split(MPI_COMM_WORLD,color,key,&MPI_COMM);
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
          for( int irank=0; irank!=SIZE/2; ++irank ) {
            idata = (irank / npart) * 2 * npart + irank % npart + i * npart;
            isend = i * SIZE / 2 + irank;
            irev[idata] = isend;
            scntd[isend] = scnt[idata];
            scnt2[i] += scnt[idata];
            for( int id=sdsp[idata]; id!=sdsp[idata]+scnt[idata]; ++id,++ic )
              send[ic] = data[id];
          }
        }
        MPI_Alltoall(scntd,SIZE/2,MPI_INT,rcnt,SIZE/2,MPI_INT,MPI_COMM);
        MPI_Alltoall(scnt2,1,MPI_INT,rcnt2,1,MPI_INT,MPI_COMM);
        sdsp2[0] = 0; sdsp2[1] = scnt2[0];
        rdsp2[0] = 0; rdsp2[1] = rcnt2[0];
        MPI_Alltoallv(send,scnt2,sdsp2,type,recv,rcnt2,rdsp2,type,MPI_COMM);
        rdsp[0] = 0;
        for( int irank=0; irank!=SIZE-1; ++irank )
          rdsp[irank+1] = rdsp[irank] + rcnt[irank];
        ic = 0;
        for( int i=0; i!=2; ++i ) {
          for( int irank=0; irank!=SIZE/2; ++irank ) {
            idata = (irank / npart) * 2 * npart + irank % npart + i * npart;
            irecv = i * SIZE / 2 + irank;
            rcntd[idata] = rcnt[irecv];
            idata = irev[irecv];
            for( int id=rdsp[idata]; id!=rdsp[idata]+rcnt[idata]; ++id,++ic )
              data[ic] = recv[id];
          }
        }
        rdsp[0] = 0;
        for( int irank=0; irank!=SIZE-1; ++irank )
          rdsp[irank+1] = rdsp[irank] + rcntd[irank];
        for( int irank=0; irank!=SIZE; ++irank ) {
          scnt[irank] = rcntd[irank];
          sdsp[irank] = rdsp[irank];
        }
#ifdef DEBUG
        print("rcnt\n",0);
        print(rcnt,0,SIZE);
        print("rcnt\n",0);
        print(rcntd,0,SIZE);
        print("send\n",0);
        print(send,0,numData);
        print("recv\n",0);
        print(recv,0,numData);
        print("data\n",0);
        print(data,0,numData);
#endif
      }
    } else {
      MPI_Alltoall(scnt,1,MPI_INT,rcnt,1,MPI_INT,MPI_COMM_WORLD);
      rdsp[0] = 0;
      for( int irank=0; irank!=SIZE-1; ++irank )
        rdsp[irank+1] = rdsp[irank] + rcnt[irank];

      MPI_Alltoallv(data,scnt,sdsp,type,recv,rcnt,rdsp,type,MPI_COMM_WORLD);
      for( int i=0; i!=numData; ++i ) data[i] = recv[i];
    }
    delete[] send;
    delete[] recv;
    delete[] sdsp;
    delete[] rcnt;
    delete[] rdsp;
  }

};

#endif
