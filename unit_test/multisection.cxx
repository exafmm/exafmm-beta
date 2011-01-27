#include "partition.h"
#include "dataset.h"
#include "bottomup.h"
#ifdef VTK
#include "vtk.h"
#endif

int main() {
  double tic,toc;
  int const numBodies(1000000);
  tic = get_time();
  Bodies bodies(numBodies);
  Bodies bodies2(numBodies);
  BottomUpTreeConstructor T(bodies);
  Dataset D(bodies);
  Partition mpi(bodies);
  toc = get_time();
  mpi.print("Allocate      : ",0);
  mpi.print(toc-tic);

  tic = get_time();
  srand(mpi.rank()+1);
  D.random();
  toc = get_time();
  mpi.print("Set bodies    : ",0);
  mpi.print(toc-tic);

  tic = get_time();
  T.setDomain();
  mpi.setDomain();
  toc = get_time();
  mpi.print("Set domain    : ",0);
  mpi.print(toc-tic);

  tic = get_time();
  T.setMorton();
  toc = get_time();
  mpi.print("Set Index     : ",0);
  mpi.print(toc-tic);

  tic = get_time();
  T.sort(T.Ibody,bodies,bodies2);
  toc = get_time();
  mpi.print("Sort Index    : ",0);
  mpi.print(toc-tic);

  tic = get_time();
  bigint nth = numBodies * mpi.size() / 2;
  nth = mpi.nth_element(T.Ibody,numBodies,nth);
  toc = get_time();
  mpi.print("Nth element   : ",0);
  mpi.print(toc-tic);
  mpi.print(nth);
  mpi.print(T.Ibody,numBodies-10,numBodies);
  for( B_iter B=bodies.begin(); B!=bodies.end(); ++B )
    T.Ibody[B-bodies.begin()] = T.Ibody[B-bodies.begin()] > nth;

  int level = int(log(mpi.size()) / M_LN2 - 1e-5) + 1;
  int nprocs[2] = {0, 0};
  int offset[2] = {0, 0};
  int  color[3] = {0, 0, 0};
  int    key[3] = {0, 0, 0};
  nprocs[0] = mpi.size();
  for( int l=0; l!=level; ++l ) {
    int oldnprocs = nprocs[0];
    int  oldcolor = color[0];
    for( int i=1; i>=0; --i ) {
      int isplit = (oldnprocs + i) / 2;
      if( mpi.rank() - offset[0] < isplit ) {
        nprocs[i] = isplit;
        offset[i] = offset[0];
         color[i] =  color[0] * 2;
      } else {
        nprocs[i] = oldnprocs - isplit;
        offset[i] = offset[0] + isplit;
         color[i] =  color[0] * 2 + 1;
      }
      key[i] = mpi.rank() - offset[i];
    }
    key[2] = color[1] % 2;
    color[2] = key[1] + oldcolor * (1 << (level - l - 1));
    MPI_Comm MPI_COMM[3];
    for( int i=0; i!=3; ++i )
      MPI_Comm_split(MPI_COMM_WORLD,color[i],key[i],&MPI_COMM[i]);
#ifdef DEBUG
    mpi.print("level : ",0);
    mpi.print(l,0);
    mpi.print("\n",0);
    mpi.print("key   : \n",0);
    mpi.print(key,0,3);
    mpi.print("color : \n",0);
    mpi.print(color,0,3);
#endif
  }

#ifdef VTK
  if( mpi.rank() == 0 ) {
    int Ncell(0);
    vtkPlot vtk;
    vtk.setDomain(T.getR0(),T.getX0());
    vtk.setGroupOfPoints(T.Ibody,bodies,Ncell);
    vtk.plot(Ncell);
  }
#endif
}
