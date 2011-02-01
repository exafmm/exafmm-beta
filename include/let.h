#ifndef let_h
#define let_h
#include "partition.h"
#include "tree.h"
#include "bottomup.h"
#include "topdown.h"

class LocalEssentialTree : public Partition,
                           public BottomUpTreeConstructor,
                           public TopDownTreeConstructor {
public:
  LocalEssentialTree(Bodies &b) : TreeStructure(b),
                                  Partition(b),
                                  BottomUpTreeConstructor(b),
                                  TopDownTreeConstructor(b) {}
  ~LocalEssentialTree() {}

  void commBodies() {
    int MPI_TYPE = getType(XMIN[0]);
    vect *xmin = new vect [SIZE];
    vect *xmax = new vect [SIZE];
    MPI_Allgather(&XMIN[0],3,MPI_TYPE,&xmin[0][0],3,MPI_TYPE,MPI_COMM_WORLD);
    MPI_Allgather(&XMAX[0],3,MPI_TYPE,&xmax[0][0],3,MPI_TYPE,MPI_COMM_WORLD);
    std::vector<int> srnks,scnts;
    std::vector<C_iter> scells;
/*
    int offset(0);
    for( int irank=0; irank!=SIZE; ++irank ) {
      int ic(0);
      for( int d=0; d!=3; ++d ) {
        if(xmin[irank][d] < XMAX[d] + EPS * abs(XMAX[d]) &&
           XMIN[d] < xmax[irank][d] + EPS * abs(xmax[irank][d]))
          ic++;
      }
      if( ic == 3 ) {
        C_iter C(C0);
        std::cout << RANK << std::endl;
        while( C->NCHILD == 0 ) {
          vect dist;
          for( int d=0; d!=3; ++d )
            dist[d] = (C->X[d] > xmax[irank][d])*
                      (C->X[d] - xmax[irank][d])+
                      (C->X[d] < xmin[irank][d])*
                      (C->X[d] - xmin[irank][d]);
          real R = std::sqrt(norm(dist));
          if( C0->R + C->R > THETA * R )
            scells.push_back(C);
        }
        srnks.push_back(irank);
        scnts.push_back(scells.size()-offset);
        offset = scells.size();
      }
    }
*/
    delete[] xmin;
    delete[] xmax;
  }

};
#endif
