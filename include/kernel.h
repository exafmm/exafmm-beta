#ifndef kernel_h
#define kernel_h
#include "tree.h"

class Kernels : public TreeStructure {
public:
  Kernels(Bodies &b) : TreeStructure(b) {}                      // Constructor
  ~Kernels() {}                                                 // Destructor

  void direct() {
    real invDist,invDistCube;
    vect dist;
    for( Bi=bodies.begin(); Bi!=bodies.end(); ++Bi ) {
      real pot(0);
      vect acc(0);
      for( Bj=bodies.begin(); Bj!=Bi; ++Bj ){
        dist = Bi->pos-Bj->pos;
        invDist = 1.0/std::sqrt(dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]+EPS);
        invDistCube = invDist*invDist*invDist;
        pot += invDist*Bj->scal;
        acc -= dist*invDistCube*Bj->scal;
        Bj->pot += invDist*Bi->scal;
        Bj->acc += dist*invDistCube*Bi->scal;
      }
      Bi->pot = pot;
      Bi->acc = acc;
    }
  }

  void direct_ij() {
    vect dist;
    real invDist,invDistCube;
    for( Bi=bodies.begin(); Bi!=bodies.end(); ++Bi ) {
      for( Bj=bodies.begin(); Bj!=bodies.end(); ++Bj ){
        dist = Bi->pos-Bj->pos;
        invDist = 1.0/sqrt(dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]+EPS);
        invDistCube = Bj->scal*invDist*invDist*invDist;
        Bi->pot += Bj->scal*invDist;
        Bi->acc -= dist*invDistCube;
      }
    }
  }

};

#endif
