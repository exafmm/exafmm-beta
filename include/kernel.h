#ifndef kernel_h
#define kernel_h

class Kernels {
public:
  void direct(B_iter B0, B_iter BN) {
    real invDist,invDistCube;
    vect dist;
    for( B_iter Bi=B0; Bi!=BN; ++Bi ) {
      real pot(0);
      vect acc(0);
      for( B_iter Bj=B0; Bj!=Bi; ++Bj ){
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

  void direct(B_iter Bi0, B_iter BiN, B_iter Bj0, B_iter BjN) {
    vect dist;
    real invDist,invDistCube;
    for( B_iter Bi=Bi0; Bi!=BiN; ++Bi ) {
      for( B_iter Bj=Bj0; Bj!=BjN; ++Bj ){
        dist = Bi->pos-Bj->pos;
        invDist = 1.0/std::sqrt(dist[0]*dist[0]+dist[1]*dist[1]+dist[2]*dist[2]+EPS);
        invDistCube = Bj->scal*invDist*invDist*invDist;
        Bi->pot += Bj->scal*invDist;
        Bi->acc -= dist*invDistCube;
      }
    }
  }

};

#endif
