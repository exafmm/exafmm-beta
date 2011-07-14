#ifndef laplace_h
#define laplace_h

void Kernel::LaplaceP2P_CPU() {
  for( B_iter BI=BI0; BI!=BIN; ++BI ) {
    for( B_iter BJ=BJ0; BJ!=BJN; ++BJ ) {
      vect dist = BI->X - BJ->X - Xperiodic;
      real R2 = norm(dist);
//      real invR = 1 / std::sqrt(norm(dist) + EPS2);
      real invR = 1 / std::sqrt(R2);
      if( R2 == 0 ) invR = 0;
      real invR3 = BJ->SRC[0] * invR * invR * invR;
      BI->TRG[0] += BJ->SRC[0] * invR;
      BI->TRG[1] -= dist[0] * invR3;
      BI->TRG[2] -= dist[1] * invR3;
      BI->TRG[3] -= dist[2] * invR3;
    }
  }
}

#endif
