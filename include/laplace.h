#ifndef laplace_h
#define laplace_h

void Kernel::LaplaceP2P_CPU() {
  for( B_iter BI=BI0; BI!=BIN; ++BI ) {
    for( B_iter BJ=BJ0; BJ!=BJN; ++BJ ) {
      vect dist = BI->X - BJ->X - Xperiodic;
      real invR = 1 / std::sqrt(norm(dist) + EPS2);
      real invR3 = BJ->Q * invR * invR * invR;
      BI->pot += BJ->Q * invR;
      BI->acc -= dist * invR3;
    }
  }
}

#endif
