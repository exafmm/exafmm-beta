#ifndef gaussian_h
#define gaussian_h

void Kernel::GaussianP2P_CPU() {
  for( B_iter BI=BI0; BI!=BIN; ++BI ) {
    for( B_iter BJ=BJ0; BJ!=BJN; ++BJ ) {
      vect dist = BI->X - BJ->X - Xperiodic;
      real S2 = 2 * BJ->SRC[3] * BJ->SRC[3];
      real R2  = norm(dist) + EPS2;
      BI->TRG[0] += BJ->SRC[0] / (M_PI * S2) / std::sqrt(M_PI * S2) * exp(-R2 / S2);
    }
  }
}

#endif
