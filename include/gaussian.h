#ifndef gaussian_h
#define gaussian_h

void Kernel::GaussianP2P_CPU() {
  for( B_iter BI=BI0; BI!=BIN; ++BI ) {
    for( B_iter BJ=BJ0; BJ!=BJN; ++BJ ) {
      vect dist = BI->X - BJ->X - Xperiodic;
      real S2 = 2 * BJ->S * BJ->S;
      real R2  = norm(dist) + EPS2;
      BI->val += BJ->Q / (M_PI * S2) / std::sqrt(M_PI * S2) * exp(-R2 / S2);
    }
  }
}

#endif
