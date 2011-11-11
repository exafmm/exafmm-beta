#ifndef gaussian_h
#define gaussian_h

void Kernel::GaussianP2P_CPU() {                                // Gaussian P2P kernel on CPU
  for( B_iter BI=BI0; BI!=BIN; ++BI ) {                         // Loop over target bodies
    for( B_iter BJ=BJ0; BJ!=BJN; ++BJ ) {                       //  Loop over source bodies
      vect dist = BI->X - BJ->X - Xperiodic;                    //   Distance vector from source to target
      real S2 = 2 * BJ->SRC[3] * BJ->SRC[3];                    //   2 * sigma^2
      real R2  = norm(dist) + EPS2;                             //   R^2 + epsilon^2
      BI->TRG[0] += BJ->SRC[0] / (M_PI * S2) / std::sqrt(M_PI * S2) * exp(-R2 / S2);// Gaussian function
    }                                                           //  End loop over source bodies
  }                                                             // End loop over target bodies
}

#endif
