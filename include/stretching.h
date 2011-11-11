#ifndef stretching_h
#define stretching_h

void Kernel::StretchingP2P_CPU() {                              // Stretching P2P kernel on CPU
  for( B_iter BI=BI0; BI!=BIN; ++BI ) {                         // Loop over target bodies
    for( B_iter BJ=BJ0; BJ!=BJN; ++BJ ) {                       //  Loop over source bodies
      vect dist = BI->X - BJ->X - Xperiodic;                    //   Distance vector from source to target
      real S2 = 2 * BJ->SRC[3] * BJ->SRC[3];                    //   2 * simga^2
      real R2  = norm(dist) + EPS2;                             //   R^2 + epsilon^2
      real RS = R2 / S2;                                        //   R^2 / (2 * sigma^2)
      real cutoff = 0.25 / M_PI / R2 / std::sqrt(R2) * (erf( std::sqrt(RS) )// cutoff function for first term
                  - std::sqrt(4 / M_PI * RS) * exp(-RS));
      BI->TRG[0] += (BI->SRC[1] * BJ->SRC[2] - BI->SRC[2] * BJ->SRC[1]) * cutoff;// x component of first term
      BI->TRG[1] += (BI->SRC[2] * BJ->SRC[0] - BI->SRC[0] * BJ->SRC[2]) * cutoff;// y component of first term
      BI->TRG[2] += (BI->SRC[0] * BJ->SRC[1] - BI->SRC[1] * BJ->SRC[0]) * cutoff;// z component of first term
      cutoff = 0.25 / M_PI / R2 / R2 / std::sqrt(R2) * (3 * erf( std::sqrt(RS) )// cutoff function for second term
             - (2 * RS + 3) * std::sqrt(4 / M_PI * RS) * exp(-RS))
             * (BI->SRC[0] * dist[0] + BI->SRC[1] * dist[1] + BI->SRC[2] * dist[2]);
      BI->TRG[0] += (BJ->SRC[1] * dist[2] - BJ->SRC[2] * dist[1]) * cutoff;// x component of second term
      BI->TRG[1] += (BJ->SRC[2] * dist[0] - BJ->SRC[0] * dist[2]) * cutoff;// y component of second term
      BI->TRG[2] += (BJ->SRC[0] * dist[1] - BJ->SRC[1] * dist[0]) * cutoff;// z component of second term
    }                                                           //  End loop over source bodies
  }                                                             // End loop over target bodies
}

#endif
