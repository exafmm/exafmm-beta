#ifndef stretching_h
#define stretching_h

void Kernel::StretchingP2P_CPU() {
  for( B_iter BI=BI0; BI!=BIN; ++BI ) {
    for( B_iter BJ=BJ0; BJ!=BJN; ++BJ ) {
      vect dist = BI->X - BJ->X - Xperiodic;
      real S2 = 2 * BJ->SRC[3] * BJ->SRC[3];
      real R2  = norm(dist) + EPS2;
      real RS = R2 / S2;
      real cutoff = 0.25 / M_PI / R2 / std::sqrt(R2) * (erf( std::sqrt(RS) )
                  - std::sqrt(4 / M_PI * RS) * exp(-RS));
      BI->TRG[0] += (BI->SRC[1] * BJ->SRC[2] - BI->SRC[2] * BJ->SRC[1]) * cutoff;
      BI->TRG[1] += (BI->SRC[2] * BJ->SRC[0] - BI->SRC[0] * BJ->SRC[2]) * cutoff;
      BI->TRG[2] += (BI->SRC[0] * BJ->SRC[1] - BI->SRC[1] * BJ->SRC[0]) * cutoff;
      cutoff = 0.25 / M_PI / R2 / R2 / std::sqrt(R2) * (3 * erf( std::sqrt(RS) )
             - (2 * RS + 3) * std::sqrt(4 / M_PI * RS) * exp(-RS))
             * (BI->SRC[0] * dist[0] + BI->SRC[1] * dist[1] + BI->SRC[2] * dist[2]);
      BI->TRG[0] += (BJ->SRC[1] * dist[2] - BJ->SRC[2] * dist[1]) * cutoff;
      BI->TRG[1] += (BJ->SRC[2] * dist[0] - BJ->SRC[0] * dist[2]) * cutoff;
      BI->TRG[2] += (BJ->SRC[0] * dist[1] - BJ->SRC[1] * dist[0]) * cutoff;
    }
  }
}

#endif
