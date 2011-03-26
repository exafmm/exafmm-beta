#ifndef stretching_h
#define stretching_h

void Kernel::StretchingP2P_CPU() {
  for( B_iter BI=BI0; BI!=BIN; ++BI ) {
    for( B_iter BJ=BJ0; BJ!=BJN; ++BJ ) {
      vect dist = BI->X - BJ->X - Xperiodic;
      real S2 = 2 * BJ->S * BJ->S;
      real R2  = norm(dist) + EPS2;
      real RS = R2 / S2;
      real cutoff = 0.25 / M_PI / R2 / std::sqrt(R2) * (erf( std::sqrt(RS) )
                  - std::sqrt(4 / M_PI * RS) * exp(-RS));
      BI->dQdt[0] += (BI->Q[1] * BJ->Q[2] - BI->Q[2] * BJ->Q[1]) * cutoff;
      BI->dQdt[1] += (BI->Q[2] * BJ->Q[0] - BI->Q[0] * BJ->Q[2]) * cutoff;
      BI->dQdt[2] += (BI->Q[0] * BJ->Q[1] - BI->Q[1] * BJ->Q[0]) * cutoff;
      cutoff = 0.25 / M_PI / R2 / R2 / std::sqrt(R2) * (3 * erf( std::sqrt(RS) )
             - (2 * RS + 3) * std::sqrt(4 / M_PI * RS) * exp(-RS))
             * (BI->Q[0] * dist[0] + BI->Q[1] * dist[1] + BI->Q[2] * dist[2]);
      BI->dQdt[0] += (BJ->Q[1] * dist[2] - BJ->Q[2] * dist[1]) * cutoff;
      BI->dQdt[1] += (BJ->Q[2] * dist[0] - BJ->Q[0] * dist[2]) * cutoff;
      BI->dQdt[2] += (BJ->Q[0] * dist[1] - BJ->Q[1] * dist[0]) * cutoff;
    }
  }
}

#endif
