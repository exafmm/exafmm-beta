#ifndef biotsavart_h
#define biotsavart_h

void Kernel::BiotSavartP2P_CPU() {
  for( B_iter BI=BI0; BI!=BIN; ++BI ) {
    for( B_iter BJ=BJ0; BJ!=BJN; ++BJ ) {
      vect dist = BI->X - BJ->X - Xperiodic;
      real S2 = 2 * BJ->S * BJ->S;
      real R2  = norm(dist) + EPS2;
      real RS = R2 / S2;
      real cutoff = 0.25 / M_PI / R2 / std::sqrt(R2) * (erf( std::sqrt(RS) )
                  - std::sqrt(4 / M_PI * RS) * exp(-RS));
      BI->vel[0] += (dist[1] * BJ->Q[2] - dist[2] * BJ->Q[1]) * cutoff;
      BI->vel[1] += (dist[2] * BJ->Q[0] - dist[0] * BJ->Q[2]) * cutoff;
      BI->vel[2] += (dist[0] * BJ->Q[1] - dist[1] * BJ->Q[0]) * cutoff;
    }
  }
}

#endif
