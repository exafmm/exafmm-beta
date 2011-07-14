#ifndef coulombvdw_h
#define coulombvdw_h

void Kernel::CoulombVdWP2P_CPU() {
  for( B_iter BI=BI0; BI!=BIN; ++BI ) {
    int atypei = BI->SRC[1];
    for( B_iter BJ=BJ0; BJ!=BJN; ++BJ ) {
      int atypej = BJ->SRC[1];
      vect dist = BI->X - BJ->X - Xperiodic;
      real R2 = norm(dist);
      if( R2 != 0 ) {
        real rs = RSCALE[atypei*ATOMS+atypej];
        real gs = GSCALE[atypei*ATOMS+atypej];
        real R2s = R2 * rs;
        real invR2 = 1.0 / R2s;
        real invR6 = invR2 * invR2 * invR2;
        real dtmp = gs * invR6 * invR2 * (2.0 * invR6 - 1.0);
        BI->TRG[1] += dist[0] * dtmp;
        BI->TRG[2] += dist[1] * dtmp;
        BI->TRG[3] += dist[2] * dtmp;
      }
    }
  }
}

#endif
