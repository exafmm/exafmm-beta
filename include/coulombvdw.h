#ifndef coulombvdw_h
#define coulombvdw_h

template<>
void Kernel<CoulombVdW>::P2P_CPU() {                            // Coulomb + Van der Waals P2P kernel on CPU
  for( B_iter BI=BI0; BI!=BIN; ++BI ) {                         // Loop over target bodies
    int atypei = BI->SRC[1];                                    //  Atom type of target
    for( B_iter BJ=BJ0; BJ!=BJN; ++BJ ) {                       //  Loop over source bodies
      int atypej = BJ->SRC[1];                                  //   Atom type of source
      vect dist = BI->X - BJ->X - Xperiodic;                    //   Distance vector from source to target
      real R2 = norm(dist);                                     //   R squared
      if( R2 != 0 ) {                                           //   Exclude self interaction
        real rs = RSCALE[atypei*ATOMS+atypej];                  //    r scale
        real gs = GSCALE[atypei*ATOMS+atypej];                  //    g scale
        real R2s = R2 * rs;                                     //    R^2 * r scale
        real invR2 = 1.0 / R2s;                                 //    1 / R^2
        real invR6 = invR2 * invR2 * invR2;                     //    1 / R^6
        real dtmp = gs * invR6 * invR2 * (2.0 * invR6 - 1.0);   //    g scale / R^2 * (2 / R^12 + 1 / R^6)
        BI->TRG[1] += dist[0] * dtmp;                           //    x component of Van der Waals force
        BI->TRG[2] += dist[1] * dtmp;                           //    y component of Van der Waals force
        BI->TRG[3] += dist[2] * dtmp;                           //    z component of Van der Waals force
      }                                                         //   End if for self interaction 
    }                                                           //  End loop over source bodies
  }                                                             // End loop over target bodies
}

#endif
