#ifndef laplace_h
#define laplace_h

void Kernel::LaplaceP2P_CPU() {                                 // Laplace P2P kernel on CPU
  for( B_iter BI=BI0; BI!=BIN; ++BI ) {                         // Loop over target bodies
    for( B_iter BJ=BJ0; BJ!=BJN; ++BJ ) {                       //  Loop over source bodies
      vect dist = BI->X - BJ->X - Xperiodic;                    //   Distance vector from source to target
//      real R2 = norm(dist);                                     //   R^2
      real invR = 1 / std::sqrt(norm(dist) + EPS2);             //   1 / sqrt(R^2 + eps^2)
//      real invR = 1 / std::sqrt(R2);                            //   1 / R
//      if( R2 == 0 ) invR = 0;                                   //   Exclude self interaction
      real invR3 = BJ->SRC[0] * invR * invR * invR;             //   charge / R^3
      BI->TRG[0] += BJ->SRC[0] * invR;                          //   potential
      BI->TRG[1] -= dist[0] * invR3;                            //   x component of force
      BI->TRG[2] -= dist[1] * invR3;                            //   y component of force
      BI->TRG[3] -= dist[2] * invR3;                            //   z component of force
    }                                                           //  End loop over source bodies
  }                                                             // End loop over target bodies
}

#endif
