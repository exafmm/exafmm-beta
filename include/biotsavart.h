#ifndef biotsavart_h
#define biotsavart_h

template<>
void Kernel<BiotSavart>::P2P_CPU() {                            // Biot-Savart P2P kernel on CPU
  for( B_iter BI=BI0; BI!=BIN; ++BI ) {                         // Loop over target bodies
    for( B_iter BJ=BJ0; BJ!=BJN; ++BJ ) {                       //  Loop over source bodies
      vect dist = BI->X - BJ->X - Xperiodic;                    //   Distance vector from source to target
      real S2 = 2 * BJ->SRC[3] * BJ->SRC[3];                    //    2 * sigma^2
      real R2  = norm(dist) + EPS2;                             //    R^2 + epsilon^2
      real RS = R2 / S2;                                        //    R^2 / (2 * simga^2)
      real cutoff = 0.25 / M_PI / R2 / std::sqrt(R2) * (erf( std::sqrt(RS) )// cutoff function
                  - std::sqrt(4 / M_PI * RS) * exp(-RS));
      BI->TRG[0] += (dist[1] * BJ->SRC[2] - dist[2] * BJ->SRC[1]) * cutoff;// x component of curl G * cutoff
      BI->TRG[1] += (dist[2] * BJ->SRC[0] - dist[0] * BJ->SRC[2]) * cutoff;// y component of curl G * cutoff
      BI->TRG[2] += (dist[0] * BJ->SRC[1] - dist[1] * BJ->SRC[0]) * cutoff;// z component of curl G * cutoff
    }                                                           //  End loop over source bodies
  }                                                             // End loop over target bodies
}

#endif
