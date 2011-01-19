#ifndef kernel_h
#define kernel_h

class Kernel {
public:
  void direct(B_iter B0, B_iter BN) {
    for( B_iter Bi=B0; Bi!=BN; ++Bi ) {
      real pot(-Bi->scal/sqrtf(EPS2));
      for( B_iter Bj=B0; Bj!=BN; ++Bj ) {
        vect dist = Bi->pos - Bj->pos;
        real r = std::sqrt(norm(dist)+EPS2);
        pot += Bj->scal / r;
      }
      Bi->pot = pot;
    }
  }

  void direct(B_iter Bi0, B_iter BiN, B_iter Bj0, B_iter BjN) {
    for( B_iter Bi=Bi0; Bi!=BiN; ++Bi ) {
      for( B_iter Bj=Bj0; Bj!=BjN; ++Bj ) {
        vect dist = Bi->pos - Bj->pos;
        real r = std::sqrt(norm(dist)+EPS2);
        Bi->pot += Bj->scal / r;
      }
    }
  }

  void P2M(C_iter C) {                                          // Evaluate P2M kernel
    C->M = 0;                                                   // Initialize Multipole expansions
    for( B_iter B=C->LEAF; B!=C->LEAF+C->NLEAF; ++B ) {         // Loop over all leafs in cell
      vect dist = C->X - B->pos;                                //  Distance between bodist[1] and cell center
      C->M[0] += B->scal;                                       //  Calculate Multipole expansions
      C->M[1] += B->scal * dist[0];
      C->M[2] += B->scal * dist[1];
      C->M[3] += B->scal * dist[2];
      C->M[4] += B->scal * dist[0] * dist[0] / 2;
      C->M[5] += B->scal * dist[1] * dist[1] / 2;
      C->M[6] += B->scal * dist[2] * dist[2] / 2;
      C->M[7] += B->scal * dist[0] * dist[1] / 2;
      C->M[8] += B->scal * dist[1] * dist[2] / 2;
      C->M[9] += B->scal * dist[2] * dist[0] / 2;
    }                                                           // End loop over all leafs in cell
  }

  void M2M(C_iter C, C_iter P) {                                // Evaluate M2M kernel
    vect dist = P->X - C->X;                                    // Distance between cell centers
    P->M[0] += C->M[0];                                         // Shift multipole expansions
    P->M[1] += C->M[1] +  dist[0]*C->M[0];
    P->M[2] += C->M[2] +  dist[1]*C->M[0];
    P->M[3] += C->M[3] +  dist[2]*C->M[0];
    P->M[4] += C->M[4] +  dist[0]*C->M[1] + dist[0] * dist[0] * C->M[0] / 2;
    P->M[5] += C->M[5] +  dist[1]*C->M[2] + dist[1] * dist[1] * C->M[0] / 2;
    P->M[6] += C->M[6] +  dist[2]*C->M[3] + dist[2] * dist[2] * C->M[0] / 2;
    P->M[7] += C->M[7] + (dist[0]*C->M[2] + dist[1] * C->M[1] + dist[0] * dist[1] * C->M[0]) / 2;
    P->M[8] += C->M[8] + (dist[1]*C->M[3] + dist[2] * C->M[2] + dist[1] * dist[2] * C->M[0]) / 2;
    P->M[9] += C->M[9] + (dist[2]*C->M[1] + dist[0] * C->M[3] + dist[2] * dist[0] * C->M[0]) / 2;
  }

  void M2P(C_iter CI, C_iter CJ) {
    if( CJ->NLEAF >= NCRIT ) {
      for( int i=0; i<CJ->NCHILD; i++ ) {
        C_iter CC = CJ->CHILD[i];
        vect dist = CI->X - CC->X;
        real R = std::sqrt(norm(dist));
        if( CI->R+CC->R > THETA*R ) {
          M2P(CI,CC);
        } else {
          for( B_iter B=CI->LEAF; B!=CI->LEAF+CI->NLEAF; ++B ) {
            dist = B->pos - CC->X;
            R = std::sqrt(norm(dist));
            real R3 = R * R * R;
            real R5 = R3 * R * R;
            B->pot += CC->M[0] / R;
            B->pot += CC->M[1] * (-dist[0] / R3);
            B->pot += CC->M[2] * (-dist[1] / R3);
            B->pot += CC->M[3] * (-dist[2] / R3);
            B->pot += CC->M[4] * (3 * dist[0] * dist[0] / R5 - 1 / R3);
            B->pot += CC->M[5] * (3 * dist[1] * dist[1] / R5 - 1 / R3);
            B->pot += CC->M[6] * (3 * dist[2] * dist[2] / R5 - 1 / R3);
            B->pot += CC->M[7] * (3 * dist[0] * dist[1] / R5);
            B->pot += CC->M[8] * (3 * dist[1] * dist[2] / R5);
            B->pot += CC->M[9] * (3 * dist[2] * dist[0] / R5);
          }
        }
      }
    } else {
      direct(CI->LEAF,CI->LEAF+CI->NLEAF,CJ->LEAF,CJ->LEAF+CJ->NLEAF);
    }
  }


};

#endif
