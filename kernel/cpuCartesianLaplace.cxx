#include "kernel.h"

void Kernel::initialize() {}

void Kernel::P2M() {
  for( B_iter B=CJ->LEAF; B!=CJ->LEAF+CJ->NLEAF; ++B ) {
    vect dist = CJ->X - B->X;
    CJ->M[0] += B->scal;
    CJ->M[1] += B->scal * dist[0];
    CJ->M[2] += B->scal * dist[1];
    CJ->M[3] += B->scal * dist[2];
    CJ->M[4] += B->scal * dist[0] * dist[0] / 2;
    CJ->M[5] += B->scal * dist[1] * dist[1] / 2;
    CJ->M[6] += B->scal * dist[2] * dist[2] / 2;
    CJ->M[7] += B->scal * dist[0] * dist[1];
    CJ->M[8] += B->scal * dist[1] * dist[2];
    CJ->M[9] += B->scal * dist[2] * dist[0];
  }
}

void Kernel::M2M_CPU() {
  vect dist = CI->X - CJ->X;
  CI->M[0] += CJ->M[0];
  CI->M[1] += CJ->M[1] +  dist[0] * CJ->M[0];
  CI->M[2] += CJ->M[2] +  dist[1] * CJ->M[0];
  CI->M[3] += CJ->M[3] +  dist[2] * CJ->M[0];
  CI->M[4] += CJ->M[4] +  dist[0] * CJ->M[1] + dist[0] * dist[0]  * CJ->M[0] / 2;
  CI->M[5] += CJ->M[5] +  dist[1] * CJ->M[2] + dist[1] * dist[1]  * CJ->M[0] / 2;
  CI->M[6] += CJ->M[6] +  dist[2] * CJ->M[3] + dist[2] * dist[2]  * CJ->M[0] / 2;
  CI->M[7] += CJ->M[7] + (dist[0] * CJ->M[2] + dist[1] * CJ->M[1] + dist[0] * dist[1] * CJ->M[0]) / 2;
  CI->M[8] += CJ->M[8] + (dist[1] * CJ->M[3] + dist[2] * CJ->M[2] + dist[1] * dist[2] * CJ->M[0]) / 2;
  CI->M[9] += CJ->M[9] + (dist[2] * CJ->M[1] + dist[0] * CJ->M[3] + dist[2] * dist[0] * CJ->M[0]) / 2;
}

void Kernel::M2L() {
  vect dist = CI->X - CJ->X - Xperiodic;
  real R = std::sqrt(norm(dist));
  real R3 = R * R * R;
  real R5 = R3 * R * R;
  CI->L[0] += CJ->M[0] / R;
  CI->L[0] += CJ->M[1] * (-dist[0] / R3);
  CI->L[0] += CJ->M[2] * (-dist[1] / R3);
  CI->L[0] += CJ->M[3] * (-dist[2] / R3);
  CI->L[0] += CJ->M[4] * (3 * dist[0] * dist[0] / R5 - 1 / R3);
  CI->L[0] += CJ->M[5] * (3 * dist[1] * dist[1] / R5 - 1 / R3);
  CI->L[0] += CJ->M[6] * (3 * dist[2] * dist[2] / R5 - 1 / R3);
  CI->L[0] += CJ->M[7] * (3 * dist[0] * dist[1] / R5);
  CI->L[0] += CJ->M[8] * (3 * dist[1] * dist[2] / R5);
  CI->L[0] += CJ->M[9] * (3 * dist[2] * dist[0] / R5);
  CI->L[1] += CJ->M[0] * (-dist[0] / R3);
  CI->L[1] += CJ->M[1] * (3 * dist[0] * dist[0] / R5 - 1 / R3);
  CI->L[1] += CJ->M[2] * (3 * dist[0] * dist[1] / R5);
  CI->L[1] += CJ->M[3] * (3 * dist[0] * dist[2] / R5);
  CI->L[2] += CJ->M[0] * (-dist[1] / R3);
  CI->L[2] += CJ->M[1] * (3 * dist[1] * dist[0] / R5);
  CI->L[2] += CJ->M[2] * (3 * dist[1] * dist[1] / R5 - 1 / R3);
  CI->L[2] += CJ->M[3] * (3 * dist[1] * dist[2] / R5);
  CI->L[3] += CJ->M[0] * (-dist[2] / R3);
  CI->L[3] += CJ->M[1] * (3 * dist[2] * dist[0] / R5);
  CI->L[3] += CJ->M[2] * (3 * dist[2] * dist[1] / R5);
  CI->L[3] += CJ->M[3] * (3 * dist[2] * dist[2] / R5 - 1 / R3);
  CI->L[4] += CJ->M[0] * (3 * dist[0] * dist[0] / R5 - 1 / R3);
  CI->L[5] += CJ->M[0] * (3 * dist[1] * dist[1] / R5 - 1 / R3);
  CI->L[6] += CJ->M[0] * (3 * dist[2] * dist[2] / R5 - 1 / R3);
  CI->L[7] += CJ->M[0] * (3 * dist[0] * dist[1] / R5);
  CI->L[8] += CJ->M[0] * (3 * dist[1] * dist[2] / R5);
  CI->L[9] += CJ->M[0] * (3 * dist[2] * dist[0] / R5);
}

void Kernel::M2P() {
  for( B_iter B=CI->LEAF; B!=CI->LEAF+CI->NLEAF; ++B ) {
    vect dist = B->X - CJ->X - Xperiodic;
    real R = std::sqrt(norm(dist));
    real R3 = R * R * R;
    real R5 = R3 * R * R;
    B->pot += CJ->M[0] / R;
    B->pot += CJ->M[1] * (-dist[0] / R3);
    B->pot += CJ->M[2] * (-dist[1] / R3);
    B->pot += CJ->M[3] * (-dist[2] / R3);
    B->pot += CJ->M[4] * (3 * dist[0] * dist[0] / R5 - 1 / R3);
    B->pot += CJ->M[5] * (3 * dist[1] * dist[1] / R5 - 1 / R3);
    B->pot += CJ->M[6] * (3 * dist[2] * dist[2] / R5 - 1 / R3);
    B->pot += CJ->M[7] * (3 * dist[0] * dist[1] / R5);
    B->pot += CJ->M[8] * (3 * dist[1] * dist[2] / R5);
    B->pot += CJ->M[9] * (3 * dist[2] * dist[0] / R5);
    B->acc -= dist * CJ->M[0] / R3;
    B->acc[0] += CJ->M[1] * (3 * dist[0] * dist[0] / R5 - 1 / R3);
    B->acc[0] += CJ->M[2] * (3 * dist[0] * dist[1] / R5);
    B->acc[0] += CJ->M[3] * (3 * dist[0] * dist[2] / R5);
    B->acc[1] += CJ->M[1] * (3 * dist[1] * dist[0] / R5);
    B->acc[1] += CJ->M[2] * (3 * dist[1] * dist[1] / R5 - 1 / R3);
    B->acc[1] += CJ->M[3] * (3 * dist[1] * dist[2] / R5);
    B->acc[2] += CJ->M[1] * (3 * dist[2] * dist[0] / R5);
    B->acc[2] += CJ->M[2] * (3 * dist[2] * dist[1] / R5);
    B->acc[2] += CJ->M[3] * (3 * dist[2] * dist[2] / R5 - 1 / R3);
  }
}

void Kernel::P2P() {
  for( B_iter BI=BI0; BI!=BIN; ++BI ) {
    for( B_iter BJ=BJ0; BJ!=BJN; ++BJ ) {
      vect dist = BI->X - BJ->X - Xperiodic;
      real invR = 1 / std::sqrt(norm(dist) + EPS2);
      real invR3 = BJ->scal * invR * invR * invR;
      BI->pot += BJ->scal * invR;
      BI->acc -= dist * invR3;
    }
  }
}

void Kernel::L2L() {
  vect dist = CI->X - CJ->X;
  for( int i=0; i<10; ++i )
    CI->L[i] += CJ->L[i];
  CI->L[0] += CJ->L[1] * dist[0];
  CI->L[0] += CJ->L[2] * dist[1];
  CI->L[0] += CJ->L[3] * dist[2];
  CI->L[0] += CJ->L[4] * dist[0] * dist[0] / 2;
  CI->L[0] += CJ->L[5] * dist[1] * dist[1] / 2;
  CI->L[0] += CJ->L[6] * dist[2] * dist[2] / 2;
  CI->L[0] += CJ->L[7] * dist[0] * dist[1];
  CI->L[0] += CJ->L[8] * dist[1] * dist[2];
  CI->L[0] += CJ->L[9] * dist[2] * dist[0];
  CI->L[1] += CJ->L[4] * dist[0];
  CI->L[1] += CJ->L[7] * dist[1];
  CI->L[1] += CJ->L[9] * dist[2];
  CI->L[2] += CJ->L[7] * dist[0];
  CI->L[2] += CJ->L[5] * dist[1];
  CI->L[2] += CJ->L[8] * dist[2];
  CI->L[3] += CJ->L[9] * dist[0];
  CI->L[3] += CJ->L[8] * dist[1];
  CI->L[3] += CJ->L[6] * dist[2];
}

void Kernel::L2P() {
  for( B_iter B=CI->LEAF; B!=CI->LEAF+CI->NLEAF; ++B ) {
    vect dist = B->X - CI->X;
    B->pot += CI->L[0];
    B->pot += CI->L[1] * dist[0];
    B->pot += CI->L[2] * dist[1];
    B->pot += CI->L[3] * dist[2];
    B->pot += CI->L[4] * dist[0] * dist[0] / 2;
    B->pot += CI->L[5] * dist[1] * dist[1] / 2;
    B->pot += CI->L[6] * dist[2] * dist[2] / 2;
    B->pot += CI->L[7] * dist[0] * dist[1];
    B->pot += CI->L[8] * dist[1] * dist[2];
    B->pot += CI->L[9] * dist[2] * dist[0];
    B->acc[0] += CI->L[1];
    B->acc[0] += CI->L[4] * dist[0];
    B->acc[0] += CI->L[7] * dist[1];
    B->acc[0] += CI->L[9] * dist[2];
    B->acc[1] += CI->L[2];
    B->acc[1] += CI->L[7] * dist[0];
    B->acc[1] += CI->L[5] * dist[1];
    B->acc[1] += CI->L[8] * dist[2];
    B->acc[2] += CI->L[3];
    B->acc[2] += CI->L[9] * dist[0];
    B->acc[2] += CI->L[8] * dist[1];
    B->acc[2] += CI->L[6] * dist[2];
  }
}

void Kernel::finalize() {}
