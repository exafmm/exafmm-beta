#include "kernel.h"
#include "laplace.h"

void Kernel::LaplaceInit() {}

#if 1
void Kernel::LaplaceP2M() {
  for( B_iter B=CI->LEAF; B!=CI->LEAF+CI->NLEAF; ++B ) {
    vect dist = CI->X - B->X;
    CI->M[0] += B->SRC[0];
    CI->M[1] += B->SRC[0] * dist[0];
    CI->M[2] += B->SRC[0] * dist[1];
    CI->M[3] += B->SRC[0] * dist[2];
    CI->M[4] += B->SRC[0] * dist[0] * dist[0] / 2;
    CI->M[5] += B->SRC[0] * dist[1] * dist[1] / 2;
    CI->M[6] += B->SRC[0] * dist[2] * dist[2] / 2;
    CI->M[7] += B->SRC[0] * dist[0] * dist[1];
    CI->M[8] += B->SRC[0] * dist[1] * dist[2];
    CI->M[9] += B->SRC[0] * dist[2] * dist[0];
  }
}

void Kernel::LaplaceM2M_CPU() {
  vect dist = CI->X - CJ->X;
  CI->M[0] += CJ->M[0];
  CI->M[1] += CJ->M[1] + dist[0] * CJ->M[0];
  CI->M[2] += CJ->M[2] + dist[1] * CJ->M[0];
  CI->M[3] += CJ->M[3] + dist[2] * CJ->M[0];
  CI->M[4] += CJ->M[4] + dist[0] * CJ->M[1] + dist[0] * dist[0]  * CJ->M[0] / 2;
  CI->M[5] += CJ->M[5] + dist[1] * CJ->M[2] + dist[1] * dist[1]  * CJ->M[0] / 2;
  CI->M[6] += CJ->M[6] + dist[2] * CJ->M[3] + dist[2] * dist[2]  * CJ->M[0] / 2;
  CI->M[7] += CJ->M[7] + dist[0] * CJ->M[2] + dist[1] * CJ->M[1] + dist[0] * dist[1] * CJ->M[0];
  CI->M[8] += CJ->M[8] + dist[1] * CJ->M[3] + dist[2] * CJ->M[2] + dist[1] * dist[2] * CJ->M[0];
  CI->M[9] += CJ->M[9] + dist[2] * CJ->M[1] + dist[0] * CJ->M[3] + dist[2] * dist[0] * CJ->M[0];
}

void Kernel::LaplaceM2L() {
  vect dist = CI->X - CJ->X - Xperiodic;
  real invR = 1 / std::sqrt(norm(dist));
  real invR3 = invR * invR * invR;
  real invR5 = invR3 * invR * invR;
  CI->L[0] += CJ->M[0] * invR;
  CI->L[0] += CJ->M[1] * (-dist[0] * invR3);
  CI->L[0] += CJ->M[2] * (-dist[1] * invR3);
  CI->L[0] += CJ->M[3] * (-dist[2] * invR3);
  CI->L[0] += CJ->M[4] * (3 * dist[0] * dist[0] * invR5 - invR3);
  CI->L[0] += CJ->M[5] * (3 * dist[1] * dist[1] * invR5 - invR3);
  CI->L[0] += CJ->M[6] * (3 * dist[2] * dist[2] * invR5 - invR3);
  CI->L[0] += CJ->M[7] * (3 * dist[0] * dist[1] * invR5);
  CI->L[0] += CJ->M[8] * (3 * dist[1] * dist[2] * invR5);
  CI->L[0] += CJ->M[9] * (3 * dist[2] * dist[0] * invR5);
  CI->L[1] += CJ->M[0] * (-dist[0] * invR3);
  CI->L[1] += CJ->M[1] * (3 * dist[0] * dist[0] * invR5 - invR3);
  CI->L[1] += CJ->M[2] * (3 * dist[0] * dist[1] * invR5);
  CI->L[1] += CJ->M[3] * (3 * dist[0] * dist[2] * invR5);
  CI->L[2] += CJ->M[0] * (-dist[1] * invR3);
  CI->L[2] += CJ->M[1] * (3 * dist[1] * dist[0] * invR5);
  CI->L[2] += CJ->M[2] * (3 * dist[1] * dist[1] * invR5 - invR3);
  CI->L[2] += CJ->M[3] * (3 * dist[1] * dist[2] * invR5);
  CI->L[3] += CJ->M[0] * (-dist[2] * invR3);
  CI->L[3] += CJ->M[1] * (3 * dist[2] * dist[0] * invR5);
  CI->L[3] += CJ->M[2] * (3 * dist[2] * dist[1] * invR5);
  CI->L[3] += CJ->M[3] * (3 * dist[2] * dist[2] * invR5 - invR3);
  CI->L[4] += CJ->M[0] * (3 * dist[0] * dist[0] * invR5 - invR3) / 2;
  CI->L[5] += CJ->M[0] * (3 * dist[1] * dist[1] * invR5 - invR3) / 2;
  CI->L[6] += CJ->M[0] * (3 * dist[2] * dist[2] * invR5 - invR3) / 2;
  CI->L[7] += CJ->M[0] * (3 * dist[0] * dist[1] * invR5);
  CI->L[8] += CJ->M[0] * (3 * dist[1] * dist[2] * invR5);
  CI->L[9] += CJ->M[0] * (3 * dist[2] * dist[0] * invR5);
}

void Kernel::LaplaceM2P() {
  for( B_iter B=CI->LEAF; B!=CI->LEAF+CI->NLEAF; ++B ) {
    vect dist = B->X - CJ->X - Xperiodic;
    real R = std::sqrt(norm(dist));
    real R3 = R * R * R;
    real R5 = R3 * R * R;
    B->TRG[0] += CJ->M[0] / R;
    B->TRG[0] += CJ->M[1] * (-dist[0] / R3);
    B->TRG[0] += CJ->M[2] * (-dist[1] / R3);
    B->TRG[0] += CJ->M[3] * (-dist[2] / R3);
    B->TRG[0] += CJ->M[4] * (3 * dist[0] * dist[0] / R5 - 1 / R3);
    B->TRG[0] += CJ->M[5] * (3 * dist[1] * dist[1] / R5 - 1 / R3);
    B->TRG[0] += CJ->M[6] * (3 * dist[2] * dist[2] / R5 - 1 / R3);
    B->TRG[0] += CJ->M[7] * (3 * dist[0] * dist[1] / R5);
    B->TRG[0] += CJ->M[8] * (3 * dist[1] * dist[2] / R5);
    B->TRG[0] += CJ->M[9] * (3 * dist[2] * dist[0] / R5);
    B->TRG[1] += CJ->M[0] * (-dist[0] / R3);
    B->TRG[1] += CJ->M[1] * (3 * dist[0] * dist[0] / R5 - 1 / R3);
    B->TRG[1] += CJ->M[2] * (3 * dist[0] * dist[1] / R5);
    B->TRG[1] += CJ->M[3] * (3 * dist[0] * dist[2] / R5);
    B->TRG[2] += CJ->M[0] * (-dist[1] / R3);
    B->TRG[2] += CJ->M[1] * (3 * dist[1] * dist[0] / R5);
    B->TRG[2] += CJ->M[2] * (3 * dist[1] * dist[1] / R5 - 1 / R3);
    B->TRG[2] += CJ->M[3] * (3 * dist[1] * dist[2] / R5);
    B->TRG[3] += CJ->M[0] * (-dist[2] / R3);
    B->TRG[3] += CJ->M[1] * (3 * dist[2] * dist[0] / R5);
    B->TRG[3] += CJ->M[2] * (3 * dist[2] * dist[1] / R5);
    B->TRG[3] += CJ->M[3] * (3 * dist[2] * dist[2] / R5 - 1 / R3);
  }
}

void Kernel::LaplaceL2L() {
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

void Kernel::LaplaceL2P() {
  for( B_iter B=CI->LEAF; B!=CI->LEAF+CI->NLEAF; ++B ) {
    vect dist = B->X - CI->X;
    B->TRG[0] += CI->L[0];
    B->TRG[0] += CI->L[1] * dist[0];
    B->TRG[0] += CI->L[2] * dist[1];
    B->TRG[0] += CI->L[3] * dist[2];
    B->TRG[0] += CI->L[4] * dist[0] * dist[0] / 2;
    B->TRG[0] += CI->L[5] * dist[1] * dist[1] / 2;
    B->TRG[0] += CI->L[6] * dist[2] * dist[2] / 2;
    B->TRG[0] += CI->L[7] * dist[0] * dist[1];
    B->TRG[0] += CI->L[8] * dist[1] * dist[2];
    B->TRG[0] += CI->L[9] * dist[2] * dist[0];
    B->TRG[1] += CI->L[1];
    B->TRG[1] += CI->L[4] * dist[0];
    B->TRG[1] += CI->L[7] * dist[1];
    B->TRG[1] += CI->L[9] * dist[2];
    B->TRG[2] += CI->L[2];
    B->TRG[2] += CI->L[7] * dist[0];
    B->TRG[2] += CI->L[5] * dist[1];
    B->TRG[2] += CI->L[8] * dist[2];
    B->TRG[3] += CI->L[3];
    B->TRG[3] += CI->L[9] * dist[0];
    B->TRG[3] += CI->L[8] * dist[1];
    B->TRG[3] += CI->L[6] * dist[2];
  }
}
#else
void Kernel::LaplaceP2M() {
  for( B_iter B=CI->LEAF; B!=CI->LEAF+CI->NLEAF; ++B ) {
    vect dist = CI->X - B->X;
    for( int n=0, nc=0; n<P; ++n ) {
      for( int nx=n; nx>=0; --nx ) {
        for( int nz=0; nz<=n-nx; ++nz, ++nc ) {
          int ny = n - nx - nz;
          CI->M[nc] += B->SRC[0] * pow(dist[0],nx) * pow(dist[1],ny) * pow(dist[2],nz)
            / factorial[nx] / factorial[ny] / factorial[nz];
        }
      }
    }
  }
}

void Kernel::LaplaceM2M_CPU() {
  vect dist = CI->X - CJ->X;
  for( int n=0, nc=0; n<P; ++n ) {
    for( int nx=n; nx>=0; --nx ) {
      for( int nz=0; nz<=n-nx; ++nz, ++nc ) {
        int ny = n - nx - nz;
        real M = 0;
        for( int k=0, kc=0; k<=n; ++k ) {
          for( int kx=k; kx>=0; --kx ) {
            for( int kz=0; kz<=k-kx; ++kz, ++kc ) {
              int ky = k - kx - kz;
              if( kx <= nx && ky <= ny && kz <= nz ) {
                M += CJ->M[kc] * pow(dist[0],nx-kx) * pow(dist[1],ny-ky) * pow(dist[2],nz-kz)
                  / factorial[nx-kx] / factorial[ny-ky] / factorial[nz-kz];
              }
            }
          }
        }
        CI->M[nc] += M;
      }
    }
  }
}

void Kernel::LaplaceM2L() {
  vect dist = CI->X - CJ->X - Xperiodic;
  real invR = 1 / std::sqrt(norm(dist));
  real invR3 = invR * invR * invR;
  real T[P+2][P+2][P+2];
  for( int nx=0; nx<P+2; ++nx ) {
    for( int ny=0; ny<P+2; ++ny ) {
      for( int nz=0; nz<P+2; ++nz ) {
        T[nx][ny][nz] = 0;
      }
    }
  }
  T[2][2][2] = invR;
  T[3][2][2] = -dist[0] * invR3;
  T[2][3][2] = -dist[1] * invR3;
  T[2][2][3] = -dist[2] * invR3;
  for( int n=2; n<P; ++n ) {
    for( int nx=n; nx>=0; --nx ) {
      for( int nz=0; nz<=n-nx; ++nz ) {
        int ny = n - nx - nz;
        int NX = nx+2, NY = ny+2, NZ = nz+2;
        T[NX][NY][NZ] = - ( (2*n-1) *
          (dist[0] * T[NX-1][NY][NZ] + dist[1] * T[NX][NY-1][NZ] + dist[2] * T[NX][NY][NZ-1])
          + (n-1) * (T[NX-2][NY][NZ] + T[NX][NY-2][NZ] + T[NX][NY][NZ-2])) / n * invR * invR;
      }
    }
  }
  for( int n=2; n<P; ++n ) {
    for( int nx=n; nx>=0; --nx ) {
      for( int nz=0; nz<=n-nx; ++nz ) {
        int ny = n - nx - nz;
        int NX = nx+2, NY = ny+2, NZ = nz+2;
        T[NX][NY][NZ] *= factorial[nx] * factorial[ny] * factorial[nz];
      }
    }
  }
  for( int n=0, nc=0; n<P; ++n ) {
    for( int nx=n; nx>=0; --nx ) {
      for( int nz=0; nz<=n-nx; ++nz, ++nc ) {
        int ny = n - nx - nz;
        real L = 0;
        for( int k=0, kc=0; k<P-n; ++k ) {
          for( int kx=k; kx>=0; --kx ) {
            for( int kz=0; kz<=k-kx; ++kz, ++kc ) {
              int ky = k - kx - kz;
              if( nx+kx < P && ny+ky < P && nz+kz < P ) {
                L += CJ->M[kc] * T[nx+kx+2][ny+ky+2][nz+kz+2];
              }
            }
          }
        }
        CI->L[nc] += L;
      }
    }
  }
}

void Kernel::LaplaceM2P() {
  for( B_iter B=CI->LEAF; B!=CI->LEAF+CI->NLEAF; ++B ) {
    vect dist = B->X - CJ->X - Xperiodic;
    real invR = 1 / std::sqrt(norm(dist));
    real invR3 = invR * invR * invR;
    real T[P+3][P+3][P+3];
    for( int nx=0; nx<P+3; ++nx ) {
      for( int ny=0; ny<P+3; ++ny ) {
        for( int nz=0; nz<P+3; ++nz ) {
          T[nx][ny][nz] = 0;
        }
      }
    }
    T[2][2][2] = invR;
    T[3][2][2] = -dist[0] * invR3;
    T[2][3][2] = -dist[1] * invR3;
    T[2][2][3] = -dist[2] * invR3;
    for( int n=2; n<P+1; ++n ) {
      for( int nx=n; nx>=0; --nx ) {
        for( int nz=0; nz<=n-nx; ++nz ) {
          int ny = n - nx - nz;
          int NX = nx+2, NY = ny+2, NZ = nz+2;
          T[NX][NY][NZ] = - ( (2*n-1) *
            (dist[0] * T[NX-1][NY][NZ] + dist[1] * T[NX][NY-1][NZ] + dist[2] * T[NX][NY][NZ-1])
            + (n-1) * (T[NX-2][NY][NZ] + T[NX][NY-2][NZ] + T[NX][NY][NZ-2])) / n * invR * invR;
        }
      }
    }
    for( int n=2; n<P+1; ++n ) {
      for( int nx=n; nx>=0; --nx ) {
        for( int nz=0; nz<=n-nx; ++nz ) {
          int ny = n - nx - nz;
          int NX = nx+2, NY = ny+2, NZ = nz+2;
          T[NX][NY][NZ] *= factorial[nx] * factorial[ny] * factorial[nz];
        }
      }
    }
    vec<4,real> TRG = 0;
    for( int n=0, nc=0; n<P; ++n ) {
      for( int nx=n; nx>=0; --nx ) {
        for( int nz=0; nz<=n-nx; ++nz, ++nc ) {
          int ny = n - nx - nz;
          TRG[0] += CJ->M[nc] * T[nx+2][ny+2][nz+2];
          TRG[1] += CJ->M[nc] * T[nx+3][ny+2][nz+2];
          TRG[2] += CJ->M[nc] * T[nx+2][ny+3][nz+2];
          TRG[3] += CJ->M[nc] * T[nx+2][ny+2][nz+3];
        }
      }
    }
    B->TRG += TRG;
  }
}

void Kernel::LaplaceL2L() {
  vect dist = CI->X - CJ->X;
  for( int n=0, nc=0; n<P; ++n ) {
    for( int nx=n; nx>=0; --nx ) {
      for( int nz=0; nz<=n-nx; ++nz, ++nc ) {
        int ny = n - nx - nz;
        real L = 0;
        for( int k=0, kc=0; k<P; ++k ) {
          for( int kx=k; kx>=0; --kx ) {
            for( int kz=0; kz<=k-kx; ++kz, ++kc ) {
              int ky = k - kx - kz;
              if( nx <= kx && ny <= ky && nz <= kz ) {
                L += CJ->L[kc] * pow(dist[0],kx-nx) * pow(dist[1],ky-ny) * pow(dist[2],kz-nz)
                  / factorial[kx-nx] / factorial[ky-ny] / factorial[kz-nz];
              }
            }
          }
        }
        CI->L[nc] += L;
      }
    }
  }
}

void Kernel::LaplaceL2P() {
  for( B_iter B=CI->LEAF; B!=CI->LEAF+CI->NLEAF; ++B ) {
    vect dist = B->X - CI->X;
    for( int n=0, nc=0; n<P; ++n ) {
      for( int nx=n; nx>=0; --nx ) {
        for( int nz=0; nz<=n-nx; ++nz, ++nc ) {
          int ny = n - nx - nz;
          B->TRG[0] += CI->L[nc] * pow(dist[0],nx) * pow(dist[1],ny) * pow(dist[2],nz)
            / factorial[nx] / factorial[ny] / factorial[nz];
          if( nx > 0 ) B->TRG[1] += CI->L[nc] * pow(dist[0],nx-1) * pow(dist[1],ny) * pow(dist[2],nz)
            / factorial[nx-1] / factorial[ny] / factorial[nz];
          if( ny > 0 ) B->TRG[2] += CI->L[nc] * ny * pow(dist[0],nx) * pow(dist[1],ny-1) * pow(dist[2],nz)
            / factorial[nx] / factorial[ny-1] / factorial[nz];
          if( nz > 0 ) B->TRG[3] += CI->L[nc] * nz * pow(dist[0],nx) * pow(dist[1],ny) * pow(dist[2],nz-1)
            / factorial[nx] / factorial[ny] / factorial[nz-1];
        }
      }
    }
  }
}
#endif

void Kernel::LaplaceFinal() {}
