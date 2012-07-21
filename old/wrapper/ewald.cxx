#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <iostream>
#include <sys/time.h>
#include "mr3.h"

int main(int, char **argv) {
  const int N = 1000;
  const double xmax = 100.0;
  const double ksize = 11.0;
  const double alpha = 0.1;

  int calcmode = atoi(argv[1]);
  int images = atoi(argv[2]);

  // allocate variables
  double *x  = new double [3*N];
  double *xj = new double [3*N];
  double *q  = new double [N];
  double *f  = new double [3*N];
  double *fd = new double [3*N];

  int knum = get_knum(ksize);
  int *ki = new int [3*knum];
  init_kvec(ksize,ki);

  srand48(2);
  // set positions and types
  for( int i=0; i!=N; ++i ) {
    for( int d=0; d!=3; ++d ) {
      x[3*i+d] = drand48() * xmax;
      f[3*i+d] = fd[3*i+d] = 0;
    }
    q[i] = drand48();
  }
  double qsum = 0;
  for( int i=0; i!=N; ++i ) qsum += q[i];
  for( int i=0; i!=N; ++i ) q[i] -= qsum / N;

  // calc with target routine
  switch( calcmode ) {
    case 0: std::cout << "GPU force" << std::endl; break;
    case 1: std::cout << "GPU potential" << std::endl; break;
    case 2: std::cout << "CPU force" << std::endl; break;
    case 3: std::cout << "CPU potential" << std::endl; break;
  }

  // Ewald wave part
  double tpot2, stress[3][3];
  double cell[3][3] = {{xmax,0.0,0.0},{0.0,xmax,0.0},{0.0,0.0,xmax}};
  if( calcmode <= 1 )
    MR3calcewald(ki,((calcmode % 2)==0 ? 1:-1)*knum,x,N,q,alpha,1.0/(M_PI*4.0),cell,f,&tpot2,stress);
  else
    MR3calcewald_host(ki,((calcmode % 2)==0 ? 1:-1)*knum,x,N,q,alpha,1.0/(M_PI*4.0),cell,f,&tpot2,stress);

  // Ewald real part
  if( calcmode <= 1 )
    MR3calccoulomb_ij(N,x,q,f,N,x,q,alpha,6+(calcmode % 2),xmax,1);
  else
    MR3calccoulomb_ij_host(N,x,q,f,N,x,q,alpha,6+(calcmode % 2),xmax,1);

  // Potential energy self term
  if( (calcmode % 2) == 1 ) {
    for( int i=0; i!=N; ++i ) {
      for( int d=0; d!=3; ++d ) {
        f[3*i+d] *= 0.5;
        f[3*i+d] -= q[i] * q[i] * alpha / sqrt(M_PI);
      }
    }
  }

  // Dipole correction
  double fc[3];
  for( int d=0; d!=3; ++d ) fc[d]=0;
  for( int i=0; i!=N; ++i ) {
    for( int d=0; d!=3; ++d ) {
      fc[d] += q[i] * (x[3*i+d] - 0.5 * xmax);
    }
  }
  if( (calcmode % 2) == 1 ) {
    for( int i=0; i!=N; ++i ) {
      for( int d=0; d!=3; ++d ) {
        f[3*i+d] += 2.0 * M_PI / (3.0 * xmax * xmax * xmax)
                  * (fc[0] * fc[0] + fc[1] * fc[1] + fc[2] * fc[2]) / N;
      }
      f[3*i+2] = f[3*i+1] = f[3*i];
    }
  } else {
    for( int i=0; i!=N; ++i ) {
      for( int d=0; d!=3; ++d ) {
        f[3*i+d] -= 4.0 * M_PI * q[i] * fc[d] / (3.0 * xmax * xmax * xmax);
      }
    }
  }

  // Direct with images
  for( int ix=-images; ix<=images; ix++ ) {
    for( int iy=-images; iy<=images; iy++ ) {
      for( int iz=-images; iz<=images; iz++ ) {
        for( int i=0; i!=N; ++i ) {
          xj[3*i+0] = x[3*i+0] + ix * xmax;
          xj[3*i+1] = x[3*i+1] + iy * xmax;
          xj[3*i+2] = x[3*i+2] + iz * xmax;
        }
        if( calcmode <= 1 )
          MR3calccoulomb_ij(N,x,q,fd,N,xj,q,1.0,0+(calcmode % 2),xmax*(2*images+2),0);
        else 
          MR3calccoulomb_ij_host(N,x,q,fd,N,xj,q,1.0,0+(calcmode % 2),xmax*(2*images+2),0);
      }
    }
  }
  for( int i=0; i!=N; ++i ) {
    for( int d=0; d!=3; ++d ) {
      fd[3*i+d] *= q[i];
    }
  }
  if( (calcmode % 2) == 1 ) {
    for( int i=0; i!=N; ++i ) {
      for( int d=0; d!=3; ++d ) {
        fd[3*i+d] *= 0.5;
      }
    }
  }

  // Calculate error
  double diff = 0, norm = 0;
  if( (calcmode % 2) == 1 ) {
    double e = 0, ed = 0;
    for( int i=0; i!=N; ++i ) {
      e += f[3*i];
      ed += fd[3*i];
    }
    diff = (e - ed) * (e - ed);
    norm = ed * ed;
  } else {
    for( int i=0; i!=N; ++i ) {
      for( int d=0; d!=3; ++d ) {
        diff += (f[3*i+d] - fd[3*i+d]) * (f[3*i+d] - fd[3*i+d]);
        norm += fd[3*i+d] * fd[3*i+d];
      }
    }
  }
  std::cout << "Error : " << std::sqrt(diff/norm) << std::endl << std::endl;

  // deallocate variables
  delete[] x;
  delete[] xj;
  delete[] q;
  delete[] f;
  delete[] fd;
  delete[] ki;
  
  return 0;
}
