/*
Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
#ifndef pregpu_h
#define pregpu_h
#include <omp.h>

static size_t  keysDevcSize = 0;                                // Size of offsets for rangeHost
static size_t  rangeDevcSize = 0;                               // Size of offsets for sourceHost
static size_t  sourceDevcSize = 0;                              // Size of sources
static size_t  targetDevcSize = 0;                              // Size of targets
static int     *keysDevc;                                       // Keys on device
static int     *rangeDevc;                                      // Ranges on device
static gpureal *sourceDevc;                                     // Sources on device
static gpureal *targetDevc;                                     // Targets on device
#pragma omp threadprivate(keysDevcSize,rangeDevcSize,sourceDevcSize,targetDevcSize)
#pragma omp threadprivate(keysDevc,rangeDevc,sourceDevc,targetDevc)
__device__ __constant__ gpureal constDevc[1];                   // Constants on device

namespace {                                                     // Limit scope of the following functions to nvcc
__device__ void cart2sph(gpureal& r, gpureal& theta, gpureal& phi,// Get r,theta,phi from x,y,z on GPU
                         gpureal dx, gpureal dy, gpureal dz) {
  r = sqrtf(dx * dx + dy * dy + dz * dz)+EPS;                   // r = sqrt(x^2 + y^2 + z^2) + eps
  theta = acosf(dz / r);                                        // theta = acos(z / r)
  if( fabs(dx) + fabs(dy) < EPS ) {                             // If |x| < eps & |y| < eps
    phi = 0;                                                    //  phi can be anything so we set it to 0
  } else if( fabs(dx) < EPS ) {                                 // If |x| < eps
    phi = dy / fabs(dy) * M_PI * 0.5;                           //  phi = sign(y) * pi / 2
  } else if( dx > 0 ) {                                         // If x > 0
    phi = atanf(dy / dx);                                       //  phi = atan(y / x)
  } else {                                                      // If x < 0
    phi = atanf(dy / dx) + M_PI;                                //  phi = atan(y / x) + pi
  }                                                             // End if for x,y cases
}

__device__ void sph2cart(gpureal r, gpureal theta, gpureal phi, // Spherical to cartesian coordinates on GPU
                         gpureal *spherical, gpureal *cartesian) {
  cartesian[0] = sinf(theta) * cosf(phi) * spherical[0]         // x component (not x itself)
               + cosf(theta) * cosf(phi) / r * spherical[1]
               - sinf(phi) / r / sinf(theta) * spherical[2];
  cartesian[1] = sinf(theta) * sinf(phi) * spherical[0]         // y component (not y itself)
               + cosf(theta) * sinf(phi) / r * spherical[1]
               + cosf(phi) / r / sinf(theta) * spherical[2];
  cartesian[2] = cosf(theta) * spherical[0]                     // z component (not z itself)
               - sinf(theta) / r * spherical[1];
}

__device__ void evalMultipole(gpureal *YnmShrd, gpureal rho,    // Evaluate solid harmonics r^n * Ynm on GPU
                              gpureal alpha, gpureal *factShrd) {
  gpureal x = cosf(alpha);                                      // x = cos(alpha)
  gpureal y = sinf(alpha);                                      // y = sin(alpha)
  gpureal fact = 1;                                             // Initialize 2 * m + 1
  gpureal pn = 1;                                               // Initialize Legendre polynomial Pn
  gpureal rhom = 1;                                             // Initialize rho^m
  for( int m=0; m<P; ++m ){                                     // Loop over m in Ynm
    gpureal p = pn;                                             //  Associate Legendre polynomial Pnm
    int npn = m * m + 2 * m;                                    //  Index of Ynm for m > 0
    int nmn = m * m;                                            //  Index of Ynm for m < 0
    YnmShrd[npn] = rhom * p / factShrd[2*m];                    //  rho^m * Ynm for m > 0
    YnmShrd[nmn] = YnmShrd[npn];                                //  Use conjugate relation for m < 0
    gpureal p1 = p;                                             //  Pnm-1
    p = x * (2 * m + 1) * p;                                    //  Pnm using recurrence relation
    rhom *= -rho;                                               //  rho^m
    gpureal rhon = rhom;                                        //  rho^n
    for( int n=m+1; n<P; ++n ){                                 //  Loop over n in Yn
      int npm = n * n + n + m;                                  //   Index of Ynm for m > 0
      int nmm = n * n + n - m;                                  //   Index of Ynm for m < 0
      YnmShrd[npm] = rhon * p / factShrd[n+m];                  //   rho^n * Ynm
      YnmShrd[nmm] = YnmShrd[npm];                              //   Use conjugate relation for m < 0
      gpureal p2 = p1;                                          //   Pnm-2
      p1 = p;                                                   //   Pnm-1
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);  //   Pnm using recurrence relation
      rhon *= -rho;                                             //   Update rho^n
    }                                                           //  End loop over n in Ynm
    pn = -pn * fact * y;                                        //  Pn
    fact += 2;                                                  //  2 * m + 1
  }                                                             // End loop over m in Ynm
}

__device__ void evalLocal(gpureal *YnmShrd, gpureal rho,        // Evaluate singular harmonics r^(-n-1) * Ynm
                          gpureal alpha, gpureal *factShrd) {
  gpureal x = cosf(alpha);                                      // x = cos(alpha)
  gpureal y = sinf(alpha);                                      // y = sin(alpha)
  gpureal rho_1 = 1 / rho;                                      // 1 / rho
  for( int l=threadIdx.x; l<(2*P+1)*P; l+=THREADS ) {           // Loop over coefficients in Ynm
    gpureal fact = 1;                                           //  Initialize 2 * m + 1
    gpureal pn = 1;                                             //  Initialize Legendre polynomial Pn
    gpureal rhom = rho_1;                                       //  Initialize rho^(-m-1)
    int nn = floor(sqrtf(2*l+0.25)-0.5);                        //  Calculate index n of Ynm
    int mm = 0;                                                 //  Initialize index m of Ynm
    gpureal Ynm;                                                //  Define temporary Ynm
    for( int i=0; i<=nn; ++i ) mm += i;                         //  Offset of m
    mm = l - mm;                                                //  Calculate index m of Ynm
    int n;                                                      //  Define temporary n
    for( int m=0; m<mm; ++m ){                                  //  Loop up to m
      rhom *= rho_1;                                            //   rho^(-m-1)
      pn = -pn * fact * y;                                      //   Pn
      fact += 2;                                                //   2 * m + 1
    }                                                           //  End loop up to m
    int m = mm;                                                 //  Define temporary m
    gpureal p = pn;                                             //  Associated Legendre polynomial Pnm
    if( mm == nn ) Ynm = rhom * p;                              //  Ynm for n == m
    gpureal p1 = p;                                             //  Pnm-1
    p = x * (2 * m + 1) * p;                                    //  Pnm
    rhom *= rho_1;                                              //  rho^(-m-1)
    gpureal rhon = rhom;                                        //  rho^(-n-1)
    for( n=m+1; n<nn; ++n ){                                    //  Loop up to n
      gpureal p2 = p1;                                          //   Pnm-2
      p1 = p;                                                   //   Pnm-1
      p = (x * (2 * n + 1) * p1 - (n + m) * p2) / (n - m + 1);  //   Pnm
      rhon *= rho_1;                                            //   rho^(-n-1)
    }                                                           //  End loop up to n
    if( n <= nn ) Ynm = rhon * p * factShrd[n-m];               //  rho^(-n-1) * Ynm
    YnmShrd[l] = Ynm;                                           //  Put Ynm in shared memory
  }                                                             // End loop over coefficients in Ynm
  __syncthreads();                                              // Syncronize threads
}
}                                                               // End anonymous namespace

#endif
