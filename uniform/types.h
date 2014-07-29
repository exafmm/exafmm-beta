#include <assert.h>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <omp.h>
#include <sys/time.h>

#define for_3d for( int d=0; d<3; d++ )
#define for_4d for( int d=0; d<4; d++ )
#define for_m for( int m=0; m<MTERM; m++ )
#define for_l for( int l=0; l<LTERM; l++ )
#define FMMMAX(a,b) (((a) > (b)) ? (a) : (b))
#define FMMMIN(a,b) (((a) < (b)) ? (a) : (b))

typedef double real;

const int P = 6;
const int MTERM = P*(P+1)*(P+2)/6;
const int LTERM = (P+1)*(P+2)*(P+3)/6;
