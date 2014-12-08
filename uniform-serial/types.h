#ifndef types_h
#define types_h

typedef double real;
const int PP = 6;
const int MTERM = PP*(PP+1)*(PP+2)/6;
const int LTERM = (PP+1)*(PP+2)*(PP+3)/6;
#define for_3 for (int d=0; d<3; d++)
#define for_4 for (int d=0; d<4; d++)
#define for_m for (int m=0; m<MTERM; m++)
#define for_l for (int l=0; l<LTERM; l++)
#define MAX(a,b) (((a) > (b)) ? (a) : (b))
#define MIN(a,b) (((a) < (b)) ? (a) : (b))

#endif
