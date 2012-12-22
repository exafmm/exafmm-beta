#ifndef _VTGRAPE_H_
#define _VTGRAPE_H_

//#define VG_STRUCT_BOUNDARY 4     // alignment of particle data

#define VG_MAX_FUNCTIONS 100

typedef struct {
  int deviceid;
  int nj;
  int ni;
  int nf;
  void *jvectors;
  int jsize;
  void *ivectors;
  int isize;
  int nati;
  int natj;
  void *matrices;
  int ssize;
  void *scalers;
  void *fvectors;
  int fsize;
  int psize;
  void *pscalers;
  //  double xmax;// only for MD_PERIODIC_FIXED==1
  double volume[3];// only for MD_PERIODIC_FIXED==1
  void (*calculate[VG_MAX_FUNCTIONS])(void *);
  float *r1;
  float *rsqrt;
  //  double *fvec;
  // following 4 variables are for pthread
  double *fthread;
  pthread_t thread;
  int ni_overlap;
  double potc;
  double potv;
  double rcut2;
  double cuton;
  // for GPU overlapping (not thread overlap)
  int gpuoverlapflag;
  int function_index;
  //for MD_QAUNION_ATYPEBIT
  //  float scaleqi_1;
  //  float scaleqj_1;
  //
  int debug_flag;
} VG_UNIT;


#endif
