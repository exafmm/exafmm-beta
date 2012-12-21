#ifndef _MD_H_
#define _MD_H_

#include "vtgrape.h"

//#define CUDA_SDK_2            // define this for CUDA SDK ver. 2

// for AMBER with Narumi's method for HPC Asica 2009 paper.
//            In addition to the AMBER settings,
//            ACCUMULATE_PAIR_FLOAT=2, MD_USE_LARGE_VAL=2
//            __fmul_rn and __fadd_rn are used for LJFB_KERNEL_CORE_SUB1 and
//            LJFB_KERNEL_CORE_SUB12, ADD_PAIR_FLOAT_SIMPLE (two times), and
//            ADD_PAIR_FLOAT.

// for AMBER, undef MD_SORT_ATYPEI, ACCUMULATE_PAIR_FLOAT=12, MD_PERIODIC_FIXED=1
//            MD_LJ_05SIGMA_8SIGMA, VG_MINIMUM_PARTICLE_BLOCK_I=64
//            VG_MINIMUM_ATYPE_BLOCK=(32*32)  

// for Paper for PDCAT08
//            define MD_SORT_ATYPEI
//            VG_MINIMUM_PARTICLE_BLOCK_I=128
//            VG_MINIMUM_ATYPE_BLOCK=256

#define MD_QAUNION_ATYPEBIT     6    // number of bits used for atom type.
                                     // if this is not defined, mixture of atype and
                                     // charge cannot be used.
                                     // MD_USE_QAUNION must be defined.

//#define MR3_MALLOC                   // print malloc and free messages for debugging
#define MD_CELLINDEX                 // support for cell-index
                                     // MD_SORT_ATYPEI must not be defined. 
#define MD_CELLINDEX_MAXCELL_PER_DIM 16
#define MD_CELLINDEX_MAXIBLOCKS      MD_CELLINDEX_MAXCELL_PER_DIM * \
                                     MD_CELLINDEX_MAXCELL_PER_DIM * \
                                     MD_CELLINDEX_MAXCELL_PER_DIM 
                                     // maximum number of cells per dimension
#define MD_USE_LARGE_VAL  20         // use LARGE_VAL as an initial value of summation.
                                     // 0: not use LARGE_VAL, 1: use LARGE_VAL for Nitadori's method
                                     // 2: use LARGE_VAL for my method
                                     // 20: LARGE_VAL for vdw, no LARGE_VAL for coulomb
                                     //     LARGE_VAL is used only for tblno=2 or 3
                                     // For 1 and 2, MD_FVEC_SIZE=3 and 
                                     // defined(ACCUMULATE_PAIR_FLOAT) are required.

//#define USE_EPS_FOR_COULOMB         // this is for gravitational calculation
                                    // which does not check r2=0
#define COULOMB_SHIFT               // support for Coulomb shift function (SHIFT in CHARMM)
#define VDW_SHIFT                   // support for vdw shift function (VSHIFT in CHARMM)
                                    // MD_LJ_05SIGMA_8SIGMA must be 2
//#define MD_SORT_ATYPEI              // sort atom type for i-particles
#define ACCUMULATE_PAIR_FLOAT 2     // define this to use pair float
                                    // 1:pair float full, 2:pair float simple, 3:double
                                    // 12:full for vdw, simple for coulomb
                                    // 10:full for vdw, no pair float for coulomb
                                    // when 3 is used, MD_USE_LARGE_VAL must be 0
                                    // and '-arch sm_13' must be added for compilation flag
#define MD_PERIODIC_FIXED 2         // support periodic boundary condition with fixed point integer
                                    // 1: uses fixed point (faster), 2: use rintf, 0: non-periodic
#define MD_LJ_05SIGMA_8SIGMA 2       // active range of lj interaction is limited
                                      // this must be defined to check with host
                                      // 0:no sigma based cutoff, 1:sigma based cutoff, 2:angstrom based cutoff
                                      // to use 2, MD_SIGMA_EPSION_IN_VGSTRUCT, MD_MATRIX_IN_SCALER and 
                                      // MD_MATRIX_COPY_TO_SHARED must not be defined
//#define MD_LJ_NOATYPE               // define this when no atom type is used
//#define MD_SIGMA_EPSION_IN_VGSTRUCT   // do not use matrix for sigma epsilon parameter
                                      // this can only be used when combination law is 
                                      // satisfied
//#define MD_MATRIX_COPY_TO_SHARED        // gscales and rscales are copyed to shared memory : a little speed-up
                                      // Defining both of MD_SIGMA_EPSION_IN_VGSTRUCT
                                      // and MD_MATRIX_COPY_TO_SHARED is not permitted
#define MD_UNROLL_KERNEL      8       // j loop is unrolled by this number
                                      // if not defined, unrolled up to maximum

#define MD_FVEC_SIZE          3       // size of fi vector (3 is enough for force)
//#define MD_FVEC_SIZE          4       // size of fi vector (3 is enough for force)

#define MD_MAX_I_PER_KERNEL   (1<<20)

#define MD_MAX_J_PER_RUN      (1024*1024) // max number of j particles per run

#define MD_USE_QAUNION                // this must be defined
//#define MD_USE_QATYPE               // this must not be defined
//#define MD_MATRIX_IN_SCALER           // gscales and rscales are in VG_SCALER: this must not be defined

#define MD_R1_EMU_FNAME             "r1.g80emu"
#define MD_RSQRT_EMU_FNAME          "rsqrt.g80emu"

//#define VG_MINIMUM_PARTICLE_BLOCK_I  8  // 64, 128 and 256 work
#define VG_MINIMUM_PARTICLE_BLOCK_I  64  // 64, 128 and 256 work
//#define VG_MINIMUM_PARTICLE_BLOCK_I  128  // 64, 128 and 256 work

//#define VG_JDIV                      1 // number of j-parallelization
#define VG_JDIV                      4 // number of j-parallelization

#define VG_MINIMUM_PARTICLE_BLOCK_J VG_MINIMUM_PARTICLE_BLOCK_I
#define VG_MINIMUM_PARTICLE_BLOCK_F VG_MINIMUM_PARTICLE_BLOCK_I

#define VG_EMU_VECSIZE 4
#define VG_RESTRICTED  __restrict__


//#define VG_MINIMUM_ATYPE_BLOCK       16
//#define VG_MINIMUM_ATYPE_BLOCK      256
//#define VG_MINIMUM_ATYPE_BLOCK       (32*32)
#define VG_MINIMUM_ATYPE_BLOCK       (64*64)

#if 0 // roundup to a little larger value
#define NI_ROUNDUP(n) (((n)+VG_MINIMUM_PARTICLE_BLOCK_I-1)/VG_MINIMUM_PARTICLE_BLOCK_I*VG_MINIMUM_PARTICLE_BLOCK_I+VG_MINIMUM_PARTICLE_BLOCK_I)
#define NJ_ROUNDUP(n) (((n)+VG_MINIMUM_PARTICLE_BLOCK_J-1)/VG_MINIMUM_PARTICLE_BLOCK_J*VG_MINIMUM_PARTICLE_BLOCK_J+VG_MINIMUM_PARTICLE_BLOCK_J)
#define NF_ROUNDUP(n) (((n)+VG_MINIMUM_PARTICLE_BLOCK_F-1)/VG_MINIMUM_PARTICLE_BLOCK_F*VG_MINIMUM_PARTICLE_BLOCK_F+VG_MINIMUM_PARTICLE_BLOCK_F)
#else
#define NI_ROUNDUP(n) (((n)+VG_MINIMUM_PARTICLE_BLOCK_I-1)/VG_MINIMUM_PARTICLE_BLOCK_I*VG_MINIMUM_PARTICLE_BLOCK_I)
#define NJ_ROUNDUP(n) (((n)+VG_MINIMUM_PARTICLE_BLOCK_J-1)/VG_MINIMUM_PARTICLE_BLOCK_J*VG_MINIMUM_PARTICLE_BLOCK_J)
#define NF_ROUNDUP(n) (((n)+VG_MINIMUM_PARTICLE_BLOCK_F-1)/VG_MINIMUM_PARTICLE_BLOCK_F*VG_MINIMUM_PARTICLE_BLOCK_F)
#endif
#define MATRIX_ROUNDUP(n) (((n)+VG_MINIMUM_ATYPE_BLOCK-1)/VG_MINIMUM_ATYPE_BLOCK*VG_MINIMUM_ATYPE_BLOCK)
//#define MATRIX_ROUNDUP(n) (n)


#ifdef MD_USE_QATYPE
#define MD_QATYPE_TO_ATYPE(qatype) ((qatype) & 0xff)
#define MD_QATYPE_TO_Q(qatype)     ((qatype) & ~0xff)
#endif

#ifdef MD_SIGMA_EPSION_IN_VGSTRUCT
#define LJ_BS_SIZE   6
#else
#define LJ_BS_SIZE   4
#endif

#ifndef MASK
#define MASK(n) ((0x1<<(n)) -1)
#endif

#define LARGE_VAL_SHIFT      21    // work of for vdw of scytalone_w24
//#define LARGE_VAL_SHIFT      20
//#define LARGE_VAL_SHIFT      10  // this does not increase the accuracy for Coulomb
#define LOWER_VAL_SHIFT       7  // work ok for scytalone_w24
//#define LOWER_VAL_SHIFT       3
#define LARGE_VAL_ORG          (3<<(LARGE_VAL_SHIFT-1))

// work 
//#define LOWER_VAL_FACTOR_ORG   (1LL<<(23-LARGE_VAL_SHIFT+32-LOWER_VAL_SHIFT))
//#define LOWER_VAL_FACTOR_1_ORG (1.0/((double)LOWER_VAL_FACTOR_ORG))

#define LOWER_VAL_FACTOR_ORG (float)( 1LL << ( 23 - LARGE_VAL_SHIFT + 32 - LOWER_VAL_SHIFT ) )
#define LOWER_VAL_FACTOR_1_ORG ( 1.0f / LOWER_VAL_FACTOR_ORG )

#if MD_USE_LARGE_VAL==1 || MD_USE_LARGE_VAL==2 || MD_USE_LARGE_VAL==20
//
// Note : LARGE_VAL, LOWER_VAL_FACTOR, and LOWER_VAL_FACTOR_1 are
//        modified before including "vtgrape_mixed.c" from "vtgrape.c".
//
#if 0 // work : Narumi's method
#define LARGE_VAL              LARGE_VAL_ORG
#define LOWER_VAL_FACTOR       LOWER_VAL_FACTOR_ORG
#define LOWER_VAL_FACTOR_1     LOWER_VAL_FACTOR_1_ORG
#else // Narumi's method v2
#define LARGE_VAL              ((float)(LARGE_VAL_ORG)*(float)(LOWER_VAL_FACTOR_ORG))
#define LOWER_VAL_FACTOR       1.0f
#define LOWER_VAL_FACTOR_1     1.0f
#endif
#endif

// for QAUNION
#ifdef MD_QAUNION_ATYPEBIT
#define Q_DOUBLE_TO_INT(q,scale)        ((((int)((q)*(scale)*(pow(2,-(MD_QAUNION_ATYPEBIT)))))<<MD_QAUNION_ATYPEBIT) \
                                        & ~(MASK(MD_QAUNION_ATYPEBIT)))
#define Q_INT_TO_FLOAT(q,scale_1)       (((q) & (int)(~(MASK(MD_QAUNION_ATYPEBIT))))*(scale_1))
#define ATYPE_PACK_TO_QATYPE(atype)     ((atype) & MASK(MD_QAUNION_ATYPEBIT))
#define ATYPE_UNPACK_FROM_QATYPE(atype) ((atype) & MASK(MD_QAUNION_ATYPEBIT))
#else
#define ATYPE_UNPACK_FROM_QATYPE(atype) (atype)
#endif


#ifndef TYPEDEF_LONGLONG
#define TYPEDEF_LONGLONG 
# ifdef BIG_ENDIAN
typedef struct { unsigned int h,l; } longlong;
# else
typedef struct { unsigned int l,h; } longlong;
# endif
#endif

#ifndef M3_UNIT
#define M3_UNIT VG_UNIT
#endif

#ifndef M3_CELL
typedef struct{
  long base;
  long size;
} M3_CELL;
#endif

#ifdef MD_USE_QAUNION
typedef union {
  float q;
  int atype;
} VG_QATYPE;
#endif


typedef struct {
  float r[3];
#ifdef MD_USE_QATYPE
  int qatype;
#elif defined(MD_USE_QAUNION)
  VG_QATYPE qatype;
#else
  int atype;
  float q;
#endif
#ifdef MD_SIGMA_EPSION_IN_VGSTRUCT
  float halfsigma;
  float sqrtepsilon;
#endif
} VG_JVEC;

typedef struct {
  float r[3];
#ifdef MD_USE_QATYPE
  int qatype;
#elif defined(MD_USE_QAUNION)
  VG_QATYPE qatype;
#else
  int atype;
  float q;
#endif
#ifdef MD_SIGMA_EPSION_IN_VGSTRUCT
  float halfsigma;
  float sqrtepsilon;
#endif
#if defined(MD_CELLINDEX) && 0
  int ci[4];
#endif
} VG_IVEC;

//#define MD_NACL

typedef struct {
  float gscale;
  float rscale;
#ifdef VDW_SHIFT
  float shift;
#endif
#ifdef MD_NACL
  float pol;
  float sigm;
  float ipotro;
  float pc;
  float pd;
  float zz;
#endif
} VG_MATRIX;

typedef struct {
  float volume[3];
  float alpha;
  float alphafac;
  float alpha3fac;
  float eps2;
#ifdef MD_MATRIX_IN_SCALER
  float gscalesf[VG_MINIMUM_ATYPE_BLOCK];
  float gscalesp[VG_MINIMUM_ATYPE_BLOCK];
  float rscales[VG_MINIMUM_ATYPE_BLOCK];
#endif
#if defined(COULOMB_SHIFT) || 1
  float rcut21;
#endif
#ifdef MD_QAUNION_ATYPEBIT
  float scaleqi_1,scaleqj_1;
#endif
} VG_SCALER;

typedef struct {
  float fi[MD_FVEC_SIZE];
#ifdef ACCUMULATE_PAIR_FLOAT
  float fi2[MD_FVEC_SIZE];
#endif
} VG_FVEC;

typedef struct {
  int base;
  int size;
} VG_CELL;

#define MD_NUM_JCELLS         32

typedef struct {
#ifdef MD_CELLINDEX
  VG_CELL ci[MD_CELLINDEX_MAXIBLOCKS][MD_NUM_JCELLS];
  int   niblock[MD_CELLINDEX_MAXIBLOCKS];
  int   nicell;
#endif
  float potc;
  float potv;
} VG_PSCALER;

typedef union {
  float f;
  int i;
} VG_UNION_FI;


#if MD_LJ_05SIGMA_8SIGMA==2 // angstrom based cutoff
#if 1 // 0.01 - 10A cutoff
#define MD_LJ_R2MIN 0.0001f
#define MD_LJ_R2MAX 100.0f
#elif 0 // 0.3 - 10A cutoff
#define MD_LJ_R2MIN 0.09f
#define MD_LJ_R2MAX 100.0f
#elif 0 // 0.01 - 6A cutoff
#define MD_LJ_R2MIN 0.0001f
#define MD_LJ_R2MAX 36.0f
#elif 0 // 0.01 - 8A cutoff
#define MD_LJ_R2MIN 0.0001f
#define MD_LJ_R2MAX 64.0f
#elif 1 // 0.01 - 12A cutoff
#define MD_LJ_R2MIN 0.0001f
#define MD_LJ_R2MAX 144.0f
#elif 1 // 0.01 - 20A cutoff
#define MD_LJ_R2MIN 0.0001f
#define MD_LJ_R2MAX 400.0f
#elif 1 // 0.01 - 100A cutoff
#define MD_LJ_R2MIN 0.0001f
#define MD_LJ_R2MAX 10000.0f
#endif
#else // sigma based cutoff
#if 1 // default of AMBER, MDM_TABLER_VER=00,01
#define MD_LJ_R2MIN (1.0f/(64.0f*64.0f))
#define MD_LJ_R2MAX (1024.0f*1024.0f)
#endif
#if 0
#define MD_LJ_R2MIN (1.0f/(64.0f*64.0f))
#define MD_LJ_R2MAX (16.0f)
#endif
#if 0 // default of MDGRAPE,-2,-3
#define MD_LJ_R2MIN 0.25f
#define MD_LJ_R2MAX 64.0f
#endif
#endif


#define MD_REAL_R2MIN 0.0001f
#define MD_REAL_R2MAX 10.0f


#ifdef MR3_MALLOC
#define MR3_malloc_pointer(x,y) MR3_my_malloc2(x,y)
#define MR3_free_pointer(x,y)   MR3_my_free2((void **)(&(x)),y)
#else
#define MR3_malloc_pointer(x,y) malloc(x)
#define MR3_free_pointer(x,y)   free(x)
#endif

//#define VG_GCC                    // for GCC hand vectorization

#if defined(VG_GCC) && __GNUC__ == 4
typedef float v4sf __attribute__ ((vector_size(16)));
typedef int   v4si __attribute__ ((vector_size(16)));
#else
typedef float* v4sf;
#endif


typedef struct {
  int ni;
  double *xi;
  double *qi;
  int *atypei;
  double *force;
  int nj;
  double *xj;
  double *qj; 
  int *atypej;
  int nat;
  double *gscalesf;
  double *gscalesp;
  double *rscales;
  double rscale;
  int tblno;
  double xmax;
  double *potc;
  double *potv;
  int periodicflag;
  int potflag;
  int changeflag;
  int *numex;
  int *natex;
} thread_arg_coulombvdwijexlist;


#define grav \
        rr2*=alpha2;\
	ftmp=jvec[jj].q*powf(rr2,-1.5f);\
	for(k=0;k<3;k++) fvec[i].fi[k]+=ftmp*dr[k];
#define gravpot \
        rr2*=alpha2;\
        ftmp=jvec[jj].q/sqrtf(rr2);\
        fvec[i].fi[0]+=ftmp;
#define lj \
        at=ivec[i].atype*natj+jvec[jj].atype;\
        rr2*=matrix[at].rscale;\
        if(rr2>=MD_LJ_R2MIN && rr2<MD_LJ_R2MAX){\
          dn2=1.0f/rr2;\
          dn6=dn2*dn2*dn2;\
          ftmp=matrix[at].gscale*dn6*dn2*(2.0f*dn6-1.0f);\
          for(k=0;k<3;k++) fvec[i].fi[k]+=ftmp*dr[k];\
	}
#define ljpot \
        at=ivec[i].atype*natj+jvec[jj].atype;\
        rr2*=matrix[at].rscale;\
        if(rr2>=MD_LJ_R2MIN && rr2<MD_LJ_R2MAX){\
          dn2=1.0f/rr2;\
          dn6=dn2*dn2*dn2;\
          ftmp=matrix[at].gscale*dn6*(dn6-1.0f);\
          fvec[i].fi[0]+=ftmp;\
	}
#define real \
        rr2*=alpha2;\
        sqdn=sqrtf(rr2);\
        ftmp=jvec[jj].q*(api2*expf(-rr2)+erfcf(sqdn)/sqdn)/rr2; \
	for(k=0;k<3;k++) fvec[i].fi[k]+=ftmp*dr[k];
#define realpot \
        rr2*=alpha2;\
        sqdn=sqrtf(rr2);\
        ftmp=jvec[jj].q*erfcf(sqdn)/sqdn;\
        fvec[i].fi[0]+=ftmp;
#define wave \
        rr2*=alpha2;
#define wavepot \
        rr2*=alpha2;
#define dummy \
        rr2*=alpha2;\
        if(rr2>=MD_REAL_R2MIN && rr2<MD_REAL_R2MAX){\
          sqdn=sqrtf(rr2);\
          ftmp=jvec[jj].q*erfcf(sqdn)/sqdn;\
          if(i==3 && rr2<10.0f){\
            printf("vg i=%d j=%d r2=%f ivec=%f %f %f jvec=%f %f %f\n",i,j,rr2,ivec[i].r[0],ivec[i].r[1],ivec[i].r[2],jvec[j].r[0],jvec[j].r[1],jvec[j].r[2]); \
	    printf("   vol=%f %f %f dr=%f %f %f ftmp=%e\n",vol[0],vol[1],vol[2],dr[0],dr[1],dr[2],ftmp); \
	    printf("   factor=%e ftmp*factor=%e\n",ivec[i].q*sqrt(alpha2),ivec[i].q*sqrt(alpha2)*ftmp); \
	    printf("   fvec=%e %e %e\n",fvec[i].fi[0],fvec[i].fi[1],fvec[i].fi[2]);\
          }\
          fvec[i].fi[0]+=ftmp;\
	  if(i==3 && rr2<10.0f) printf("   fvec=%e %e %e\n",fvec[i].fi[0],fvec[i].fi[1],fvec[i].fi[2]); \
	}

//  for(i=0;i<VG_MINIMUM_PARTICLE_BLOCK_I;i++){	
//    for(j=0;j<VG_MINIMUM_PARTICLE_BLOCK_J;j++){	

#endif
