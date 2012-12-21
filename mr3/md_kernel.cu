/*
  Following part is crucial to FMA operations which are not 
  IEEE754 compliant.

  #define LJFB_KERNEL_CORE_SUB1 
      rij=dx*dx+dy*dy+dz*dz; -- original 
      rij=__fadd_rn(__fadd_rn(__fmul_rn(dx,dx),__fmul_rn(dy,dy)),__fmul_rn(dz,dz));
                             -- modified to prevent fma
  #define LJFB_KERNEL_CORE_SUB12
      similar modification with LJFB_KERNEL_CORE_SUB1 

  #define ADD_PAIR_FLOAT(xv,a,axch,axcl,tmp2v,rinvv,mrinvv,rinv2v) 
      tmp2v=xv*a;            -- original
      tmp2v=__fmul_rn(xv,a); -- modified
  #define ADD_PAIR_FLOAT_SIMPLE
      similar modification with ADD_PAIR_FLOAT

  This seems not effective
  #define LJFB_KERNEL_CORE_SUB3 
      tmp=d_matrix[itype].gscale*dn6*cutoff*(2.0f*dn6-1.0f); 
                             -- original
      tmp=d_matrix[itype].gscale*dn6*cutoff*(__fmul_rn(2.0f,dn6)-1.0f);
                             -- modified
 */

#include <stdio.h>

#include "md.h"

/* 
   SHIFT_FUNC   1: vdw force shift   - VSHIFT
                2: vdw force switch  - VSWITCH
              101: vdw energy shift  - VSHIFT
              102: vdw energy switch - VSWITCH
*/
template <int SHIFT_FUNC>
__inline__ __device__ float vdw_shift(float rr, float r2max, float rcut21, float r2rscale_1, float shift, float gscale)
{
  float tmp,dn6,rij;
  dn6=r2rscale_1*r2rscale_1*r2rscale_1;
  switch(SHIFT_FUNC){
  case 0:   // normal Lennard Jones force
    tmp=VDW_FORCE(r2rscale_1,dn6,1.0f,2.0f,__fmul_rn,__fadd_rn);
    break;
  case 100: // normal Lennard Jones energy
    tmp=VDW_ENERGY(dn6,1.0f,__fmul_rn,__fadd_rn);
    break;
  case 1:   // VSHIFT force
    tmp=VDW_FORCE_SHIFT(r2rscale_1,dn6,rr,rcut21,shift,1.0f,2.0f,3.0f,__fmul_rn,__fadd_rn,rij);
    break;
  case 101: // VSHIFT energy
    tmp=VDW_ENERGY_SHIFT(r2rscale_1,dn6,rr,rcut21,shift,1.0f,2.0f,3.0f,__fmul_rn,__fadd_rn,rij);
    break;
  case 2:   // VSWITCH force
    if(rr<r2max){
      if(rr>=shift){
	tmp = VDW_FORCE_SWITCH(r2rscale_1,dn6,rr,rcut21,shift,1.0f,2.0f,3.0f,__fmul_rn,__fadd_rn,tmp,rij);
      }
      else{
	tmp = VDW_FORCE(r2rscale_1,dn6,1.0f,2.0f,__fmul_rn,__fadd_rn);
      }
    }
    break;
  case 102: // VSWITCH energy
    if(rr<r2max){
      if(rr>=shift){
	tmp = VDW_ENERGY_SWITCH(r2rscale_1,dn6,rr,shift,1.0f,2.0f,3.0f,__fmul_rn,__fadd_rn,rij);
      }
      else{
	tmp = VDW_ENERGY(dn6,1.0f,__fmul_rn,__fadd_rn);
      }
    }
    break;
  }
  
  if(rr<MD_LJ_R2MIN)  tmp=0.0f;
  if(rr>=r2max)       tmp=0.0f;
  tmp*=gscale;
  
  return tmp;
}


#define CHECK_BANK_CONFLICTS 0
#if CHECK_BANK_CONFLICTS
#define AS(i, j) CUT_BANK_CHECKER(((float*)&As[0][0]), (BLOCK_SIZE * i + j))
#define BS(i, j) CUT_BANK_CHECKER(((float*)&Bs[0][0]), (BLOCK_SIZE * i + j))
#define DS(i, j) CUT_BANK_CHECKER(((float*)&Ds[0][0]), (BLOCK_SIZE * i + j))
#define ES(i, j) CUT_BANK_CHECKER(((float*)&Es[0][0]), (BLOCK_SIZE * i + j))
#define FS(i, j) CUT_BANK_CHECKER(((float*)&Fs[0][0]), (BLOCK_SIZE * i + j))
#define GS(i, j) CUT_BANK_CHECKER(((float*)&Gs[0][0]), (BLOCK_SIZE * i + j))
#define QIS(i, j) CUT_BANK_CHECKER(((float*)&Qis[0][0]), (BLOCK_SIZE * i + j))
#define QJS(i, j) CUT_BANK_CHECKER(((float*)&Qjs[0][0]), (BLOCK_SIZE * i + j))
#else
#define AS(i, j) As[i][j]
#define BS(i, j) Bs[i][j]
#define DS(i, j) Ds[i][j]
#define ES(i, j) Es[i][j]
#define FS(i, j) Fs[i][j]
#define GS(i, j) Gs[i][j]
#define QIS(i, j) Qis[i][j]
#define QJS(i, j) Qjs[i][j]
#endif


#ifdef MD_USE_CONSTANT
__device__ __constant__ VG_MATRIX d_matrix[VG_MINIMUM_ATYPE_BLOCK];
__device__ __constant__ VG_SCALER d_scalers[1];
#endif

#if 1

#define LJ_KERNEL_CORE \
      dn2 = (dx * dx + dy * dy + dz * dz) * d_matrix[itype].rscale;\
      if(j*VG_MINIMUM_PARTICLE_BLOCK_I+jj<nj && dn2 >= MD_LJ_R2MIN && dn2<MD_LJ_R2MAX && dn2!=0.0f){\
	float r2inv;\
        r2inv=1.0f/dn2;\
	dn6=r2inv*r2inv*r2inv;\
	tmp = d_matrix[itype].gscale * dn6 * r2inv * (2.e0 * dn6 - 1.e0);\
	Csub[0] += tmp * dx;\
	Csub[1] += tmp * dy;\
	Csub[2] += tmp * dz;\
      }
#define LJPOT_KERNEL_CORE \
      dn2 = (dx * dx + dy * dy + dz * dz) * d_matrix[itype].rscale;\
      if(j*VG_MINIMUM_PARTICLE_BLOCK_I+jj<nj && dn2 >= MD_LJ_R2MIN && dn2<MD_LJ_R2MAX && dn2!=0.0f){\
	float r2inv;\
        r2inv=1.0f/dn2;\
	dn6=r2inv*r2inv*r2inv;\
	tmp = d_matrix[itype].gscale * dn6 * (dn6 - 1.e0);\
	Csub[0] += tmp;\
      }
#if 1 // normal grav 
#define GRAV_KERNEL_CORE \
      dn2 = (dx * dx + dy * dy + dz * dz + eps2) * alpha2;\
      if(j*VG_MINIMUM_PARTICLE_BLOCK_I+jj<nj && dn2!=0.0f){\
        float rinv;\
        rinv=rsqrtf(dn2);\
	tmp = QJS(0,jj)*rinv*rinv*rinv;\
	Csub[0] += tmp * dx;\
	Csub[1] += tmp * dy;\
	Csub[2] += tmp * dz;\
      }
#else // only work for gravitational simulation
#define GRAV_KERNEL_CORE \
      dn2 = (dx * dx + dy * dy + dz * dz) + eps2;\
      float rinv;\
      rinv=rsqrtf(dn2);\
      tmp = QJS(0,jj)*rinv*rinv*rinv;\
      Csub[0] += tmp * dx;\
      Csub[1] += tmp * dy;\
      Csub[2] += tmp * dz;
#endif
#define GRAVPOT_KERNEL_CORE \
      dn2 = (dx * dx + dy * dy + dz * dz) * alpha2;\
      if(j*VG_MINIMUM_PARTICLE_BLOCK_I+jj<nj && dn2!=0.0f){\
        float rinv;\
        rinv=rsqrtf(dn2);\
	tmp = QJS(0,jj)*rinv;\
	Csub[0] += tmp;\
      }
#if MD_LJ_05SIGMA_8SIGMA==2
#define REAL_KERNEL_CORE_ISR2ZERO \
      if(r2<MD_LJ_R2MIN)  dn2=0.0f;\
      if(r2>=MD_LJ_R2MAX) dn2=0.0f;
#define REALFB_KERNEL_CORE_ISR2ZERO \
      if(rij<MD_LJ_R2MIN)  dn6=0.0f;\
      if(rij>=r2max)       dn6=0.0f;
#else
#define REAL_KERNEL_CORE_ISR2ZERO \
      if(r2==0.0f) dn2=0.0f;
#define REALFB_KERNEL_CORE_ISR2ZERO \
      if(rij==0.0f) dn6=0.0f;
#endif
#define REAL_KERNEL_CORE \
      dn2 = dx * dx + dy * dy + dz * dz;\
      if(j*VG_MINIMUM_PARTICLE_BLOCK_I+jj<nj && dn2!=0.0f){\
        float sqdn,r2=dn2;\
        dn2*=alpha2;\
        sqdn=sqrtf(dn2);\
	tmp = QJS(0,jj)*(api2*expf(-dn2)+erfcf(sqdn)/sqdn);\
        dn2=1.0f/dn2;\
        REAL_KERNEL_CORE_ISR2ZERO;\
        tmp*=dn2;\
	Csub[0] += tmp * dx;\
	Csub[1] += tmp * dy;\
	Csub[2] += tmp * dz;\
      }
#define REALPOT_KERNEL_CORE \
      dn2 = dx * dx + dy * dy + dz * dz;\
      if(j*VG_MINIMUM_PARTICLE_BLOCK_I+jj<nj){\
        float sqdn,r2=dn2;\
        dn2*=alpha2;\
        sqdn=sqrtf(dn2);\
	tmp = QJS(0,jj)*erfcf(sqdn);\
        dn2=1.0f/sqdn;\
        REAL_KERNEL_CORE_ISR2ZERO;\
	tmp*=dn2;\
	Csub[0] += tmp;\
      }
#define WAVE_KERNEL_CORE \
      dn2 = 0.0f;
#define WAVEPOT_KERNEL_CORE \
      dn2 = 0.0f;
#define NACL_KERNEL_CORE \
      dn2 = (dx * dx + dy * dy + dz * dz);\
      if(j*VG_MINIMUM_PARTICLE_BLOCK_I+jj<nj && dn2!=0.0f){\
	float r,inr,inr2,inr4,inr8,d3,pb=0.338e-19/(14.39*1.60219e-19),dphir; \
        r=sqrtf(dn2);\
	inr=1.0f/r;\
	inr2=inr*inr;\
	inr4=inr2*inr2;\
	inr8=inr4*inr4;\
	d3=pb*d_matrix[itype].pol*expf((d_matrix[itype].sigm-r)*d_matrix[itype].ipotro);\
	dphir=(d3*d_matrix[itype].ipotro*inr\
			 - 6.0f*d_matrix[itype].pc*inr8\
			 - 8.0f*d_matrix[itype].pd*inr8*inr2\
			 + inr*inr*inr*d_matrix[itype].zz);\
	Csub[0] += dphir * dx;\
	Csub[1] += dphir * dy;\
	Csub[2] += dphir * dz;\
      }
#define R1_KERNEL_CORE \
      dn2 = dx;\
      if(j*VG_MINIMUM_PARTICLE_BLOCK_I+jj<nj && dn2!=0.0f){\
	float r2inv;\
        r2inv=QJS(0,jj)/dn2;\
	Csub[0] += r2inv;\
      }
#define RSQRT_KERNEL_CORE \
      dn2 = dx;\
      if(j*VG_MINIMUM_PARTICLE_BLOCK_I+jj<nj && dn2!=0.0f){\
	float r2inv;\
        r2inv=QJS(0,jj)*rsqrtf(dn2);\
	Csub[0] += r2inv;\
      }
#define MUL_KERNEL_CORE \
      dn2 = dx;\
      if(j*VG_MINIMUM_PARTICLE_BLOCK_I+jj<nj && dn2!=0.0f){\
	float r2inv;\
        r2inv=QJS(0,jj)+dn2*dn2;\
	Csub[0] += r2inv;\
      }

#if MD_PERIODIC_FIXED==1
#if 0 // use al3
#define DISTANCE_PERIODIC \
      dx = xi[0] - __float_as_int(BS(0,jj4  ));\
      dy = xi[1] - __float_as_int(BS(0,jj4+1));\
      dz = xi[2] - __float_as_int(BS(0,jj4+2));\
      dx *= al3;\
      dy *= al3;\
      dz *= al3;
#define SCALE_R2 ;
#else // else of use al3
#define DISTANCE_PERIODIC \
      dx = xi[0] - __float_as_int(BS(0,jj4  ));\
      dy = xi[1] - __float_as_int(BS(0,jj4+1));\
      dz = xi[2] - __float_as_int(BS(0,jj4+2));\
      dx *= al2[0];\
      dy *= al2[1];\
      dz *= al2[2];
#define SCALE_R2 ;
#endif // end of use al3
#elif MD_PERIODIC_FIXED==2
#define DISTANCE_PERIODIC \
      dx = xi[0] - BS(0,jj4  );\
      dy = xi[1] - BS(0,jj4+1);\
      dz = xi[2] - BS(0,jj4+2);\
      dx = dx - rintf(dx * al2[0]) * size[0];\
      dy = dy - rintf(dy * al2[1]) * size[1];\
      dz = dz - rintf(dz * al2[2]) * size[2];
#define SCALE_R2 ;
#elif MD_PERIODIC_FIXED==0
#define DISTANCE_PERIODIC \
      dx = xi[0] - BS(0,jj4  );\
      dy = xi[1] - BS(0,jj4+1);\
      dz = xi[2] - BS(0,jj4+2);
#define SCALE_R2 ;
#endif

#define ADD_PAIR_FLOAT_SIMPLE(xv,a,axch,axcl,tmp2v,rinvv,mrinvv,rinv2v) \
      /*tmp2v=xv*a;*/\
      tmp2v=__fmul_rn(xv,a);\
      mrinvv=axch;\
      axch=axch+tmp2v;\
      rinvv=axch-mrinvv;\
      rinvv=tmp2v-rinvv;\
      axcl=axcl+rinvv;
#define ADD_PAIR_FLOAT_NARUMI2(xv,a,axch,axcl,tmp2v,rinvv,mrinvv,rinv2v) \
      /*tmp2v=xv*a;*/\
      tmp2v=__fmul_rn(xv,a);\
      mrinvv=axch;\
      axch=axch+tmp2v;\
      rinvv=axch-mrinvv;\
      rinvv=tmp2v-rinvv;\
      axcl=__int_as_float(__float_as_int(axcl)+__float2int_rn(rinvv*((float)LOWER_VAL_FACTOR)));
#define ADD_PAIR_FLOAT(xv,a,axch,axcl,tmp2v,rinvv,mrinvv,rinv2v) \
      /* tmp2v, rinvv, mrinvv, rinv2v are temporary variables*/\
      /*tmp2v=xv*a;*/\
      tmp2v=__fmul_rn(xv,a);\
      rinvv=axch+tmp2v;\
      mrinvv=rinvv-axch;\
      rinv2v=rinvv-mrinvv;\
      rinv2v=axch-rinv2v;\
      mrinvv=tmp2v-mrinvv;\
      rinv2v=mrinvv+rinv2v;\
      rinv2v=rinv2v+axcl;\
      axch=rinvv+rinv2v;\
      axcl=axch-rinvv;\
      axcl=rinv2v-axcl;
#define ADD_PAIR_FLOAT_DOUBLE(xv,a,axch,tmp2v,rinvv,mrinvv,rinv2v) \
      /*tmp2v=xv*a;*/\
      tmp2v=__fmul_rn(xv,a);\
      axch+=tmp2v;\
      /*axch+=1.0;*/\
      /*axch=0.0*/;


#ifdef ACCUMULATE_PAIR_FLOAT
#if ACCUMULATE_PAIR_FLOAT==1 
#define CLFB_KERNEL_CORE_ACCUM \
      ADD_PAIR_FLOAT(dx,tmp,Csub[0],Csub2[0],dx,tmp2,tmp3,tmp4);\
      ADD_PAIR_FLOAT(dy,tmp,Csub[1],Csub2[1],dy,tmp2,tmp3,tmp4);\
      ADD_PAIR_FLOAT(dz,tmp,Csub[2],Csub2[2],dz,tmp2,tmp3,tmp4);
#define CLFB_KERNEL_CORE_ACCUMPOT \
      ADD_PAIR_FLOAT(1.0f,tmp,Csub[0],Csub2[0],dx,tmp2,tmp3,tmp4);
#elif ACCUMULATE_PAIR_FLOAT==2 || ACCUMULATE_PAIR_FLOAT==12
#if MD_USE_LARGE_VAL==2
#define CLFB_KERNEL_CORE_ACCUM \
      ADD_PAIR_FLOAT_NARUMI2(dx,tmp,Csub[0],Csub2[0],dx,tmp2,tmp3,tmp4);\
      ADD_PAIR_FLOAT_NARUMI2(dy,tmp,Csub[1],Csub2[1],dy,tmp2,tmp3,tmp4);\
      ADD_PAIR_FLOAT_NARUMI2(dz,tmp,Csub[2],Csub2[2],dz,tmp2,tmp3,tmp4);
#define CLFB_KERNEL_CORE_ACCUMPOT \
      ADD_PAIR_FLOAT_NARUMI2(1.0f,tmp,Csub[0],Csub2[0],dx,tmp2,tmp3,tmp4);
#else
#define CLFB_KERNEL_CORE_ACCUM \
      ADD_PAIR_FLOAT_SIMPLE(dx,tmp,Csub[0],Csub2[0],dx,tmp2,tmp3,tmp4);\
      ADD_PAIR_FLOAT_SIMPLE(dy,tmp,Csub[1],Csub2[1],dy,tmp2,tmp3,tmp4);\
      ADD_PAIR_FLOAT_SIMPLE(dz,tmp,Csub[2],Csub2[2],dz,tmp2,tmp3,tmp4);
#define CLFB_KERNEL_CORE_ACCUMPOT \
      ADD_PAIR_FLOAT_SIMPLE(1.0f,tmp,Csub[0],Csub2[0],dx,tmp2,tmp3,tmp4);
#endif
#elif ACCUMULATE_PAIR_FLOAT==10
#define CLFB_KERNEL_CORE_ACCUM \
      Csub[0]+=tmp*dx;\
      Csub[1]+=tmp*dy;\
      Csub[2]+=tmp*dz;
#define CLFB_KERNEL_CORE_ACCUMPOT \
      Csub[0]+=tmp;
#elif ACCUMULATE_PAIR_FLOAT==3
#define CLFB_KERNEL_CORE_ACCUM \
      ADD_PAIR_FLOAT_DOUBLE(dx,tmp,Csubd[0],dx,tmp2,tmp3,tmp4);\
      ADD_PAIR_FLOAT_DOUBLE(dy,tmp,Csubd[1],dy,tmp2,tmp3,tmp4);\
      ADD_PAIR_FLOAT_DOUBLE(dz,tmp,Csubd[2],dz,tmp2,tmp3,tmp4);
#define CLFB_KERNEL_CORE_ACCUMPOT \
      ADD_PAIR_FLOAT_DOUBLE(1.0f,tmp,Csubd[0],dx,tmp2,tmp3,tmp4);
#endif
#else
#define CLFB_KERNEL_CORE_ACCUM \
      Csub[0]+=tmp*dx;\
      Csub[1]+=tmp*dy;\
      Csub[2]+=tmp*dz;
#define CLFB_KERNEL_CORE_ACCUMPOT \
      Csub[0]+=tmp;
#endif

// Coulomb kernel
#if MD_LJ_05SIGMA_8SIGMA==2
#define GRAVFB_KERNEL_CORE_EPS 
#define GRAVFB_KERNEL_CORE_ISR2ZERO \
      if(rij<MD_LJ_R2MIN)  cutoff=0.0f;\
      if(rij>=r2max)       cutoff=0.0f;
#else // else of MD_LJ_05SIGMA_8SIGMA==2
#ifdef USE_EPS_FOR_COULOMB
#define GRAVFB_KERNEL_CORE_EPS rij+=eps
#define GRAVFB_KERNEL_CORE_ISR2ZERO 
#else
#define GRAVFB_KERNEL_CORE_EPS 
#define GRAVFB_KERNEL_CORE_ISR2ZERO \
      if(rij==0.0f) cutoff=0.0f;\
      /*if(rij>=400.0f) cutoff=0.0f;*/
#endif
#endif // end of MD_LJ_05SIGMA_8SIGMA==2

#ifdef MD_QAUNION_ATYPEBIT
#define GRAVFB_KERNEL_CORE_SHIFT_QJ(x) (Q_INT_TO_FLOAT(__float_as_int(x),d_scalers->scaleqj_1))
#else
#define GRAVFB_KERNEL_CORE_SHIFT_QJ(x) (x)
#endif

#ifdef COULOMB_SHIFT
#define GRAVFB_KERNEL_CORE_SHIFT \
      tmp=GRAVFB_KERNEL_CORE_SHIFT_QJ(BS(0,jj4+3))*cutoff;	\
      rij=1.0f-__fmul_rn(rij,rcut21);\
      /*tmp=__fmul_rn(tmp*rij*rij*cutoff,cutoff)+__fmul_rn(4.0f*rcut21*tmp,rij);*/ \
      tmp*=rij*(__fmul_rn(rij*cutoff,cutoff)+__fmul_rn(4.0f,rcut21));
#define GRAVPOTFB_KERNEL_CORE_SHIFT \
      tmp=GRAVFB_KERNEL_CORE_SHIFT_QJ(BS(0,jj4+3))*cutoff;	\
      rij=1.0f-__fmul_rn(rij,rcut21);\
      tmp*=rij*rij;
#else
#define GRAVFB_KERNEL_CORE_SHIFT \
      tmp=GRAVFB_KERNEL_CORE_SHIFT_QJ(BS(0,jj4+3))*cutoff*cutoff*cutoff;
#define GRAVPOTFB_KERNEL_CORE_SHIFT \
      tmp=GRAVFB_KERNEL_CORE_SHIFT_QJ(BS(0,jj4+3))*cutoff;
#endif

#define GRAVFB_KERNEL_CORE \
      DISTANCE_PERIODIC;\
      rij=dx*dx+dy*dy+dz*dz;\
      GRAVFB_KERNEL_CORE_EPS;\
      SCALE_R2;\
      cutoff=rsqrtf(rij);\
      GRAVFB_KERNEL_CORE_ISR2ZERO;\
      GRAVFB_KERNEL_CORE_SHIFT;\
      jj4+=LJ_BS_SIZE;\
      CLFB_KERNEL_CORE_ACCUM;

#define GRAVPOTFB_KERNEL_CORE \
      DISTANCE_PERIODIC;\
      rij=dx*dx+dy*dy+dz*dz;\
      GRAVFB_KERNEL_CORE_EPS;\
      SCALE_R2;\
      cutoff=rsqrtf(rij);\
      GRAVFB_KERNEL_CORE_ISR2ZERO;\
      GRAVPOTFB_KERNEL_CORE_SHIFT;\
      jj4+=LJ_BS_SIZE;\
      CLFB_KERNEL_CORE_ACCUMPOT;

#define REALFB_KERNEL_CORE \
      DISTANCE_PERIODIC;\
      rij=dx*dx+dy*dy+dz*dz;\
      GRAVFB_KERNEL_CORE_EPS;\
      SCALE_R2;\
      /*cutoff=rsqrtf(rij);*/\
      dn6=rij*alpha2;\
      cutoff=sqrtf(dn6);\
      tmp=BS(0,jj4+3)*(api2*expf(-dn6)+erfcf(cutoff)/cutoff);\
      dn6=1.0f/dn6;\
      REALFB_KERNEL_CORE_ISR2ZERO;\
      tmp*=dn6;\
      jj4+=LJ_BS_SIZE;\
      CLFB_KERNEL_CORE_ACCUM;

#define REALPOTFB_KERNEL_CORE \
      DISTANCE_PERIODIC;\
      rij=dx*dx+dy*dy+dz*dz;\
      GRAVFB_KERNEL_CORE_EPS;\
      SCALE_R2;\
      /*cutoff=rsqrtf(rij);*/\
      dn6=rij*alpha2;\
      cutoff=sqrtf(dn6);\
      tmp=BS(0,jj4+3)*erfcf(cutoff);\
      dn6=1.0f/cutoff;\
      REALFB_KERNEL_CORE_ISR2ZERO;\
      tmp*=dn6;\
      jj4+=LJ_BS_SIZE;\
      CLFB_KERNEL_CORE_ACCUMPOT;

// LJ kernel
#define LJFB_KERNEL_CORE_DISTANCE \
      DISTANCE_PERIODIC;
#if defined(MD_SIGMA_EPSION_IN_VGSTRUCT)
#define LJFB_KERNEL_CORE_SUB1 \
      halfsigmaj=BS(0,jj4+4)+halfsigmai;\
      sqrtepsilonj=BS(0,jj4+5);\
      /*rij=dx*dx+dy*dy+dz*dz;*/\
      rij=__fadd_rn(__fadd_rn(__fmul_rn(dx,dx),__fmul_rn(dy,dy)),__fmul_rn(dz,dz));\
      SCALE_R2;\
      rij*=halfsigmaj*halfsigmaj;
#define LJFB_KERNEL_CORE_SUB3 \
      dn6=cutoff*cutoff*cutoff;\
      tmp=sqrtepsiloni*sqrtepsilonj*dn6*cutoff*(2.0f*dn6-1.0f);
#define LJFB_KERNEL_CORE_SUB3POT \
      dn6=cutoff*cutoff*cutoff;\
      tmp=sqrtepsiloni*sqrtepsilonj*dn6*(dn6-1.0f);
#elif defined(MD_MATRIX_IN_SCALER) || defined(MD_MATRIX_COPY_TO_SHARED)
#define LJFB_KERNEL_CORE_SUB1 \
      itype=__float_as_int(BS(0,jj4+3));\
      itype+=atypei;\
      /*rij=dx*dx+dy*dy+dz*dz;*/\
      rij=__fadd_rn(__fadd_rn(__fmul_rn(dx,dx),__fmul_rn(dy,dy)),__fmul_rn(dz,dz));\
      SCALE_R2;\
      rij*=AS(2,itype);
      /*rij*=d_scalers->rscales[itype];*/
#define LJFB_KERNEL_CORE_SUB3 \
      dn6=cutoff*cutoff*cutoff;\
      tmp=AS(0,itype)*dn6*cutoff*(2.0f*dn6-1.0f);
      /*tmp=d_scalers->gscalesf[itype]*dn6*cutoff*(2.0f*dn6-1.0f);*/
      /*tmp=dn6*cutoff*(2.0f*dn6-1.0f);*/
#define LJFB_KERNEL_CORE_SUB3POT \
      dn6=cutoff*cutoff*cutoff;\
      tmp=AS(0,itype)*dn6*(dn6-1.0f);
      /*tmp=d_scalers->gscalesf[itype]*dn6*cutoff*(2.0f*dn6-1.0f);*/
      /*tmp=dn6*cutoff*(2.0f*dn6-1.0f);*/
#elif defined(MD_LJ_NOATYPE)
#define LJFB_KERNEL_CORE_SUB1 \
      rij=dx*dx+dy*dy+dz*dz;
#define LJFB_KERNEL_CORE_SUB12 \
      rij=dx*dx+dy*dy+dz*dz;
#define LJFB_KERNEL_CORE_SUB3 \
      dn6=cutoff*cutoff*cutoff;\
      tmp=dn6*cutoff*(2.0f*dn6-1.0f);
#define LJFB_KERNEL_CORE_SUB3POT \
      dn6=cutoff*cutoff*cutoff;\
      tmp=dn6*(dn6-1.0f);
#elif MD_LJ_05SIGMA_8SIGMA==2 // rscale is multiplyed after cutoff
#define LJFB_KERNEL_CORE_SUB1 \
      /*itype=__float_as_int(BS(0,jj4+3));*/	\
      itype=ATYPE_UNPACK_FROM_QATYPE(__float_as_int(BS(0,jj4+3)));	\
      /*itype&=MASK(MD_QAUNION_ATYPEBIT);*/\
      itype+=atypei;\
      /*itype+=atypei & (MASK(MD_QAUNION_ATYPEBIT)<<MASK(MD_QAUNION_ATYPEBIT));*/\
      /*rij=dx*dx+dy*dy+dz*dz;*/\
      rij=__fadd_rn(__fadd_rn(__fmul_rn(dx,dx),__fmul_rn(dy,dy)),__fmul_rn(dz,dz));\
      SCALE_R2;\
      /*rij*=d_matrix[itype].rscale;*/
#define LJFB_KERNEL_CORE_SUB12 \
      itype=ATYPE_UNPACK_FROM_QATYPE(__float_as_int(BS(0,jj4+3)));\
      /*itype+=atypei;*/\
      /*rij=dx*dx+dy*dy+dz*dz;*/\
      rij=__fadd_rn(__fadd_rn(__fmul_rn(dx,dx),__fmul_rn(dy,dy)),__fmul_rn(dz,dz));\
      SCALE_R2;\
      /*rij*=d_matrix[itype].rscale;*/
#if VDW_SHIFT>=1
#if VDW_SHIFT==1
#define LJFB_KERNEL_CORE_SUB3 \
      tmp=vdw_shift<1>(rij,r2max,rcut21,cutoff,d_matrix[itype].shift,d_matrix[itype].gscale);// \
      dn6=cutoff*cutoff*cutoff;\
      tmp=cutoff;\
      cutoff=d_matrix[itype].shift;\
      dn2=rij;\
      rij*=rcut21;\
      rij*=rij*rij;\
      /* tmp*=__fmul_rn(dn6,__fadd_rn(2.0f*dn6,-1.0f))-__fmul_rn(rij*cutoff,__fadd_rn(2.0f*cutoff,-1.0f));*/ /* not correct until 120706 */ \
      tmp*=__fmul_rn(dn6,__fadd_rn(2.0f*dn6,-1.0f))-__fmul_rn(cutoff,__fadd_rn(3.0f*cutoff,-2.0f)); \
      if(dn2<MD_LJ_R2MIN)  tmp=0.0f;\
      if(dn2>=r2max)       tmp=0.0f;\
      tmp*=d_matrix[itype].gscale;
#define LJFB_KERNEL_CORE_SUB3POT \
      tmp=vdw_shift<101>(rij,r2max,rcut21,cutoff,d_matrix[itype].shift,d_matrix[itype].gscale);// \
      dn6=cutoff*cutoff*cutoff;\
      tmp=__fmul_rn(dn6,__fadd_rn(dn6,-1.0f));\
      cutoff=d_matrix[itype].shift;\
      dn2=rij;\
      rij*=rcut21;\
      rij*=rij*rij;\
      tmp+=rij*cutoff*__fadd_rn(2.0f*cutoff,-1.0f)-cutoff*__fadd_rn(3.0f*cutoff,-2.0f);\
      if(dn2<MD_LJ_R2MIN)  tmp=0.0f;					\
      if(dn2>=r2max)       tmp=0.0f;					\
      tmp*=d_matrix[itype].gscale;
#elif VDW_SHIFT==2
#define LJFB_KERNEL_CORE_SUB3 \
      tmp=vdw_shift<2>(rij,r2max,rcut21,cutoff,d_matrix[itype].shift,d_matrix[itype].gscale);// \
      dn6=cutoff*cutoff*cutoff;\
      tmp=cutoff*__fmul_rn(dn6,__fadd_rn(2.0f*dn6,-1.0f));\
      dn2=rij;\
      if(dn2<MD_LJ_R2MIN)  tmp=0.0f;\
      if(dn2>=r2max)       tmp=0.0f;\
      tmp*=d_matrix[itype].gscale;\
      if(dn2>=d_matrix[itype].shift && dn2<r2max){\
        tmp = r2max - dn2;\
        tmp /= (r2max - d_matrix[itype].shift) * (r2max - d_matrix[itype].shift) * (r2max - d_matrix[itype].shift);\
	rij = tmp * (r2max - dn2) * (r2max - 3.0f * d_matrix[itype].shift + 2.0f * dn2);\
        tmp = d_matrix[itype].gscale * ( dn6 * (dn6 - 1.0f) * 12.0f * (d_matrix[itype].shift - dn2) * tmp \
			 - 6.0f * dn6 * (dn6 + (dn6 - 1.0f) * rij) * rij/dn2 );\
      }
#define LJFB_KERNEL_CORE_SUB3POT \
      tmp=vdw_shift<102>(rij,r2max,rcut21,cutoff,d_matrix[itype].shift,d_matrix[itype].gscale);// \
      dn6=cutoff*cutoff*cutoff;\
      tmp=__fmul_rn(dn6,__fadd_rn(dn6,-1.0f));\
      dn2=rij;\
      if(dn2<MD_LJ_R2MIN)  tmp=0.0f;					\
      if(dn2>=r2max)       tmp=0.0f;					\
      tmp*=d_matrix[itype].gscale;\
      if(dn2>=d_matrix[itype].shift && dn2<r2max){\
	rij = r2max - dn2;\
        rij = rij * rij * (r2max - 3.0f * d_matrix[itype].shift + 2.0f * dn2);\
        rij /= (r2max - d_matrix[itype].shift) * (r2max - d_matrix[itype].shift) * (r2max - d_matrix[itype].shift);\
        tmp*= rij;\
      }
#endif
#else // else of VDW_SHIFT
#define LJFB_KERNEL_CORE_SUB3 \
      dn6=cutoff*cutoff*cutoff;\
      tmp=d_matrix[itype].gscale*dn6*cutoff*(__fmul_rn(2.0f,dn6)-1.0f);
#define LJFB_KERNEL_CORE_SUB3POT \
      dn6=cutoff*cutoff*cutoff;\
      tmp=d_matrix[itype].gscale*dn6*(dn6-1.0f);\
      /*tmp=rij*LOWER_VAL_FACTOR_ORG;*/\
      /*if(rij<MD_LJ_R2MIN)  tmp=0.0f;*/\
      /*if(rij>=r2max)       tmp=0.0f;*/
#endif // end of VDW_SHIFT
#else  // else of MD_LJ_05SIGMA_8SIGMA==2
#define LJFB_KERNEL_CORE_SUB1 \
      itype=__float_as_int(BS(0,jj4+3));\
      itype+=atypei;\
      /*rij=dx*dx+dy*dy+dz*dz;*/\
      rij=__fadd_rn(__fadd_rn(__fmul_rn(dx,dx),__fmul_rn(dy,dy)),__fmul_rn(dz,dz));\
      SCALE_R2;\
      rij*=d_matrix[itype].rscale;
#define LJFB_KERNEL_CORE_SUB12 \
      itype=__float_as_int(BS(0,jj4+3));\
      /*itype+=atypei;*/\
      /*rij=dx*dx+dy*dy+dz*dz;*/\
      rij=__fadd_rn(__fadd_rn(__fmul_rn(dx,dx),__fmul_rn(dy,dy)),__fmul_rn(dz,dz));\
      SCALE_R2;\
      rij*=d_matrix[itype].rscale;
#define LJFB_KERNEL_CORE_SUB3 \
      dn6=cutoff*cutoff*cutoff;\
      tmp=d_matrix[itype].gscale*dn6*cutoff*(2.0f*dn6-1.0f);\
      /*tmp=d_matrix[itype].gscale*dn6*cutoff*(__fmul_rn(2.0f,dn6)-1.0f);*/
#define LJFB_KERNEL_CORE_SUB3POT \
      dn6=cutoff*cutoff*cutoff;\
      tmp=d_matrix[itype].gscale*dn6*(dn6-1.0f);
#endif
#if MD_LJ_05SIGMA_8SIGMA==1       // active range is [0.5 sigma, 8 sigma)
#define LJFB_KERNEL_CORE_ISR2ZERO \
      /*if(rij>=MD_LJ_R2MIN) cutoff=1.0f/rij;*/\
      /*else                 cutoff=0.0f;*/\
      cutoff=1.0f/rij;\
      if(rij<MD_LJ_R2MIN)  cutoff=0.0f;\
      if(rij>=r2max)       cutoff=0.0f;
#elif MD_LJ_05SIGMA_8SIGMA==0
#define LJFB_KERNEL_CORE_ISR2ZERO \
      /*if(rij!=0.0f) cutoff=1.0f/rij;*/\
      /*else          cutoff=0.0f;*/\
      cutoff=1.0f/rij;\
      if(rij==0.0f) cutoff=0.0f;
#elif MD_LJ_05SIGMA_8SIGMA==2
#if VDW_SHIFT>=1 // same definition is used currently
#define LJFB_KERNEL_CORE_ISR2ZERO \
      cutoff=1.0f/(rij*d_matrix[itype].rscale);\
      if(rij<MD_LJ_R2MIN)  cutoff=0.0f;\
      if(rij>=r2max)       cutoff=0.0f;
#else
#define LJFB_KERNEL_CORE_ISR2ZERO \
      cutoff=1.0f/(rij*d_matrix[itype].rscale);\
      if(rij<MD_LJ_R2MIN)  cutoff=0.0f;\
      if(rij>=r2max)       cutoff=0.0f;
#endif // end of VDW_SHIFT 
#endif 

#ifdef ACCUMULATE_PAIR_FLOAT
#if ACCUMULATE_PAIR_FLOAT==1 || ACCUMULATE_PAIR_FLOAT==12 || ACCUMULATE_PAIR_FLOAT==10
#define LJFB_KERNEL_CORE_ACCUM \
      ADD_PAIR_FLOAT(dx,tmp,Csub[0],Csub2[0],tmp5,tmp2,tmp3,tmp4);\
      ADD_PAIR_FLOAT(dy,tmp,Csub[1],Csub2[1],tmp5,tmp2,tmp3,tmp4);\
      ADD_PAIR_FLOAT(dz,tmp,Csub[2],Csub2[2],tmp5,tmp2,tmp3,tmp4);
#define LJFB_KERNEL_CORE_ACCUMPOT \
      ADD_PAIR_FLOAT(1.0f,tmp,Csub[0],Csub2[0],tmp5,tmp2,tmp3,tmp4);
#define LJFB_KERNEL_CORE_ACCUMPOT2 \
      ADD_PAIR_FLOAT(1.0f,tmp,Csub[1],Csub2[1],tmp5,tmp2,tmp3,tmp4);
#elif ACCUMULATE_PAIR_FLOAT==2
#if MD_USE_LARGE_VAL==2 || MD_USE_LARGE_VAL==20
#define LJFB_KERNEL_CORE_ACCUM \
      ADD_PAIR_FLOAT_NARUMI2(dx,tmp,Csub[0],Csub2[0],tmp5,tmp2,tmp3,tmp4);\
      ADD_PAIR_FLOAT_NARUMI2(dy,tmp,Csub[1],Csub2[1],tmp5,tmp2,tmp3,tmp4);\
      ADD_PAIR_FLOAT_NARUMI2(dz,tmp,Csub[2],Csub2[2],tmp5,tmp2,tmp3,tmp4);
#define LJFB_KERNEL_CORE_ACCUMPOT \
      ADD_PAIR_FLOAT_NARUMI2(1.0f,tmp,Csub[0],Csub2[0],tmp5,tmp2,tmp3,tmp4);
#define LJFB_KERNEL_CORE_ACCUMPOT2 \
      ADD_PAIR_FLOAT_NARUMI2(1.0f,tmp,Csub[1],Csub2[1],tmp5,tmp2,tmp3,tmp4);
#else
#define LJFB_KERNEL_CORE_ACCUM \
      ADD_PAIR_FLOAT_SIMPLE(dx,tmp,Csub[0],Csub2[0],tmp5,tmp2,tmp3,tmp4);\
      ADD_PAIR_FLOAT_SIMPLE(dy,tmp,Csub[1],Csub2[1],tmp5,tmp2,tmp3,tmp4);\
      ADD_PAIR_FLOAT_SIMPLE(dz,tmp,Csub[2],Csub2[2],tmp5,tmp2,tmp3,tmp4);
#define LJFB_KERNEL_CORE_ACCUMPOT \
      ADD_PAIR_FLOAT_SIMPLE(1.0f,tmp,Csub[0],Csub2[0],tmp5,tmp2,tmp3,tmp4);
#define LJFB_KERNEL_CORE_ACCUMPOT2 \
      ADD_PAIR_FLOAT_SIMPLE(1.0f,tmp,Csub[1],Csub2[1],tmp5,tmp2,tmp3,tmp4);
#endif
#elif ACCUMULATE_PAIR_FLOAT==3
#define LJFB_KERNEL_CORE_ACCUM \
      ADD_PAIR_FLOAT_DOUBLE(dx,tmp,Csubd[0],tmp5,tmp2,tmp3,tmp4);\
      ADD_PAIR_FLOAT_DOUBLE(dy,tmp,Csubd[1],tmp5,tmp2,tmp3,tmp4);\
      ADD_PAIR_FLOAT_DOUBLE(dz,tmp,Csubd[2],tmp5,tmp2,tmp3,tmp4);
#define LJFB_KERNEL_CORE_ACCUMPOT \
      ADD_PAIR_FLOAT_DOUBLE(1.0f,tmp,Csubd[0],tmp5,tmp2,tmp3,tmp4);
#define LJFB_KERNEL_CORE_ACCUMPOT2 \
      ADD_PAIR_FLOAT_DOUBLE(1.0f,tmp,Csubd[1],tmp5,tmp2,tmp3,tmp4);
#endif
#else
#define LJFB_KERNEL_CORE_ACCUM \
      Csub[0]+=tmp*dx;\
      Csub[1]+=tmp*dy;\
      Csub[2]+=tmp*dz;
#define LJFB_KERNEL_CORE_ACCUMPOT \
      Csub[0]+=tmp;
#define LJFB_KERNEL_CORE_ACCUMPOT2 \
      Csub[1]+=tmp;
#endif
#define LJFB_KERNEL_CORE \
      LJFB_KERNEL_CORE_DISTANCE;\
      LJFB_KERNEL_CORE_SUB1;\
      LJFB_KERNEL_CORE_ISR2ZERO;\
      LJFB_KERNEL_CORE_SUB3;\
      jj4+=LJ_BS_SIZE;\
      /*Csub[0]+=rij;*/\
      /*Csub[0]+=dx*dx+dy*dy;*/\
      /*Csub[0]+=__fadd_rn(__fmul_rn(dx,dx),__fmul_rn(dy,dy));*/\
      LJFB_KERNEL_CORE_ACCUM;

#define LJPOTFB_KERNEL_CORE \
      LJFB_KERNEL_CORE_DISTANCE;\
      LJFB_KERNEL_CORE_SUB1;\
      LJFB_KERNEL_CORE_ISR2ZERO;\
      LJFB_KERNEL_CORE_SUB3POT;\
      jj4+=LJ_BS_SIZE;\
      LJFB_KERNEL_CORE_ACCUMPOT;
#ifdef MD_SORT_ATYPEI
#define LJFB_KERNEL_CORE2 \
      LJFB_KERNEL_CORE_DISTANCE;\
      LJFB_KERNEL_CORE_SUB12;\
      LJFB_KERNEL_CORE_ISR2ZERO;\
      LJFB_KERNEL_CORE_SUB3;\
      jj4+=LJ_BS_SIZE;\
      LJFB_KERNEL_CORE_ACCUM;
#define LJPOTFB_KERNEL_CORE2 \
      LJFB_KERNEL_CORE_DISTANCE;\
      LJFB_KERNEL_CORE_SUB12;\
      LJFB_KERNEL_CORE_ISR2ZERO;\
      LJFB_KERNEL_CORE_SUB3POT;\
      jj4+=LJ_BS_SIZE;\
      LJFB_KERNEL_CORE_ACCUMPOT;
#endif

#define REALLJFB_KERNEL_CORE_SUB3 \
      dn2=rij;tmp=vdw_shift<VDW_SHIFT>(rij,r2max,rcut21,cutoff,d_matrix[itype].shift,d_matrix[itype].gscale); // \
      dn6=cutoff*cutoff*cutoff;						\
      tmp=cutoff;								\
      cutoff=d_matrix[itype].shift;						\
      if(rij<MD_LJ_R2MIN)  cutoff=0.0f;					\
      if(rij>=r2max)       cutoff=0.0f;					\
      dn2=rij;								\
      rij*=rcut21;								\
      rij*=rij*rij;								\
      tmp*=__fmul_rn(dn6,__fadd_rn(2.0f*dn6,-1.0f))-__fmul_rn(rij*cutoff,__fadd_rn(2.0f*cutoff,-1.0f)); \
      tmp*=d_matrix[itype].gscale;						\
      /*if(dn2!=0.0f && sqrtf(dn2)<1.9f) tmp=1.0f; else tmp=0.0f;*/
#define REALLJFB_KERNEL_CORE_SUB3POT \
      dn2=rij;tmp=vdw_shift<VDW_SHIFT+100>(rij,r2max,rcut21,cutoff,d_matrix[itype].shift,d_matrix[itype].gscale); // \
      dn6=cutoff*cutoff*cutoff;\
      tmp=__fmul_rn(dn6,__fadd_rn(dn6,-1.0f));\
      cutoff=d_matrix[itype].shift;\
      dn2=rij;\
      rij*=rcut21;\
      rij*=rij*rij;\
      tmp+=rij*cutoff*__fadd_rn(2.0f*cutoff,-1.0f)-cutoff*__fadd_rn(3.0f*cutoff,-2.0f);\
      if(dn2<MD_LJ_R2MIN)  tmp=0.0f;\
      if(dn2>=r2max)       tmp=0.0f;\
      tmp*=d_matrix[itype].gscale;\
      /*if(dn2>=MD_LJ_R2MIN && dn2<r2max) tmp=dn2;*/              \
      /*if(dn2!=0.0f && sqrtf(dn2)<1.0f) tmp=tmp; else tmp=0.0f;*/      \
      /*if(dn2!=0.0f && sqrtf(dn2)<1.9f) tmp=1.0f; else tmp=0.0f;*/
#define REALLJFB_KERNEL_CORE_SUB4 \
      {\
        float sqdn;\
        float qj;\
        float alpha2=(d_scalers->alpha)*(d_scalers->alpha);\
	float api2=(float)(2.e0/sqrt(M_PI));				\
        qj=Q_INT_TO_FLOAT(__float_as_int(BS(0,jj4+3)),d_scalers->scaleqj_1); \
        rij=dn2;\
        dn2*=alpha2;\
        sqdn=sqrtf(dn2);\
	tmp = qj*(api2*expf(-dn2)+erfcf(sqdn)/sqdn);\
        dn2=1.0f/dn2;\
        tmp*=dn2;\
        tmp*=LOWER_VAL_FACTOR_ORG;\
        if(rij<MD_LJ_R2MIN)  tmp=0.0f;	\
        if(rij>=r2max)       tmp=0.0f;	\
        tmp*=Csub[3]*d_scalers->alpha3fac;	\
        /*Csub[0]+=tmp*dx;Csub[1]+=tmp*dy;Csub[2]+=tmp*dz;*/	\
      }
#define REALLJFB_KERNEL_CORE_SUB4POT \
      {\
        float sqdn;\
        float qj;\
        float alpha2=(d_scalers->alpha)*(d_scalers->alpha);\
	float api2=(float)(2.e0/sqrt(M_PI));				\
        qj=Q_INT_TO_FLOAT(__float_as_int(BS(0,jj4+3)),d_scalers->scaleqj_1); \
        rij=dn2;\
        dn2*=alpha2;\
        sqdn=sqrtf(dn2);\
	tmp = qj*erfcf(sqdn);\
        dn2=1.0f/sqdn;\
        tmp*=dn2;\
        tmp*=LOWER_VAL_FACTOR_ORG;\
        if(rij<MD_LJ_R2MIN)  tmp=0.0f;	\
        if(rij>=r2max)       tmp=0.0f;	\
        tmp*=Csub[3]*d_scalers->alphafac;			\
      }

#define MYERFCF(z,t,r) /* z:input, t:tmporary, r:return */\
	t=1.0f/(1.0f+0.5f*z);\
        r=t*expf(-z*z-1.26551223f+t*(1.00002368f+t*(0.37409196f+t*(0.09678418f+\
             t*(-0.18628806f+t*(0.27886807f+t*(-1.13520398f+t*(1.48851587f+\
		t*(-0.82215223f+t*0.17087277f)))))))));

#define REALLJFB_KERNEL_CORE_SNDHALF_DEBUG \
      if(rij>=MD_LJ_R2MIN && rij<r2max){\
        cutoff=1.0f/(rij*d_matrix[itype].rscale);\
        dn6=cutoff*cutoff*cutoff;\
        tmp=d_matrix[itype].gscale*cutoff*dn6*(2.0f-dn6);\
        LJFB_KERNEL_CORE_ACCUM;\
        \
        qj=Q_INT_TO_FLOAT(__float_as_int(BS(0,jj4+3)),scaleqj_1); \
        tmp = LOWER_VAL_FACTOR_ORG*Csub[3]*alpha3fac*qj/sqrtf(rij);\
        LJFB_KERNEL_CORE_ACCUM;\
      }\
      jj4+=LJ_BS_SIZE;

#define REALLJFB_KERNEL_CORE_SNDHALF_DEBUG2 \
      if(rij>=MD_LJ_R2MIN && rij<r2max){\
        cutoff=1.0f/(rij*d_matrix[itype].rscale);\
        dn6=cutoff*cutoff*cutoff;\
        tmp=d_matrix[itype].gscale*cutoff*dn6*(2.0f-dn6);\
        LJFB_KERNEL_CORE_ACCUM;\
        \
        qj=Q_INT_TO_FLOAT(__float_as_int(BS(0,jj4+3)),scaleqj_1); \
        MYERFCF(sqdn,tmp3,tmp2);\
        tmp = qj*(api2*expf(-dn2)+tmp2/sqdn);\
        LJFB_KERNEL_CORE_ACCUM;\
      }\
      jj4+=LJ_BS_SIZE;

#define REALLJFB_KERNEL_CORE_SNDHALF_DEBUG3 \
      if(rij>=MD_LJ_R2MIN && rij<r2max){\
        cutoff=1.0f;\
        dn6=cutoff*cutoff*cutoff;\
        tmp=cutoff*dn6*(2.0f-dn6);\
        LJFB_KERNEL_CORE_ACCUM;\
        \
        qj=Q_INT_TO_FLOAT(__float_as_int(BS(0,jj4+3)),scaleqj_1); \
        MYERFCF(sqdn,tmp3,tmp2);\
        tmp = qj*(api2*expf(-dn2)+tmp2/sqdn);\
        LJFB_KERNEL_CORE_ACCUM;\
      }\
      jj4+=LJ_BS_SIZE;

#define REALLJFB_KERNEL_CORE_SNDHALF_DEBUG4 \
      if(rij>=MD_LJ_R2MIN && rij<r2max){\
/*      if(rij!=0.0f){*/\
        cutoff=1.0f/(rij*d_matrix[itype].rscale);\
        dn6=cutoff*cutoff*cutoff;\
        tmp=d_matrix[itype].gscale*cutoff*dn6*(2.0f-dn6);\
/*        cutoff=rij;\
        dn6=cutoff*cutoff*cutoff;\
        tmp=cutoff*dn6*(2.0f-dn6);*/\
      Csub[0]+=tmp*dx;\
      Csub[1]+=tmp*dy;\
      Csub[2]+=tmp*dz;\
 /*        LJFB_KERNEL_CORE_ACCUM;*/\
      }\
      jj4+=LJ_BS_SIZE;

#define REALLJFB_KERNEL_CORE_SNDHALF_WORK3 \
      if(rij>=MD_LJ_R2MIN && rij<r2max){\
	tmp=vdw_shift<VDW_SHIFT>(rij,r2max,rcut21,1.0f/(rij*d_matrix[itype].rscale),\
                         d_matrix[itype].shift,d_matrix[itype].gscale);	  \
        /* cutoff=1.0f/(rij*d_matrix[itype].rscale); */			\
        /* dn6=cutoff*cutoff*cutoff; */					\
        /* tmp=cutoff; */						\
        /* cutoff=d_matrix[itype].shift; */				\
        /* dn2=rij; */							\
        /* dn2*=rcut21; */						\
        /* dn2*=dn2*dn2; */						\
        /* tmp*=__fmul_rn(dn6,__fadd_rn(2.0f*dn6,-1.0f))-__fmul_rn(dn2*cutoff,__fadd_rn(2.0f*cutoff,-1.0f)); */ \
        /* tmp*=d_matrix[itype].gscale; */			\
        LJFB_KERNEL_CORE_ACCUM;\
        \
        qj=Q_INT_TO_FLOAT(__float_as_int(BS(0,jj4+3)),scaleqj_1); \
        /*tmp = LOWER_VAL_FACTOR_ORG*Csub[3]*alpha3fac*qj/sqrtf(rij);*/\
        dn2=rij*alpha2;\
        sqdn=sqrtf(dn2);\
        /*tmp = qj*(api2*expf(-dn2)+erfcf(sqdn)/sqdn);*/\
        MYERFCF(sqdn,tmp3,tmp2);\
        tmp = qj*(api2*expf(-dn2)+tmp2/sqdn);\
        dn2=1.0f/dn2;\
        tmp*=dn2;\
        tmp*=LOWER_VAL_FACTOR_ORG;\
        tmp*=Csub[3]*alpha3fac;	\
        LJFB_KERNEL_CORE_ACCUM;\
      }\
      jj4+=LJ_BS_SIZE;

#define REALLJFB_KERNEL_CORE_SNDHALF_WORK2 \
      /* this part is not modified to use vdw_shift because this routine is not used */ \
      if(rij>=MD_LJ_R2MIN && rij<r2max){\
        cutoff=1.0f/(rij*d_matrix[itype].rscale);\
        dn6=cutoff*cutoff*cutoff;\
        tmp=cutoff;\
        cutoff=d_matrix[itype].shift;\
        dn2=rij;\
        dn2*=rcut21;\
        dn2*=dn2*dn2;\
        tmp*=__fmul_rn(dn6,__fadd_rn(2.0f*dn6,-1.0f))-__fmul_rn(dn2*cutoff,__fadd_rn(2.0f*cutoff,-1.0f));\
        tmp*=d_matrix[itype].gscale;\
        \
        LJFB_KERNEL_CORE_ACCUM;\
        \
        qj=Q_INT_TO_FLOAT(__float_as_int(BS(0,jj4+3)),scaleqj_1); \
        dn2=rij*alpha2;\
        sqdn=sqrtf(dn2);\
        tmp = qj*(api2*expf(-dn2)+erfcf(sqdn)/sqdn);\
        dn2=1.0f/dn2;\
        tmp*=dn2;\
        tmp*=LOWER_VAL_FACTOR_ORG;\
        tmp*=Csub[3]*alpha3fac;	\
        LJFB_KERNEL_CORE_ACCUM;\
      }\
      jj4+=LJ_BS_SIZE;

#define REALLJFB_KERNEL_CORE_SNDHALF_WORK \
      LJFB_KERNEL_CORE_ISR2ZERO;\
      REALLJFB_KERNEL_CORE_SUB3;\
      LJFB_KERNEL_CORE_ACCUM;		\
      REALLJFB_KERNEL_CORE_SUB4;\
      jj4+=LJ_BS_SIZE;\
      LJFB_KERNEL_CORE_ACCUM;

#define REALLJFB_KERNEL_CORE_SNDHALF \
        REALLJFB_KERNEL_CORE_SNDHALF_WORK3
//        REALLJFB_KERNEL_CORE_SNDHALF_DEBUG4

#if VG_JDIV==1 || 0
#define REALLJFB_KERNEL_CORE \
      LJFB_KERNEL_CORE_DISTANCE;\
      LJFB_KERNEL_CORE_SUB1;\
      REALLJFB_KERNEL_CORE_SNDHALF;\
      /*LJFB_KERNEL_CORE_ISR2ZERO;*/\
      /*REALLJFB_KERNEL_CORE_SUB3;*/\
      /*LJFB_KERNEL_CORE_ACCUM;*/		\
      /*REALLJFB_KERNEL_CORE_SUB4;*/\
      /*jj4+=LJ_BS_SIZE;*/\
      /*LJFB_KERNEL_CORE_ACCUM;*/
#else // else of VG_JDIV==1
#define REALLJFB_KERNEL_CORE \
      /* LJFB_KERNEL_CORE_DISTANCE */\
      /*iii=(ty/VG_JDIV*LJ_BS_SIZE)+ii4;*/\
      jjj=jj4+(ty/(VG_MINIMUM_PARTICLE_BLOCK_I/VG_JDIV))*LJ_BS_SIZE;\
      /*jjj=jj4;*/\
/*      iii=ty4;\
      atypei=ATYPE_UNPACK_FROM_QATYPE(__float_as_int(Xis[iii+3]));\
      atypei*=natj;\
      Csub[3]=Q_INT_TO_FLOAT(__float_as_int(Xis[iii+3]),d_scalers->scaleqi_1);\
      dx = Xis[iii+0] - BS(0,jjj+0); \
      dy = Xis[iii+1] - BS(0,jjj+1); \
      dz = Xis[iii+2] - BS(0,jjj+2); */\
      dx = xi[0] - BS(0,jjj+0); \
      dy = xi[1] - BS(0,jjj+1); \
      dz = xi[2] - BS(0,jjj+2); \
      dx = dx - rintf(dx * al2[0]) * size[0];\
      dy = dy - rintf(dy * al2[1]) * size[1];\
      dz = dz - rintf(dz * al2[2]) * size[2];\
      /*LJFB_KERNEL_CORE_SUB1;		     */	\
      itype=ATYPE_UNPACK_FROM_QATYPE(__float_as_int(BS(0,jjj+3)));	\
      itype+=atypei;\
      rij=__fadd_rn(__fadd_rn(__fmul_rn(dx,dx),__fmul_rn(dy,dy)),__fmul_rn(dz,dz));\
      /*REALLJFB_KERNEL_CORE_SNDHALF;*/\
      if(rij>=MD_LJ_R2MIN && rij<r2max){\
	tmp=vdw_shift<VDW_SHIFT>(rij,r2max,rcut21,1.0f/(rij*d_matrix[itype].rscale),  \
                          d_matrix[itype].shift,d_matrix[itype].gscale);  \
        /* cutoff=1.0f/(rij*d_matrix[itype].rscale); */			\
        /* dn6=cutoff*cutoff*cutoff; */					\
        /* tmp=cutoff; */						\
        /* cutoff=d_matrix[itype].shift; */				\
        /* dn2=rij; */							\
        /* dn2*=rcut21;	*/						\
        /* dn2*=dn2*dn2; */						\
        /* tmp*=__fmul_rn(dn6,__fadd_rn(2.0f*dn6,-1.0f))-__fmul_rn(dn2*cutoff,__fadd_rn(2.0f*cutoff,-1.0f)); */ \
        /* tmp*=d_matrix[itype].gscale;	*/				\
        LJFB_KERNEL_CORE_ACCUM;\
        \
        qj=Q_INT_TO_FLOAT(__float_as_int(BS(0,jjj+3)),scaleqj_1); \
        dn2=rij*alpha2;\
        sqdn=sqrtf(dn2);\
        MYERFCF(sqdn,tmp3,tmp2);\
        tmp = qj*(api2*expf(-dn2)+tmp2/sqdn);\
        dn2=1.0f/dn2;\
        tmp*=dn2;\
        tmp*=LOWER_VAL_FACTOR_ORG;\
        tmp*=Csub[3]*alpha3fac;	\
        LJFB_KERNEL_CORE_ACCUM;\
      }\
      jj4+=VG_JDIV*LJ_BS_SIZE;\
      /*jj4+=LJ_BS_SIZE;*/

#endif // end of VG_JDIV==1

// for real and lj potential (currently tblno=47)
#if VG_JDIV==1 || 0
#define REALLJPOTFB_KERNEL_CORE \
      LJFB_KERNEL_CORE_DISTANCE;\
      LJFB_KERNEL_CORE_SUB1;\
      REALLJPOTFB_KERNEL_CORE_SNDHALF; /* not defined yet */    \
      /*LJFB_KERNEL_CORE_ISR2ZERO;*/\
      /*REALLJFB_KERNEL_CORE_SUB3;*/\
      /*LJFB_KERNEL_CORE_ACCUM;*/		\
      /*REALLJFB_KERNEL_CORE_SUB4;*/\
      /*jj4+=LJ_BS_SIZE;*/\
      /*LJFB_KERNEL_CORE_ACCUM;*/
#else // else of VG_JDIV==1
#define REALLJPOTFB_KERNEL_CORE \
      jjj=jj4+(ty/(VG_MINIMUM_PARTICLE_BLOCK_I/VG_JDIV))*LJ_BS_SIZE;\
      dx = xi[0] - BS(0,jjj+0); \
      dy = xi[1] - BS(0,jjj+1); \
      dz = xi[2] - BS(0,jjj+2); \
      dx = dx - rintf(dx * al2[0]) * size[0];\
      dy = dy - rintf(dy * al2[1]) * size[1];\
      dz = dz - rintf(dz * al2[2]) * size[2];\
      itype=ATYPE_UNPACK_FROM_QATYPE(__float_as_int(BS(0,jjj+3)));	\
      itype+=atypei;\
      rij=__fadd_rn(__fadd_rn(__fmul_rn(dx,dx),__fmul_rn(dy,dy)),__fmul_rn(dz,dz));\
      /*REALLJPOTFB_KERNEL_CORE_SNDHALF;*/\
      if(rij>=MD_LJ_R2MIN && rij<r2max){\
	tmp=vdw_shift<VDW_SHIFT+100>(rij,r2max,rcut21,1.0f/(rij*d_matrix[itype].rscale),\
                         d_matrix[itype].shift,d_matrix[itype].gscale);	  \
        /* cutoff=1.0f/(rij*d_matrix[itype].rscale); */			\
        /* dn6=cutoff*cutoff*cutoff; */					\
        /* tmp=__fmul_rn(dn6,__fadd_rn(dn6,-1.0f)); */			\
        /* cutoff=d_matrix[itype].shift; */				\
        /* dn2=rij; */							\
        /* dn2*=rcut21; */						\
        /* dn2*=dn2*dn2; */						\
        /* tmp+=__fmul_rn(dn2,cutoff*__fadd_rn(2.0f*cutoff,-1.0f))-cutoff*__fadd_rn(3.0f*cutoff,-2.0f); */ \
        /* tmp*=d_matrix[itype].gscale; */				\
        LJFB_KERNEL_CORE_ACCUMPOT2;\
        \
        qj=Q_INT_TO_FLOAT(__float_as_int(BS(0,jjj+3)),scaleqj_1); \
        dn2=rij*alpha2;\
        sqdn=sqrtf(dn2);\
        MYERFCF(sqdn,tmp3,tmp2);\
	tmp = qj*tmp2;\
        dn2=1.0f/sqdn;\
        tmp*=dn2;\
        tmp*=LOWER_VAL_FACTOR_ORG;\
        tmp*=Csub[3]*d_scalers->alphafac;			\
        LJFB_KERNEL_CORE_ACCUMPOT;\
      }\
      jj4+=VG_JDIV*LJ_BS_SIZE;\
      /*jj4+=LJ_BS_SIZE;*/

#endif // end of VG_JDIV==1



#define REALLJPOTFB_NOJDIV_KERNEL_CORE \
      LJFB_KERNEL_CORE_DISTANCE;\
      LJFB_KERNEL_CORE_SUB1;\
      LJFB_KERNEL_CORE_ISR2ZERO;\
      REALLJFB_KERNEL_CORE_SUB3POT;\
      LJFB_KERNEL_CORE_ACCUMPOT2;		\
      REALLJFB_KERNEL_CORE_SUB4POT;\
      jj4+=LJ_BS_SIZE;\
      LJFB_KERNEL_CORE_ACCUMPOT;

#define REALLJFB_KERNEL_CORE_TMP \
      dn2 = dx * dx + dy * dy + dz * dz;\
      if(j*VG_MINIMUM_PARTICLE_BLOCK_I+jj<nj && dn2!=0.0f){\
	float sqdn;\
	float qj;\
        float r2=dn2;\
        float alpha2=(d_scalers->alpha)*(d_scalers->alpha);\
	float api2=(float)(2.e0/sqrt(M_PI));				\
        dn2*=alpha2;\
        sqdn=sqrtf(dn2);\
        qj=Q_INT_TO_FLOAT(__float_as_int(BS(0,jj4+3)),d_scalers->scaleqj_1); \
	tmp = qj*(api2*expf(-dn2)+erfcf(sqdn)/sqdn);\
        dn2=1.0f/dn2;\
        REAL_KERNEL_CORE_ISR2ZERO;\
        tmp*=dn2;\
	Csub[0] += tmp * dx;\
	Csub[1] += tmp * dy;\
	Csub[2] += tmp * dz;\
      }

/*
#define KERNEL_CORE_COULOMB        GRAVFB_KERNEL_CORE
#define KERNEL_CORE_LJ             LJFB_KERNEL_CORE
#define KERNEL_CORE_LJ_SAMEATYPEI  LJFB_KERNEL_CORE2

#define KERNEL_CORE_COULOMB4 \
      KERNEL_CORE_COULOMB;\
      KERNEL_CORE_COULOMB;\
      KERNEL_CORE_COULOMB;\
      KERNEL_CORE_COULOMB;
#define KERNEL_CORE_COULOMB8 \
      KERNEL_CORE_COULOMB4;\
      KERNEL_CORE_COULOMB4;
#define KERNEL_CORE_COULOMB16 \
      KERNEL_CORE_COULOMB8;\
      KERNEL_CORE_COULOMB8;
#define KERNEL_CORE_COULOMB32 \
      KERNEL_CORE_COULOMB16;\
      KERNEL_CORE_COULOMB16;
#define KERNEL_CORE_COULOMB64 \
      KERNEL_CORE_COULOMB32;\
      KERNEL_CORE_COULOMB32;
#define KERNEL_CORE_COULOMB96 \
      KERNEL_CORE_COULOMB64;\
      KERNEL_CORE_COULOMB32;
#define KERNEL_CORE_COULOMB128 \
      KERNEL_CORE_COULOMB64;\
      KERNEL_CORE_COULOMB64;
#define KERNEL_CORE_COULOMB192 \
      KERNEL_CORE_COULOMB128;\
      KERNEL_CORE_COULOMB64;
#define KERNEL_CORE_COULOMB256 \
      KERNEL_CORE_COULOMB128;\
      KERNEL_CORE_COULOMB128;

#ifdef MD_UNROLL_KERNEL
#if MD_UNROLL_KERNEL==4
#define KERNEL_CORE_COULOMB_UNROLL KERNEL_CORE_COULOMB4
#elif MD_UNROLL_KERNEL==8
#define KERNEL_CORE_COULOMB_UNROLL KERNEL_CORE_COULOMB8
#elif MD_UNROLL_KERNEL==16
#define KERNEL_CORE_COULOMB_UNROLL KERNEL_CORE_COULOMB16
#elif MD_UNROLL_KERNEL==32
#define KERNEL_CORE_COULOMB_UNROLL KERNEL_CORE_COULOMB32
#elif MD_UNROLL_KERNEL==64
#define KERNEL_CORE_COULOMB_UNROLL KERNEL_CORE_COULOMB64
#elif MD_UNROLL_KERNEL==128
#define KERNEL_CORE_COULOMB_UNROLL KERNEL_CORE_COULOMB128
#endif
#else // else of MD_UNROLL_KERNEL
#if VG_MINIMUM_PARTICLE_BLOCK_I==64
#define KERNEL_CORE_COULOMB_UNROLL KERNEL_CORE_COULOMB64
#elif VG_MINIMUM_PARTICLE_BLOCK_I==96
#define KERNEL_CORE_COULOMB_UNROLL KERNEL_CORE_COULOMB96
#elif VG_MINIMUM_PARTICLE_BLOCK_I==128
#define KERNEL_CORE_COULOMB_UNROLL KERNEL_CORE_COULOMB128
#elif VG_MINIMUM_PARTICLE_BLOCK_I==192
#define KERNEL_CORE_COULOMB_UNROLL KERNEL_CORE_COULOMB192
#elif VG_MINIMUM_PARTICLE_BLOCK_I==256
#define KERNEL_CORE_COULOMB_UNROLL KERNEL_CORE_COULOMB256
#endif
#endif // end of MD_UNROLL_KERNEL


#define KERNEL_CORE_LJ4 \
      KERNEL_CORE_LJ;\
      KERNEL_CORE_LJ;\
      KERNEL_CORE_LJ;\
      KERNEL_CORE_LJ;
#define KERNEL_CORE_LJ8 \
      KERNEL_CORE_LJ4;\
      KERNEL_CORE_LJ4;
#define KERNEL_CORE_LJ16 \
      KERNEL_CORE_LJ8;\
      KERNEL_CORE_LJ8;
#define KERNEL_CORE_LJ64 \
      KERNEL_CORE_LJ16;\
      KERNEL_CORE_LJ16;\
      KERNEL_CORE_LJ16;\
      KERNEL_CORE_LJ16;
#define KERNEL_CORE_LJ96 \
      KERNEL_CORE_LJ64;\
      KERNEL_CORE_LJ16;\
      KERNEL_CORE_LJ16;
#define KERNEL_CORE_LJ128 \
      KERNEL_CORE_LJ64;\
      KERNEL_CORE_LJ64;
#define KERNEL_CORE_LJ192 \
      KERNEL_CORE_LJ128;\
      KERNEL_CORE_LJ64;
#define KERNEL_CORE_LJ256 \
      KERNEL_CORE_LJ128;\
      KERNEL_CORE_LJ128;

#ifdef MD_UNROLL_KERNEL
#if MD_UNROLL_KERNEL==4
#define KERNEL_CORE_LJ_UNROLL KERNEL_CORE_LJ4
#elif MD_UNROLL_KERNEL==8
#define KERNEL_CORE_LJ_UNROLL KERNEL_CORE_LJ8
#elif MD_UNROLL_KERNEL==16
#define KERNEL_CORE_LJ_UNROLL KERNEL_CORE_LJ16
#elif MD_UNROLL_KERNEL==32
#define KERNEL_CORE_LJ_UNROLL KERNEL_CORE_LJ32
#elif MD_UNROLL_KERNEL==64
#define KERNEL_CORE_LJ_UNROLL KERNEL_CORE_LJ64
#elif MD_UNROLL_KERNEL==128
#define KERNEL_CORE_LJ_UNROLL KERNEL_CORE_LJ128
#endif
#else // else of MD_UNROLL_KERNEL
#if VG_MINIMUM_PARTICLE_BLOCK_I==64
#define KERNEL_CORE_LJ_UNROLL KERNEL_CORE_LJ64
#elif VG_MINIMUM_PARTICLE_BLOCK_I==96
#define KERNEL_CORE_LJ_UNROLL KERNEL_CORE_LJ96
#elif VG_MINIMUM_PARTICLE_BLOCK_I==128
#define KERNEL_CORE_LJ_UNROLL KERNEL_CORE_LJ128
#elif VG_MINIMUM_PARTICLE_BLOCK_I==192
#define KERNEL_CORE_LJ_UNROLL KERNEL_CORE_LJ192
#elif VG_MINIMUM_PARTICLE_BLOCK_I==256
#define KERNEL_CORE_LJ_UNROLL KERNEL_CORE_LJ256
#endif
#endif // end of MD_UNROLL_KERNEL


#define KERNEL_CORE_LJ_SAMEATYPEI4 \
      KERNEL_CORE_LJ_SAMEATYPEI;\
      KERNEL_CORE_LJ_SAMEATYPEI;\
      KERNEL_CORE_LJ_SAMEATYPEI;\
      KERNEL_CORE_LJ_SAMEATYPEI;
#define KERNEL_CORE_LJ_SAMEATYPEI8 \
      KERNEL_CORE_LJ_SAMEATYPEI4;\
      KERNEL_CORE_LJ_SAMEATYPEI4;
#define KERNEL_CORE_LJ_SAMEATYPEI16 \
      KERNEL_CORE_LJ_SAMEATYPEI8;\
      KERNEL_CORE_LJ_SAMEATYPEI8;
#define KERNEL_CORE_LJ_SAMEATYPEI64 \
      KERNEL_CORE_LJ_SAMEATYPEI16;\
      KERNEL_CORE_LJ_SAMEATYPEI16;\
      KERNEL_CORE_LJ_SAMEATYPEI16;\
      KERNEL_CORE_LJ_SAMEATYPEI16;
#define KERNEL_CORE_LJ_SAMEATYPEI96 \
      KERNEL_CORE_LJ_SAMEATYPEI64;\
      KERNEL_CORE_LJ_SAMEATYPEI16;\
      KERNEL_CORE_LJ_SAMEATYPEI16;
#define KERNEL_CORE_LJ_SAMEATYPEI128 \
      KERNEL_CORE_LJ_SAMEATYPEI64;\
      KERNEL_CORE_LJ_SAMEATYPEI64;
#define KERNEL_CORE_LJ_SAMEATYPEI192 \
      KERNEL_CORE_LJ_SAMEATYPEI128;\
      KERNEL_CORE_LJ_SAMEATYPEI64;
#define KERNEL_CORE_LJ_SAMEATYPEI256 \
      KERNEL_CORE_LJ_SAMEATYPEI128;\
      KERNEL_CORE_LJ_SAMEATYPEI128;

#ifdef MD_UNROLL_KERNEL
#if MD_UNROLL_KERNEL==4
#define KERNEL_CORE_LJ_SAMEATYPEI_UNROLL KERNEL_CORE_LJ_SAMEATYPEI4
#elif MD_UNROLL_KERNEL==8
#define KERNEL_CORE_LJ_SAMEATYPEI_UNROLL KERNEL_CORE_LJ_SAMEATYPEI8
#elif MD_UNROLL_KERNEL==16
#define KERNEL_CORE_LJ_SAMEATYPEI_UNROLL KERNEL_CORE_LJ_SAMEATYPEI16
#elif MD_UNROLL_KERNEL==32
#define KERNEL_CORE_LJ_SAMEATYPEI_UNROLL KERNEL_CORE_LJ_SAMEATYPEI32
#elif MD_UNROLL_KERNEL==64
#define KERNEL_CORE_LJ_SAMEATYPEI_UNROLL KERNEL_CORE_LJ_SAMEATYPEI64
#elif MD_UNROLL_KERNEL==128
#define KERNEL_CORE_LJ_SAMEATYPEI_UNROLL KERNEL_CORE_LJ_SAMEATYPEI128
#endif
#else // else of MD_UNROLL_KERNEL
#if VG_MINIMUM_PARTICLE_BLOCK_I==64
#define KERNEL_CORE_LJ_SAMEATYPEI_UNROLL KERNEL_CORE_LJ_SAMEATYPEI64
#elif VG_MINIMUM_PARTICLE_BLOCK_I==96
#define KERNEL_CORE_LJ_SAMEATYPEI_UNROLL KERNEL_CORE_LJ_SAMEATYPEI96
#elif VG_MINIMUM_PARTICLE_BLOCK_I==128
#define KERNEL_CORE_LJ_SAMEATYPEI_UNROLL KERNEL_CORE_LJ_SAMEATYPEI128
#elif VG_MINIMUM_PARTICLE_BLOCK_I==192
#define KERNEL_CORE_LJ_SAMEATYPEI_UNROLL KERNEL_CORE_LJ_SAMEATYPEI192
#elif VG_MINIMUM_PARTICLE_BLOCK_I==256
#define KERNEL_CORE_LJ_SAMEATYPEI_UNROLL KERNEL_CORE_LJ_SAMEATYPEI256
#endif
#endif // end of MD_UNROLL_KERNEL
*/

#if MD_PERIODIC_FIXED==1
#define DEFINE_AS_BS \
  __shared__ int As[3][VG_MINIMUM_PARTICLE_BLOCK_I];\
  __shared__ int Bs[3][VG_MINIMUM_PARTICLE_BLOCK_I];
#define LOAD_AS \
  for(k=0;k<3;k++) AS(k,ty)=__float_as_int(ivec[aBegin].r[k]);
#define LOAD_BS \
    for(k=0;k<3;k++) BS(k,ty)=__float_as_int(jvec[jty].r[k]);
#else
#define DEFINE_AS_BS \
  __shared__ float As[3][VG_MINIMUM_PARTICLE_BLOCK_I];\
  __shared__ float Bs[3][VG_MINIMUM_PARTICLE_BLOCK_I];
#define LOAD_AS \
  for(k=0;k<3;k++) AS(k,ty)=ivec[aBegin].r[k];
#define LOAD_BS \
    for(k=0;k<3;k++) BS(k,ty)=jvec[jty].r[k];
#endif

#ifdef MD_QAUNION_ATYPEBIT
#define LOAD_QATYPEI \
  DS(0,ty)=ATYPE_UNPACK_FROM_QATYPE(ivec[aBegin].qatype.atype);\
  QIS(0,ty)=Q_INT_TO_FLOAT(ivec[aBegin].qatype.atype,d_scalers->scaleqi_1);
#define LOAD_QATYPEJ \
  ES(0,ty)=ATYPE_UNPACK_FROM_QATYPE(jvec[jty].qatype.atype);\
  QJS(0,ty)=Q_INT_TO_FLOAT(jvec[jty].qatype.atype,d_scalers->scaleqj_1);
#else
#define LOAD_QATYPEI \
  DS(0,ty)=ivec[aBegin].qatype.atype;\
  QIS(0,ty)=ivec[aBegin].qatype.q;
#define LOAD_QATYPEJ \
  ES(0,ty)=jvec[jty].qatype.atype;\
  QJS(0,ty)=jvec[jty].qatype.q;
#endif

#if MD_PERIODIC_FIXED==1
#define DISTANCE_PERIODIC_OLD \
      dx = AS(0,ty) - BS(0,jj);\
      dy = AS(1,ty) - BS(1,jj);\
      dz = AS(2,ty) - BS(2,jj);\
      dx *= al2[0];\
      dy *= al2[1];\
      dz *= al2[2];
#elif MD_PERIODIC_FIXED==2
#define DISTANCE_PERIODIC_OLD \
      dx = AS(0,ty) - BS(0,jj);\
      dy = AS(1,ty) - BS(1,jj);\
      dz = AS(2,ty) - BS(2,jj);\
      dx = dx - rintf(dx * al) * l;\
      dy = dy - rintf(dy * al) * l;\
      dz = dz - rintf(dz * al) * l;
#elif MD_PERIODIC_FIXED==0
#define DISTANCE_PERIODIC_OLD \
      dx = AS(0,ty) - BS(0,jj);\
      dy = AS(1,ty) - BS(1,jj);\
      dz = AS(2,ty) - BS(2,jj);
#endif


#ifdef MD_USE_QAUNION
#define MD_KERNEL(gpu_funcname,gpu_kernel_core) \
extern "C" __global__ \
void gpu_funcname(int nii, int nj, int natj,\
		  VG_IVEC *ivec, VG_JVEC *jvec,\
		  VG_FVEC *fvec, VG_PSCALER *pscal)\
{\
  int by = blockIdx.x;\
  int ty = threadIdx.x;\
  int aBegin = VG_MINIMUM_PARTICLE_BLOCK_I*by+ty;\
  int jEnd = (nj+VG_MINIMUM_PARTICLE_BLOCK_J-1)/VG_MINIMUM_PARTICLE_BLOCK_J;\
  int jty;\
  int itype,k;\
  \
  float l = d_scalers->volume[0];\
  float al = 1.e0 / d_scalers->volume[0];\
  float dn2, dn6, alpha2=(d_scalers->alpha)*(d_scalers->alpha);	\
  float dx,dy,dz,tmp, eps2=d_scalers->eps2;\
  float api2=(float)(2.e0/sqrt(M_PI));\
  \
  float Csub[3],al2[3];\
  \
  DEFINE_AS_BS;\
  __shared__ int Ds[1][VG_MINIMUM_PARTICLE_BLOCK_I];\
  __shared__ int Es[1][VG_MINIMUM_PARTICLE_BLOCK_I];\
  __shared__ float Qis[1][VG_MINIMUM_PARTICLE_BLOCK_I];\
  __shared__ float Qjs[1][VG_MINIMUM_PARTICLE_BLOCK_I];\
  \
  al2[0]=scalbnf(d_scalers->volume[0],-32);\
  al2[1]=scalbnf(d_scalers->volume[1],-32);\
  al2[2]=scalbnf(d_scalers->volume[2],-32);\
  /*al2[0]*=(float)(1.0/((double)(1LL<<32)));*/\
  /*al2[1]*=(float)(1.0/((double)(1LL<<32)));*/\
  /*al2[2]*=(float)(1.0/((double)(1LL<<32)));*/\
  Csub[0] = Csub[1] = Csub[2] = 0.e0;\
  LOAD_AS;\
  /*DS(0,ty)=ivec[aBegin].qatype.atype;*/\
  /*QIS(0,ty)=ivec[aBegin].qatype.q;*/\
  LOAD_QATYPEI;\
  for (int j = 0; j < jEnd; j++){\
    jty=j*VG_MINIMUM_PARTICLE_BLOCK_I+ty;\
    LOAD_BS;\
    /*ES(0,ty)=jvec[jty].qatype.atype;*/\
    /*QJS(0,ty)=jvec[jty].qatype.q;*/\
    LOAD_QATYPEJ;\
    __syncthreads();\
    for (int jj = 0; jj < VG_MINIMUM_PARTICLE_BLOCK_J; jj++){\
      DISTANCE_PERIODIC_OLD;\
      itype = natj * DS(0,ty) + ES(0,jj);\
      gpu_kernel_core;\
    }\
    __syncthreads();\
  }\
  for(k=0;k<3;k++) fvec[aBegin].fi[k]=Csub[k];\
}

#if 0 // this is not used 
#define MD_KERNEL_CI_JLOOP_ONE(kernel) \
      DISTANCE_PERIODIC_OLD;\
      itype = natj * DS(0,ty) + ES(0,jj);\
      kernel;\
      jj++;
#if MD_UNROLL_KERNEL==8 // Please modify 
#define MD_KERNEL_CI_JLOOP(kernel) \
      MD_KERNEL_CI_JLOOP_ONE(kernel);\
      MD_KERNEL_CI_JLOOP_ONE(kernel);\
      MD_KERNEL_CI_JLOOP_ONE(kernel);\
      MD_KERNEL_CI_JLOOP_ONE(kernel);\
      MD_KERNEL_CI_JLOOP_ONE(kernel);\
      MD_KERNEL_CI_JLOOP_ONE(kernel);\
      MD_KERNEL_CI_JLOOP_ONE(kernel);\
      MD_KERNEL_CI_JLOOP_ONE(kernel);
#endif
#endif

#define MD_KERNEL_CI(gpu_funcname,gpu_kernel_core) \
extern "C" __global__ \
void gpu_funcname(int nii, int nj, int natj,\
		  VG_IVEC *ivec, VG_JVEC *jvec,\
		  VG_FVEC *fvec, VG_PSCALER *pscal)\
{\
  int by = blockIdx.x;\
  int ty = threadIdx.x;\
  int aBegin = VG_MINIMUM_PARTICLE_BLOCK_I*by+ty;\
  int jEnd = (nj+VG_MINIMUM_PARTICLE_BLOCK_J-1)/VG_MINIMUM_PARTICLE_BLOCK_J;\
  int jty;\
  int itype,k;\
  \
  float l = d_scalers->volume[0];\
  float al = 1.e0 / d_scalers->volume[0];\
  float dn2, dn6, alpha2=(d_scalers->alpha)*(d_scalers->alpha);	\
  float dx,dy,dz,tmp, eps2=d_scalers->eps2;\
  float api2=(float)(2.e0/sqrt(M_PI));\
  \
  float Csub[3],al2[3];\
  \
  DEFINE_AS_BS;\
  __shared__ int Ds[1][VG_MINIMUM_PARTICLE_BLOCK_I];\
  __shared__ int Es[1][VG_MINIMUM_PARTICLE_BLOCK_I];\
  __shared__ float Qis[1][VG_MINIMUM_PARTICLE_BLOCK_I];\
  __shared__ float Qjs[1][VG_MINIMUM_PARTICLE_BLOCK_I];\
  \
  al2[0]=scalbnf(d_scalers->volume[0],-32);\
  al2[1]=scalbnf(d_scalers->volume[1],-32);\
  al2[2]=scalbnf(d_scalers->volume[2],-32);\
  Csub[0] = Csub[1] = Csub[2] = 0.e0;\
  LOAD_AS;\
  /*DS(0,ty)=ivec[aBegin].qatype.atype;*/\
  /*QIS(0,ty)=ivec[aBegin].qatype.q;*/\
  LOAD_QATYPEI;\
  for(int jci=0;jci<27;jci++){\
   int joffset=pscal->ci[by][jci].base;\
   nj=pscal->ci[by][jci].size;\
   jEnd=(nj+VG_MINIMUM_PARTICLE_BLOCK_J-1)/VG_MINIMUM_PARTICLE_BLOCK_J;\
   for (int j = 0; j < jEnd; j++){\
    jty=j*VG_MINIMUM_PARTICLE_BLOCK_I+joffset+ty;\
    LOAD_BS;\
    /*ES(0,ty)=jvec[jty].qatype.atype;*/\
    /*QJS(0,ty)=jvec[jty].qatype.q;*/\
    LOAD_QATYPEJ;\
    __syncthreads();\
    for (int jj = 0; jj < VG_MINIMUM_PARTICLE_BLOCK_J;){\
      DISTANCE_PERIODIC_OLD;itype = natj * DS(0,ty) + ES(0,jj);gpu_kernel_core;jj++;\
      /* unrolling is not implemented yet */\
    }\
    __syncthreads();\
   }\
  }\
  for(k=0;k<3;k++) fvec[aBegin].fi[k]=Csub[k];\
}


#if defined(MD_MATRIX_IN_SCALER)
#define LJ_LOAD_MATRIX \
  if(ty<VG_MINIMUM_ATYPE_BLOCK){\
    AS(0,ty)=d_scalers->gscalesf[ty];\
    AS(1,ty)=d_scalers->gscalesp[ty];\
    AS(2,ty)=d_scalers->rscales[ty];\
  }\
  __syncthreads();
#elif defined(MD_MATRIX_COPY_TO_SHARED)
#define LJ_LOAD_MATRIX \
  if(ty<VG_MINIMUM_ATYPE_BLOCK){\
    AS(0,ty)=d_matrix[ty].gscale;\
    AS(2,ty)=d_matrix[ty].rscale;\
  }\
  __syncthreads();
#else
#define LJ_LOAD_MATRIX ;
#endif

#if MD_PERIODIC_FIXED==1
#define DEFINE_XI \
  int xi[3];\
  /*float al3=(float)(1.0 /((double)(1LL<<32)) * d_scalers->volume[0]);*/\
  float al3=scalbnf(d_scalers->volume[0],-32);\
  /*float al2[3];*/\
  al2[0]=scalbnf(d_scalers->volume[0],-32);\
  al2[1]=scalbnf(d_scalers->volume[1],-32);\
  al2[2]=scalbnf(d_scalers->volume[2],-32);\
  /*al2[0]*=(float)(1.0/((double)(1LL<<32)));*/\
  /*al2[1]*=(float)(1.0/((double)(1LL<<32)));*/\
  /*al2[2]*=(float)(1.0/((double)(1LL<<32)));*/
#if 0
  float al2[3]={scalbnf(d_scalers->volume[0],-32),\
                scalbnf(d_scalers->volume[1],-32),\
                scalbnf(d_scalers->volume[2],-32)};
#endif
#define LOAD_I \
  xi[0]=__float_as_int(jvecpoint[jty4  ]);\
  xi[1]=__float_as_int(jvecpoint[jty4+1]);\
  xi[2]=__float_as_int(jvecpoint[jty4+2]);
#elif MD_PERIODIC_FIXED==2
#define DEFINE_XI float xi[3];\
  al2[0]=1.0f/d_scalers->volume[0];\
  al2[1]=1.0f/d_scalers->volume[1];\
  al2[2]=1.0f/d_scalers->volume[2];
#define LOAD_I \
  xi[0]=jvecpoint[jty4  ];\
  xi[1]=jvecpoint[jty4+1];\
  xi[2]=jvecpoint[jty4+2];
#else
#define DEFINE_XI float xi[3];
#define LOAD_I \
  xi[0]=jvecpoint[jty4  ];\
  xi[1]=jvecpoint[jty4+1];\
  xi[2]=jvecpoint[jty4+2];
#endif

#ifdef MD_SIGMA_EPSION_IN_VGSTRUCT
#define LJ_LOAD_I \
  LOAD_I;\
  halfsigmai=jvecpoint[jty4+4];\
  sqrtepsiloni=jvecpoint[jty4+5];
#define LJ_LOAD_J \
    BS(0,ty4  )=jvecpoint[jty4  ];				\
    BS(0,ty4+1)=jvecpoint[jty4+1];				\
    BS(0,ty4+2)=jvecpoint[jty4+2];				\
    BS(0,ty4+3)=jvecpoint[jty4+3];				\
    BS(0,ty4+4)=jvecpoint[jty4+4];				\
    BS(0,ty4+5)=jvecpoint[jty4+5];				
#else // else of MD_SIGMA_EPSION_IN_VGSTRUCT
#ifdef MD_SORT_ATYPEI
#define LJ_LOAD_I \
  LOAD_I;\
  atypei=__float_as_int(jvecpoint[jty4+3]);\
  sameatypei=atypei & 0x100;\
  if(sameatypei!=0){\
    atypei=__float_as_int(jvecpoint[VG_MINIMUM_PARTICLE_BLOCK_I*by*LJ_BS_SIZE+3]) & 0xff;\
  }\
  atypei*=natj;
#define LJ_LOAD_J2 \
    BS(0,ty4  )=jvecpoint[jty4  ];				\
    BS(0,ty4+1)=jvecpoint[jty4+1];				\
    BS(0,ty4+2)=jvecpoint[jty4+2];				\
    BS(0,ty4+3)=jvecpoint[jty4+3];				\
    BS(0,ty4+3)=__int_as_float(__float_as_int(BS(0,ty4+3))+atypei);
#else // else of MD_SORT_ATYPEI
#if defined(MD_QAUNION_ATYPEBIT) && 1
#define LJ_LOAD_I \
  LOAD_I;\
  atypei=ATYPE_UNPACK_FROM_QATYPE(__float_as_int(jvecpoint[jty4+3]));\
  atypei*=natj;\
  Csub[3]=Q_INT_TO_FLOAT(__float_as_int(jvecpoint[jty4+3]),d_scalers->scaleqi_1);
#else // else of MD_QAUNION_ATYPEBIT
#define LJ_LOAD_I \
  LOAD_I;\
  atypei=__float_as_int(jvecpoint[jty4+3]);\
  atypei*=natj;
#endif // end of MD_QAUNION_ATYPEBIT
#endif // end of MD_SORT_ATYPEI
#ifdef MD_QAUNION_ATYPEBIT
#define LJ_LOAD_J \
    BS(0,ty4  )=jvecpoint[jty4  ];				\
    BS(0,ty4+1)=jvecpoint[jty4+1];				\
    BS(0,ty4+2)=jvecpoint[jty4+2];				\
    /*BS(0,ty4+3)=__int_as_float(ATYPE_UNPACK_FROM_QATYPE(__float_as_int(jvecpoint[jty4+3]))); */	\
    BS(0,ty4+3)=jvecpoint[jty4+3];
#define LJ_LOAD_J_ATYPEBIT \
    BS(0,ty4  )=jvecpoint[jty4  ];				\
    BS(0,ty4+1)=jvecpoint[jty4+1];				\
    BS(0,ty4+2)=jvecpoint[jty4+2];				\
    BS(0,ty4+3)=__int_as_float(ATYPE_UNPACK_FROM_QATYPE(__float_as_int(jvecpoint[jty4+3])));
#else
#define LJ_LOAD_J \
    BS(0,ty4  )=jvecpoint[jty4  ];				\
    BS(0,ty4+1)=jvecpoint[jty4+1];				\
    BS(0,ty4+2)=jvecpoint[jty4+2];				\
    BS(0,ty4+3)=jvecpoint[jty4+3];
#define LJ_LOAD_J_ATYPEBIT LJ_LOAD_J
#endif
#endif // end of MD_SIGMA_EPSION_IN_VGSTRUCT

#ifdef ACCUMULATE_PAIR_FLOAT
#if ACCUMULATE_PAIR_FLOAT==3
#define WRITE_BACK_TO_FVEC \
  for(k=0;k<MD_FVEC_SIZE;k++){((double *)(fvec[aBegin].fi))[k]=Csubd[k];}
#else
#define WRITE_BACK_TO_FVEC \
  for(k=0;k<MD_FVEC_SIZE;k++){fvec[aBegin].fi[k]=Csub[k];fvec[aBegin].fi2[k]=Csub2[k];}
#endif
#else
#define WRITE_BACK_TO_FVEC \
  for(k=0;k<MD_FVEC_SIZE;k++) fvec[aBegin].fi[k]=Csub[k];
#endif

#if MD_FVEC_SIZE==4
#ifdef ACCUMULATE_PAIR_FLOAT
#define INITIALIZE_CSUB \
  Csub[0] = Csub[1] = Csub[2] = Csub[3] = 0.0f;\
  Csub2[0] = Csub2[1] = Csub2[2] = Csub2[3] = 0.0f;
#else
#define INITIALIZE_CSUB \
  Csub[0] = Csub[1] = Csub[2] = Csub[3] = 0.0f;
#endif
#elif MD_FVEC_SIZE==3
#ifdef ACCUMULATE_PAIR_FLOAT
#if ACCUMULATE_PAIR_FLOAT==3
#define INITIALIZE_CSUB \
  Csubd[0] = Csubd[1] = Csubd[2] = 0.0;
#else // else of ACCUMULATE_PAIR_FLOAT==3
#if MD_USE_LARGE_VAL==1 || MD_USE_LARGE_VAL==2 || MD_USE_LARGE_VAL==20
#define INITIALIZE_CSUB \
  Csub[0] = Csub[1] = Csub[2] = LARGE_VAL;\
  Csub2[0] = Csub2[1] = Csub2[2] = 0.0f;
#else
#define INITIALIZE_CSUB \
  Csub[0] = Csub[1] = Csub[2] = 0.0f;\
  Csub2[0] = Csub2[1] = Csub2[2] = 0.0f;
#endif
#endif // end of ACCUMULATE_PAIR_FLOAT==3
#else // else of ACCUMULATE_PAIR_FLOAT
#define INITIALIZE_CSUB \
  Csub[0] = Csub[1] = Csub[2] = 0.0f;
#endif // end of ACCUMULATE_PAIR_FLOAT
#endif

#if MD_USE_LARGE_VAL==20
#define INITIALIZE_CSUB_COULOMB \
  Csub[0] = Csub[1] = Csub[2] = 0.0f;\
  Csub2[0] = Csub2[1] = Csub2[2] = 0.0f;
#define INITIALIZE_CSUB_VDW     INITIALIZE_CSUB
#else
#define INITIALIZE_CSUB_COULOMB INITIALIZE_CSUB
#define INITIALIZE_CSUB_VDW     INITIALIZE_CSUB
#endif

#ifdef ACCUMULATE_PAIR_FLOAT
#if ACCUMULATE_PAIR_FLOAT==3
#define DEFINE_CSUB \
  double Csubd[MD_FVEC_SIZE];\
  float *Csub=(float *)Csubd,*Csub2=((float *)Csubd)+MD_FVEC_SIZE;
#else // else of ACCUMULATE_PAIR_FLOAT==3
#define DEFINE_CSUB \
  float Csub[MD_FVEC_SIZE+1],Csub2[MD_FVEC_SIZE];\
  double *Csubd=NULL;
#endif // end of ACCUMULATE_PAIR_FLOAT==3
#else // else of ACCUMULATE_PAIR_FLOAT
#define DEFINE_CSUB \
  float Csub[MD_FVEC_SIZE+1],Csub2[MD_FVEC_SIZE];\
  double *Csubd=NULL;
#endif // end of ACCUMULATE_PAIR_FLOAT

/* work
#define UPDATE_FORCE_CORE \
          for(int k=0;k<3;k++){;\
            Csub[k]+=(float)(__float_as_int(Csub2[k]) & (MASK(LOWER_VAL_SHIFT)\
              <<(32-LOWER_VAL_SHIFT)))*((float)LOWER_VAL_FACTOR_1);\
            Csub2[k]=__int_as_float(__float_as_int(Csub2[k]) & MASK(32-LOWER_VAL_SHIFT));\
	  }
*/
#define UPDATE_FORCE_CORE \
          for(int k=0;k<3;k++){;\
            Csub[k]+=(float)(__float_as_int(Csub2[k]) & (MASK(LOWER_VAL_SHIFT)\
              <<(32-LOWER_VAL_SHIFT)))*(float)LOWER_VAL_FACTOR_1;\
            Csub2[k]=__int_as_float(__float_as_int(Csub2[k]) & MASK(32-LOWER_VAL_SHIFT));\
	  }
#if MD_USE_LARGE_VAL==2 
#define UPDATE_FORCE_COULOMB UPDATE_FORCE_CORE
#define UPDATE_FORCE_VDW     UPDATE_FORCE_CORE
#elif MD_USE_LARGE_VAL==20
#define UPDATE_FORCE_COULOMB
#define UPDATE_FORCE_VDW     UPDATE_FORCE_CORE
#else
#define UPDATE_FORCE_COULOMB
#define UPDATE_FORCE_VDW
#endif

#define UNROLL_KERNEL4(kernel) \
  kernel;\
  kernel;\
  kernel;\
  kernel;
#define UNROLL_KERNEL8(kernel) \
  UNROLL_KERNEL4(kernel);\
  UNROLL_KERNEL4(kernel);
#define UNROLL_KERNEL16(kernel) \
  UNROLL_KERNEL8(kernel);\
  UNROLL_KERNEL8(kernel);
#define UNROLL_KERNEL32(kernel) \
  UNROLL_KERNEL16(kernel);\
  UNROLL_KERNEL16(kernel);
#define UNROLL_KERNEL64(kernel) \
  UNROLL_KERNEL132(kernel);\
  UNROLL_KERNEL132(kernel);
#define UNROLL_KERNEL128(kernel) \
  UNROLL_KERNEL164(kernel);\
  UNROLL_KERNEL164(kernel);
#define JJMAX VG_MINIMUM_PARTICLE_BLOCK_I/MD_UNROLL_KERNEL

#ifdef MD_UNROLL_KERNEL
#if MD_UNROLL_KERNEL==4
#define OUTER_UNROLL_KERNEL_CORE(kernel) for(jj=0;jj<JJMAX;jj++){UNROLL_KERNEL4(kernel);}
#elif MD_UNROLL_KERNEL==8
#define OUTER_UNROLL_KERNEL_CORE(kernel) for(jj=0;jj<JJMAX;jj++){UNROLL_KERNEL8(kernel);}
#elif MD_UNROLL_KERNEL==16
#define OUTER_UNROLL_KERNEL_CORE(kernel) for(jj=0;jj<JJMAX;jj++){UNROLL_KERNEL16(kernel);}
#elif MD_UNROLL_KERNEL==32
#define OUTER_UNROLL_KERNEL_CORE(kernel) for(jj=0;jj<JJMAX;jj++){UNROLL_KERNEL32(kernel);}
#elif MD_UNROLL_KERNEL==64
#define OUTER_UNROLL_KERNEL_CORE(kernel) for(jj=0;jj<JJMAX;jj++){UNROLL_KERNEL64(kernel);}
#elif MD_UNROLL_KERNEL==128
#define OUTER_UNROLL_KERNEL_CORE(kernel) for(jj=0;jj<JJMAX;jj++){UNROLL_KERNEL128(kernel);}
#endif
#else // else of MD_UNROLL_KERNEL
#define OUTER_UNROLL_KERNEL_CORE(kernel) for(jj=0;jj<VG_MINIMUM_PARTICLE_BLOCK_I;jj++){kernel;}
#endif // end of MD_UNROLL_KERNEL

#define OUTER_UNROLL_KERNEL_COULOMB(kernel) \
          OUTER_UNROLL_KERNEL_CORE(kernel);\
          UPDATE_FORCE_COULOMB;
#define OUTER_UNROLL_KERNEL_VDW(kernel) \
          OUTER_UNROLL_KERNEL_CORE(kernel);\
          UPDATE_FORCE_VDW;

#if MD_PERIODIC_FIXED==1
#define XITYPE int
#else
#define XITYPE float
#endif

#define MD_KERNEL_FB_SUB(gpu_funcname_sub,gpu_kernel_core,lj_load_j,flag) \
__device__ \
void gpu_funcname_sub(int ty4, int nj, int jEnd, int atypei,\
		      float *jvecpoint,\
		      float Bs[1][VG_MINIMUM_PARTICLE_BLOCK_I*LJ_BS_SIZE],\
		      float As[3][VG_MINIMUM_ATYPE_BLOCK], float al2[3],\
                      float size[3], float r2max, float rcut21,          \
		      XITYPE xi[3], float Csub[], float Csub2[], double Csubd[])\
{\
  int j,jj,jty4,jj4,itype;\
  float cutoff,dn2,dn6,tmp,rij,dx,dy,dz,tmp2,tmp3,tmp4,tmp5;	\
  if(flag){\
   for (j = 0; j < jEnd; j++){		\
    jty4=j*VG_MINIMUM_PARTICLE_BLOCK_I*LJ_BS_SIZE+ty4;	\
    lj_load_j;\
    __syncthreads();\
    jj4=0;\
    OUTER_UNROLL_KERNEL_VDW(gpu_kernel_core);\
    __syncthreads();\
   }\
  }\
  else{\
   for (j = 0; j < jEnd; j++){		\
    jty4=j*VG_MINIMUM_PARTICLE_BLOCK_I*LJ_BS_SIZE+ty4;	\
    lj_load_j;\
    __syncthreads();\
    jj4=0;\
    OUTER_UNROLL_KERNEL_VDW(gpu_kernel_core);\
    __syncthreads();\
   }\
    jty4=j*VG_MINIMUM_PARTICLE_BLOCK_I*LJ_BS_SIZE+ty4;	\
    lj_load_j;\
    __syncthreads();\
    jj4=0;\
    for(jj=0;jj<nj-j*VG_MINIMUM_PARTICLE_BLOCK_J;jj++){\
      gpu_kernel_core;\
    }\
    __syncthreads();\
  }\
}


#define MD_KERNEL_FB_SUB_CI(gpu_funcname_sub,gpu_kernel_core,lj_load_j) \
__device__ \
void gpu_funcname_sub(int ty4, int nj, int jEnd,\
                      VG_PSCALER *pscal, int atypei,\
		      float *jvecpoint,\
		      float Bs[1][VG_MINIMUM_PARTICLE_BLOCK_I*LJ_BS_SIZE],\
		      float As[3][VG_MINIMUM_ATYPE_BLOCK], float al2[3],\
                      float size[3], float r2max, float rcut21,          \
		      XITYPE xi[3], float Csub[], float Csub2[], double Csubd[])\
{\
  int j,jj,jty4,jj4,itype;\
  float cutoff,dn2,dn6,tmp,rij,dx,dy,dz,tmp2,tmp3,tmp4,tmp5;	\
  int by = blockIdx.x;\
  for(int jci=0;jci<27;jci++){\
   int joffset=pscal->ci[by][jci].base;\
   nj=pscal->ci[by][jci].size;\
   jEnd=(nj+VG_MINIMUM_PARTICLE_BLOCK_J-1)/VG_MINIMUM_PARTICLE_BLOCK_J;\
   for (j = 0; j < jEnd; j++){		\
    jty4=(j*VG_MINIMUM_PARTICLE_BLOCK_I+joffset)*LJ_BS_SIZE+ty4;	\
    lj_load_j;\
    __syncthreads();\
    jj4=0;\
    OUTER_UNROLL_KERNEL_VDW(gpu_kernel_core);\
    __syncthreads();\
   }\
  }\
}


#define MD_KERNEL_CVFB_SUB_CI(gpu_funcname_sub,gpu_kernel_core,lj_load_j) \
__device__ \
void gpu_funcname_sub(int ty4, int nj, int jEnd,\
                      VG_PSCALER *pscal, int atypei,\
		      float *jvecpoint,\
		      float Bs[1][VG_MINIMUM_PARTICLE_BLOCK_I*LJ_BS_SIZE],\
		      float As[3][VG_MINIMUM_ATYPE_BLOCK], float al2[3],\
                      float size[3], float r2max, float rcut21,          \
                      float alpha2, float api2, float scaleqj_1, float alpha3fac,\
		      XITYPE xi[3], float Csub[], float Csub2[], double Csubd[])\
{\
  int j,jj,jty4,jj4,itype;\
  float cutoff,dn2,dn6,tmp,rij,dx,dy,dz,tmp2,tmp3,tmp4,tmp5;	\
  float sqdn,qj;\
  int by = blockIdx.x;\
  for(int jci=0;jci<27;jci++){\
   int joffset=pscal->ci[by][jci].base;\
   nj=pscal->ci[by][jci].size;\
   jEnd=(nj+VG_MINIMUM_PARTICLE_BLOCK_J-1)/VG_MINIMUM_PARTICLE_BLOCK_J;\
   for (j = 0; j < jEnd; j++){		\
    jty4=(j*VG_MINIMUM_PARTICLE_BLOCK_I+joffset)*LJ_BS_SIZE+ty4;	\
    lj_load_j;\
    __syncthreads();\
    jj4=0;\
    OUTER_UNROLL_KERNEL_VDW(gpu_kernel_core);\
    __syncthreads();\
   }\
  }\
}


#define MD_KERNEL_FB(gpu_funcname,gpu_kernel_core,flag) \
extern "C" __global__ \
void gpu_funcname(int nii, int nj, int natj,\
		  VG_IVEC *ivec, VG_JVEC *jvec,\
		  VG_FVEC *fvec, VG_PSCALER *pscal)\
{\
  int by = blockIdx.x;\
  int ty = threadIdx.x;\
  int aBegin = VG_MINIMUM_PARTICLE_BLOCK_I*by+ty;;\
  /*int jEnd = (nj+VG_MINIMUM_PARTICLE_BLOCK_J-1)/VG_MINIMUM_PARTICLE_BLOCK_J;*/\
  int jEnd = nj/VG_MINIMUM_PARTICLE_BLOCK_J;\
  int j,jj,jty,jty4,jty8,ty4,ty8,jj4;					\
  int k,atypei,itype;							\
  float dx, dy, dz, eps=d_scalers->eps2;	\
  float api2=(float)(2.e0/sqrt(M_PI)),kkpi1=1.0/(2.0*M_PI*sqrt(2.0*M_PI)); \
  float rij,sij,pi4=4.0/M_PI,pi14=0.25/M_PI;\
  float l = d_scalers->volume[0];\
  float al = 1.e0 / d_scalers->volume[0];\
  float size[3]={d_scalers->volume[0],d_scalers->volume[1],d_scalers->volume[2]};\
  float dn6,al2[3],tmp2,tmp3,tmp4;\
  float alpha2=(d_scalers->alpha)*(d_scalers->alpha);\
  float r2max=d_scalers->r2max;\
  float rcut21=d_scalers->rcut21;\
  DEFINE_CSUB;\
  /*double Csubd[MD_FVEC_SIZE];*/\
  /*float *Csub=(float *)Csubd,*Csub2=((float *)Csubd)+MD_FVEC_SIZE;*/\
  /*float Csub[MD_FVEC_SIZE],Csub2[MD_FVEC_SIZE];*/\
  /*double *Csubd=NULL;*/\
  float cutoff,kk1=0.5f,pi=(float)M_PI; \
  float tmp,*jvecpoint,halfsigmai,halfsigmaj,sqrtepsiloni,sqrtepsilonj;\
  __shared__ float Bs[1][VG_MINIMUM_PARTICLE_BLOCK_I*LJ_BS_SIZE];\
  __shared__ float As[3][VG_MINIMUM_ATYPE_BLOCK];\
  DEFINE_XI;\
  ty4=ty*LJ_BS_SIZE;\
  INITIALIZE_CSUB_COULOMB;\
  jvecpoint=(float *)ivec;\
  jty4=aBegin*LJ_BS_SIZE;\
  LOAD_I;\
  jvecpoint=(float *)jvec;\
  if(flag){\
   for (j = 0; j < jEnd; j++){		\
    jty4=j*VG_MINIMUM_PARTICLE_BLOCK_I*LJ_BS_SIZE+ty4;	\
    LJ_LOAD_J;\
    __syncthreads();\
    jj4=0;\
    OUTER_UNROLL_KERNEL_COULOMB(gpu_kernel_core);\
    /*gpu_kernel_core;*/\
    __syncthreads();\
   }\
  }\
  else{\
   for (j = 0; j < jEnd; j++){		\
    jty4=j*VG_MINIMUM_PARTICLE_BLOCK_I*LJ_BS_SIZE+ty4;	\
    LJ_LOAD_J;\
    __syncthreads();\
    jj4=0;\
    OUTER_UNROLL_KERNEL_COULOMB(gpu_kernel_core);\
    /*gpu_kernel_core;*/\
    __syncthreads();\
   }\
    jty4=j*VG_MINIMUM_PARTICLE_BLOCK_I*LJ_BS_SIZE+ty4;	\
    LJ_LOAD_J;\
    __syncthreads();\
    jj4=0;\
    for(jj=0;jj<nj-j*VG_MINIMUM_PARTICLE_BLOCK_J;jj++){\
      gpu_kernel_core;\
    }\
    /*gpu_kernel_core;*/\
    __syncthreads();\
  }\
  WRITE_BACK_TO_FVEC;\
  /*for(k=0;k<3;k++) fvec[aBegin].fi[k]=Csub[k];*/\
}


#define MD_KERNEL_FB_CI(gpu_funcname,gpu_kernel_core) \
extern "C" __global__ \
void gpu_funcname(int nii, int nj, int natj,\
		  VG_IVEC *ivec, VG_JVEC *jvec,\
		  VG_FVEC *fvec, VG_PSCALER *pscal)\
{\
  int by = blockIdx.x;\
  int ty = threadIdx.x;\
  int aBegin = VG_MINIMUM_PARTICLE_BLOCK_I*by+ty;;\
  /*int jEnd = (nj+VG_MINIMUM_PARTICLE_BLOCK_J-1)/VG_MINIMUM_PARTICLE_BLOCK_J;*/\
  int jEnd = nj/VG_MINIMUM_PARTICLE_BLOCK_J;\
  int j,jj,jty,jty4,jty8,ty4,ty8,jj4;					\
  int k,atypei,itype;							\
  float dx, dy, dz, eps=d_scalers->eps2;	\
  float api2=(float)(2.e0/sqrt(M_PI)),kkpi1=1.0/(2.0*M_PI*sqrt(2.0*M_PI)); \
  float rij,sij,pi4=4.0/M_PI,pi14=0.25/M_PI;\
  float l = d_scalers->volume[0];\
  float al = 1.e0 / d_scalers->volume[0];\
  float size[3]={d_scalers->volume[0],d_scalers->volume[1],d_scalers->volume[2]};\
  float dn6,al2[3],tmp2,tmp3,tmp4;\
  float alpha2=(d_scalers->alpha)*(d_scalers->alpha);\
  float r2max=d_scalers->r2max;\
  float rcut21=d_scalers->rcut21;\
  DEFINE_CSUB;\
  /*double Csubd[MD_FVEC_SIZE];*/\
  /*float *Csub=(float *)Csubd,*Csub2=((float *)Csubd)+MD_FVEC_SIZE;*/\
  /*float Csub[MD_FVEC_SIZE],Csub2[MD_FVEC_SIZE];*/\
  /*double *Csubd=NULL;*/\
  float cutoff,kk1=0.5f,pi=(float)M_PI; \
  float tmp,*jvecpoint,halfsigmai,halfsigmaj,sqrtepsiloni,sqrtepsilonj;\
  __shared__ float Bs[1][VG_MINIMUM_PARTICLE_BLOCK_I*LJ_BS_SIZE];\
  __shared__ float As[3][VG_MINIMUM_ATYPE_BLOCK];\
  DEFINE_XI;\
  ty4=ty*LJ_BS_SIZE;\
  INITIALIZE_CSUB_COULOMB;\
  jvecpoint=(float *)ivec;\
  jty4=aBegin*LJ_BS_SIZE;\
  LOAD_I;\
  jvecpoint=(float *)jvec;\
  for(int jci=0;jci<27;jci++){\
   int joffset=pscal->ci[by][jci].base;\
   nj=pscal->ci[by][jci].size;\
   jEnd=(nj+VG_MINIMUM_PARTICLE_BLOCK_J-1)/VG_MINIMUM_PARTICLE_BLOCK_J;\
   for (j = 0; j < jEnd; j++){		\
    jty4=(j*VG_MINIMUM_PARTICLE_BLOCK_I+joffset)*LJ_BS_SIZE+ty4;	\
    LJ_LOAD_J;\
    __syncthreads();\
    jj4=0;\
    OUTER_UNROLL_KERNEL_COULOMB(gpu_kernel_core);\
    __syncthreads();\
   }\
  }\
  WRITE_BACK_TO_FVEC;\
}


#define MD_KERNEL_FBLJ(gpu_funcname,gpu_kernel_core,gpu_kernel_core2) \
extern "C" __global__ \
void gpu_funcname(int nii, int nj, int natj,\
		  VG_IVEC *ivec, VG_JVEC *jvec,\
		  VG_FVEC *fvec, VG_PSCALER *pscal)\
{\
  int by = blockIdx.x;\
  int ty = threadIdx.x;\
  int aBegin = VG_MINIMUM_PARTICLE_BLOCK_I*by+ty;;\
  /*int jEnd = (nj+VG_MINIMUM_PARTICLE_BLOCK_J-1)/VG_MINIMUM_PARTICLE_BLOCK_J;*/\
  int jEnd = nj/VG_MINIMUM_PARTICLE_BLOCK_J;\
  int jty,jty4,jty8,ty4,ty8,jj4,jj;					\
  int k,atypei,itype;							\
  float dx, dy, dz, eps=d_scalers->eps2;	\
  float api2=(float)(2.e0/sqrt(M_PI)),kkpi1=1.0/(2.0*M_PI*sqrt(2.0*M_PI)); \
  float rij,sij,pi4=4.0/M_PI,pi14=0.25/M_PI;\
  float l = d_scalers->volume[0];\
  float al = 1.e0 / d_scalers->volume[0];\
  float size[3]={d_scalers->volume[0],d_scalers->volume[1],d_scalers->volume[2]};\
  float r2max=d_scalers->r2max;\
  float rcut21=d_scalers->rcut21;\
  float dn6;\
  DEFINE_CSUB;\
  /*double Csubd[MD_FVEC_SIZE];*/\
  /*float *Csub=(float *)Csubd,*Csub2=((float *)Csubd)+MD_FVEC_SIZE;*/\
  /*float Csub[MD_FVEC_SIZE],Csub2[MD_FVEC_SIZE];*/\
  /*double *Csubd=NULL;*/\
  float cutoff,kk1=0.5f,pi=(float)M_PI; \
  float tmp,*jvecpoint,halfsigmai,halfsigmaj,sqrtepsiloni,sqrtepsilonj;\
  __shared__ float Bs[1][VG_MINIMUM_PARTICLE_BLOCK_I*LJ_BS_SIZE];\
  __shared__ float As[3][VG_MINIMUM_ATYPE_BLOCK];\
  float tmp2,tmp3,tmp4,al2[3];\
  int sameatypei=0;\
  DEFINE_XI;\
  \
  \
  ty4=ty*LJ_BS_SIZE;\
  \
  INITIALIZE_CSUB_VDW;\
  \
  jvecpoint=(float *)ivec;\
  jty4=aBegin*LJ_BS_SIZE;\
  LJ_LOAD_I;\
  LJ_LOAD_MATRIX;\
  \
  jvecpoint=(float *)jvec;\
  if(sameatypei==0){\
  /*if(by<10){*/\
  /*if(1){*/\
    gpu_kernel_core(ty4,nj,jEnd,atypei,jvecpoint,Bs,As,al2,size,r2max,rcut21,xi,Csub,Csub2,Csubd); \
  }\
  else{\
    gpu_kernel_core2(ty4,nj,jEnd,atypei,jvecpoint,Bs,As,al2,size,r2max,rcut21,xi,Csub,Csub2,Csubd); \
  }\
  WRITE_BACK_TO_FVEC;\
}

#define MD_KERNEL_FBLJ_CI(gpu_funcname,gpu_kernel_core,gpu_kernel_core2) \
extern "C" __global__ \
void gpu_funcname(int nii, int nj, int natj,\
		  VG_IVEC *ivec, VG_JVEC *jvec,\
		  VG_FVEC *fvec, VG_PSCALER *pscal)\
{\
  int by = blockIdx.x;\
  int ty = threadIdx.x;\
  int aBegin = VG_MINIMUM_PARTICLE_BLOCK_I*by+ty;;\
  /*int jEnd = (nj+VG_MINIMUM_PARTICLE_BLOCK_J-1)/VG_MINIMUM_PARTICLE_BLOCK_J;*/\
  int jEnd = nj/VG_MINIMUM_PARTICLE_BLOCK_J;\
  int jty,jty4,jty8,ty4,ty8,jj4,jj;					\
  int k,atypei,itype;							\
  float dx, dy, dz, eps=d_scalers->eps2;	\
  float api2=(float)(2.e0/sqrt(M_PI)),kkpi1=1.0/(2.0*M_PI*sqrt(2.0*M_PI)); \
  float rij,sij,pi4=4.0/M_PI,pi14=0.25/M_PI;\
  float l = d_scalers->volume[0];\
  float al = 1.e0 / d_scalers->volume[0];\
  float size[3]={d_scalers->volume[0],d_scalers->volume[1],d_scalers->volume[2]};\
  float r2max=d_scalers->r2max;\
  float rcut21=d_scalers->rcut21;\
  float dn6;\
  DEFINE_CSUB;\
  /*double Csubd[MD_FVEC_SIZE];*/\
  /*float *Csub=(float *)Csubd,*Csub2=((float *)Csubd)+MD_FVEC_SIZE;*/\
  /*float Csub[MD_FVEC_SIZE],Csub2[MD_FVEC_SIZE];*/\
  /*double *Csubd=NULL;*/\
  float cutoff,kk1=0.5f,pi=(float)M_PI; \
  float tmp,*jvecpoint,halfsigmai,halfsigmaj,sqrtepsiloni,sqrtepsilonj;\
  __shared__ float Bs[1][VG_MINIMUM_PARTICLE_BLOCK_I*LJ_BS_SIZE];\
  __shared__ float As[3][VG_MINIMUM_ATYPE_BLOCK];\
  float tmp2,tmp3,tmp4,al2[3];\
  int sameatypei=0;\
  DEFINE_XI;\
  \
  \
  ty4=ty*LJ_BS_SIZE;\
  \
  INITIALIZE_CSUB_VDW;\
  \
  jvecpoint=(float *)ivec;\
  jty4=aBegin*LJ_BS_SIZE;\
  LJ_LOAD_I;\
  LJ_LOAD_MATRIX;\
  \
  jvecpoint=(float *)jvec;\
  if(sameatypei==0){\
  /*if(by<10){*/\
  /*if(1){*/\
    gpu_kernel_core(ty4,nj,jEnd,pscal,atypei,jvecpoint,Bs,As,al2,size,r2max,rcut21,xi,Csub,Csub2,Csubd); \
  }\
  else{\
    gpu_kernel_core2(ty4,nj,jEnd,pscal,atypei,jvecpoint,Bs,As,al2,size,r2max,rcut21,xi,Csub,Csub2,Csubd); \
  }\
  WRITE_BACK_TO_FVEC;\
}

#define MD_KERNEL_CVFB_CI(gpu_funcname,gpu_kernel_core,gpu_kernel_core2) \
extern "C" __global__ \
void gpu_funcname(int nii, int nj, int natj,\
		  VG_IVEC *ivec, VG_JVEC *jvec,\
		  VG_FVEC *fvec, VG_PSCALER *pscal)\
{\
  int by = blockIdx.x;\
  int ty = threadIdx.x;\
  int aBegin = VG_MINIMUM_PARTICLE_BLOCK_I*by+ty;;\
  /*int jEnd = (nj+VG_MINIMUM_PARTICLE_BLOCK_J-1)/VG_MINIMUM_PARTICLE_BLOCK_J;*/\
  int jEnd = nj/VG_MINIMUM_PARTICLE_BLOCK_J;\
  int jty,jty4,jty8,ty4,ty8,jj4,jj;					\
  int k,atypei,itype;							\
  float dx, dy, dz, eps=d_scalers->eps2;	\
  float api2=(float)(2.e0/sqrt(M_PI)),kkpi1=1.0/(2.0*M_PI*sqrt(2.0*M_PI)); \
  float rij,sij,pi4=4.0/M_PI,pi14=0.25/M_PI;\
  float l = d_scalers->volume[0];\
  float al = 1.e0 / d_scalers->volume[0];\
  float size[3]={d_scalers->volume[0],d_scalers->volume[1],d_scalers->volume[2]};\
  float r2max=d_scalers->r2max;\
  float rcut21=d_scalers->rcut21;\
  float alpha2=(d_scalers->alpha)*(d_scalers->alpha);\
  float scaleqj_1=d_scalers->scaleqj_1;\
  float alpha3fac=d_scalers->alpha3fac;\
  float dn6;\
  DEFINE_CSUB;\
  /*double Csubd[MD_FVEC_SIZE];*/\
  /*float *Csub=(float *)Csubd,*Csub2=((float *)Csubd)+MD_FVEC_SIZE;*/\
  /*float Csub[MD_FVEC_SIZE],Csub2[MD_FVEC_SIZE];*/\
  /*double *Csubd=NULL;*/\
  float cutoff,kk1=0.5f,pi=(float)M_PI; \
  float tmp,*jvecpoint,halfsigmai,halfsigmaj,sqrtepsiloni,sqrtepsilonj;\
  __shared__ float Bs[1][VG_MINIMUM_PARTICLE_BLOCK_I*LJ_BS_SIZE];\
  __shared__ float As[3][VG_MINIMUM_ATYPE_BLOCK];\
  float tmp2,tmp3,tmp4,al2[3];\
  int sameatypei=0;\
  DEFINE_XI;\
  \
  \
  ty4=ty*LJ_BS_SIZE;\
  \
  INITIALIZE_CSUB_VDW;\
  \
  jvecpoint=(float *)ivec;\
  jty4=aBegin*LJ_BS_SIZE;\
  LJ_LOAD_I;\
  LJ_LOAD_MATRIX;\
  \
  jvecpoint=(float *)jvec;\
  if(sameatypei==0){\
  /*if(by<10){*/\
  /*if(1){*/\
    gpu_kernel_core(ty4,nj,jEnd,pscal,atypei,jvecpoint,Bs,As,al2,size,r2max,rcut21,alpha2,api2,scaleqj_1,alpha3fac,xi,Csub,Csub2,Csubd); \
  }\
  else{\
    gpu_kernel_core2(ty4,nj,jEnd,pscal,atypei,jvecpoint,Bs,As,al2,size,r2max,rcut21,alpha2,api2,scaleqj_1,alpha3fac,xi,Csub,Csub2,Csubd); \
  }\
  WRITE_BACK_TO_FVEC;\
}


#if VG_MINIMUM_PARTICLE_BLOCK_J==64 && VG_JDIV==4 && MD_UNROLL_KERNEL==8
#define MD_KERNEL_CVFB_SUB_CI_JDIV(gpu_funcname_sub,gpu_kernel_core,lj_load_j) \
__device__ \
void gpu_funcname_sub(int ty4, int nj, int natj,	\
                      VG_PSCALER *pscal, int atypei,\
		      float *Xis, float *jvecpoint,			\
		      float Bs[1][VG_MINIMUM_PARTICLE_BLOCK_I*LJ_BS_SIZE], \
		      float As[3][VG_MINIMUM_ATYPE_BLOCK], float al2[3], \
		      float size[3], float r2max, float rcut21,          \
		      float alpha2, float api2, float scaleqj_1, float alpha3fac, \
		      XITYPE xi[3], float Csub[], float Csub2[], double Csubd[]) \
{\
  int j,jj,jty4,jj4,itype,jEnd,ii,ii4,iii,jjj;				\
  float cutoff,dn2,dn6,tmp,rij,dx,dy,dz,tmp2,tmp3,tmp4,tmp5;	\
  float sqdn,qj;\
  int by = blockIdx.x;\
  int ty = threadIdx.x, jty;		\
  int ni=pscal->niblock[by],nimax;			\
  int aBegin = VG_MINIMUM_PARTICLE_BLOCK_I*by+ty;\
  float fi[VG_JDIV][3],fi2[VG_JDIV][3];\
  __shared__ VG_CELL cis[MD_NUM_JCELLS];\
  nimax=(ni+VG_MINIMUM_PARTICLE_BLOCK_J/VG_JDIV-1)/(VG_MINIMUM_PARTICLE_BLOCK_J/VG_JDIV);  \
  jty4=aBegin*LJ_BS_SIZE;\
  __syncthreads();\
  if(ty<MD_NUM_JCELLS*2){\
    if((ty % 2)==0) cis[ty/2].base=pscal->ci[by][ty/2].base;\
    else            cis[ty/2].size=pscal->ci[by][ty/2].size;\
  } __syncthreads();\
  for(int idiv=0;idiv<nimax;idiv++){\
    iii=(ty % (VG_MINIMUM_PARTICLE_BLOCK_J/VG_JDIV))+idiv*VG_MINIMUM_PARTICLE_BLOCK_J/VG_JDIV;\
    ii4=iii*LJ_BS_SIZE;\
    xi[0]=Xis[ii4+0];\
    xi[1]=Xis[ii4+1];\
    xi[2]=Xis[ii4+2];\
    atypei=ATYPE_UNPACK_FROM_QATYPE(__float_as_int(Xis[ii4+3]));\
    atypei*=natj;\
    for(int k=0;k<3;k++){Csub[k]=LARGE_VAL;Csub2[k]=0.0f;}\
    Csub[3]=Q_INT_TO_FLOAT(__float_as_int(Xis[ii4+3]),d_scalers->scaleqi_1);\
    /*for(int jci=0;jci<27;jci++){*/					\
    for(int jci=0;jci<cis[27].size;jci++){				\
    /*for(int jci=0;jci<1;jci++){*/					\
      int joffset=cis[jci].base;\
      nj=cis[jci].size;\
      jEnd=(nj+VG_MINIMUM_PARTICLE_BLOCK_J/VG_JDIV-1)/(VG_MINIMUM_PARTICLE_BLOCK_J/VG_JDIV);  \
      for (j = 0; j < jEnd; j++){		\
        jty=(j*VG_MINIMUM_PARTICLE_BLOCK_I/VG_JDIV+joffset)*LJ_BS_SIZE+ty;  \
        BS(0,ty)=jvecpoint[jty];						\
        __syncthreads();\
/*        jj4=0;for(jj=0;jj<VG_MINIMUM_PARTICLE_BLOCK_I/MD_UNROLL_KERNEL/VG_JDIV;jj++){UNROLL_KERNEL8(gpu_kernel_core);UPDATE_FORCE_VDW;} __syncthreads();*/\
\
        jj4=0;\
/*        UNROLL_KERNEL16(gpu_kernel_core);*/\
        UNROLL_KERNEL4(gpu_kernel_core);\
        UPDATE_FORCE_VDW;\
        __syncthreads();\
      } \
    }\
/*    if(ty/(VG_MINIMUM_PARTICLE_BLOCK_I/VG_JDIV)!=idiv) for(int k=0;k<3;k++){Csub[k]=LARGE_VAL;Csub2[k]=0.0f;}*/\
    for(int k=0;k<3;k++){ fi[idiv][k]=Csub[k];fi2[idiv][k]=Csub2[k];}\
  }\
  /* add 4 partial sums */\
  for(int idiv=0;idiv<nimax;idiv++){\
    __syncthreads();\
    for(int k=0;k<3;k++) Xis[ty4+k]=fi[idiv][k];\
    __syncthreads();\
    if(ty/(VG_MINIMUM_PARTICLE_BLOCK_I/VG_JDIV)==idiv){\
      for(int k=0;k<3;k++) Csub[k]=LARGE_VAL;\
      for(iii=0;iii<VG_JDIV;iii++){\
        for(int k=0;k<3;k++) Csub[k]+=Xis[(iii*VG_MINIMUM_PARTICLE_BLOCK_J/VG_JDIV+(ty % (VG_MINIMUM_PARTICLE_BLOCK_J/VG_JDIV)))*4+k]-LARGE_VAL;\
      }      \
/*      for(int k=0;k<3;k++) Csub[k]=Xis[(idiv*VG_MINIMUM_PARTICLE_BLOCK_J/VG_JDIV+(ty % (VG_MINIMUM_PARTICLE_BLOCK_J/VG_JDIV)))*4+k];*/\
    }\
    __syncthreads();\
    for(int k=0;k<3;k++) Xis[ty4+k]=fi2[idiv][k];\
    __syncthreads();\
    if(ty/(VG_MINIMUM_PARTICLE_BLOCK_I/VG_JDIV)==idiv){\
      for(int k=0;k<3;k++) Csub2[k]=__int_as_float(0);\
      for(iii=0;iii<VG_JDIV;iii++){\
        for(int k=0;k<3;k++) Csub2[k]=__int_as_float(__float_as_int(Csub2[k])+__float_as_int(Xis[(iii*VG_MINIMUM_PARTICLE_BLOCK_J/VG_JDIV+(ty % (VG_MINIMUM_PARTICLE_BLOCK_J/VG_JDIV)))*4+k]));\
      }     \
      /*for(int k=0;k<3;k++) Csub2[k]=Xis[(idiv*VG_MINIMUM_PARTICLE_BLOCK_J/VG_JDIV+(ty % (VG_MINIMUM_PARTICLE_BLOCK_J/VG_JDIV)))*4+k];*/\
    }\
  }\
\
/*  for(int k=0;k<3;k++){\
    iii=ty/(VG_MINIMUM_PARTICLE_BLOCK_J/VG_JDIV);\
    Csub[k]=fi[iii][k];Csub2[k]=fi2[iii][k];\
    }*/\
}


#define MD_KERNEL_CVFB_CI_JDIV(gpu_funcname,gpu_kernel_core,gpu_kernel_core2) \
extern "C" __global__ \
void gpu_funcname(int nii, int nj, int natj,\
		  VG_IVEC *ivec, VG_JVEC *jvec,\
		  VG_FVEC *fvec, VG_PSCALER *pscal)\
{\
  int by = blockIdx.x;\
  int ty = threadIdx.x;\
  int aBegin = VG_MINIMUM_PARTICLE_BLOCK_I*by+ty;;\
  /*int jEnd = (nj+VG_MINIMUM_PARTICLE_BLOCK_J-1)/VG_MINIMUM_PARTICLE_BLOCK_J;*/\
  int jEnd = nj/VG_MINIMUM_PARTICLE_BLOCK_J;\
  int jty,jty4,jty8,ty4,ty8,jj4,jj;					\
  int k,atypei,itype;							\
  float dx, dy, dz, eps=d_scalers->eps2;	\
  float api2=(float)(2.e0/sqrt(M_PI)),kkpi1=1.0/(2.0*M_PI*sqrt(2.0*M_PI)); \
  float rij,sij,pi4=4.0/M_PI,pi14=0.25/M_PI;\
  float l = d_scalers->volume[0];\
  float al = 1.e0 / d_scalers->volume[0];\
  float size[3]={d_scalers->volume[0],d_scalers->volume[1],d_scalers->volume[2]};\
  float r2max=d_scalers->r2max;\
  float rcut21=d_scalers->rcut21;\
  float alpha2=(d_scalers->alpha)*(d_scalers->alpha);\
  float scaleqj_1=d_scalers->scaleqj_1;\
  float alpha3fac=d_scalers->alpha3fac;\
  float dn6;\
  DEFINE_CSUB;\
  /*double Csubd[MD_FVEC_SIZE];*/\
  /*float *Csub=(float *)Csubd,*Csub2=((float *)Csubd)+MD_FVEC_SIZE;*/\
  /*float Csub[MD_FVEC_SIZE],Csub2[MD_FVEC_SIZE];*/\
  /*double *Csubd=NULL;*/\
  float cutoff,kk1=0.5f,pi=(float)M_PI; \
  float tmp,*jvecpoint,halfsigmai,halfsigmaj,sqrtepsiloni,sqrtepsilonj;\
  __shared__ float Xis[VG_MINIMUM_PARTICLE_BLOCK_I*LJ_BS_SIZE];\
  __shared__ float Bs[1][VG_MINIMUM_PARTICLE_BLOCK_I*LJ_BS_SIZE];\
  __shared__ float As[3][VG_MINIMUM_ATYPE_BLOCK];\
  float tmp2,tmp3,tmp4,al2[3];\
  int sameatypei=0;\
  DEFINE_XI;\
  \
  \
  ty4=ty*LJ_BS_SIZE;\
  \
  INITIALIZE_CSUB_VDW;\
  \
  /*jvecpoint=(float *)ivec;*/			\
  /*jty4=aBegin*LJ_BS_SIZE;*/			\
  /*LJ_LOAD_I;*/				\
  /*LJ_LOAD_MATRIX;*/				\
  \
  jvecpoint=(float *)ivec;			\
  jty4=aBegin*LJ_BS_SIZE;			\
  Xis[ty4+0]=jvecpoint[jty4+0];\
  Xis[ty4+1]=jvecpoint[jty4+1];\
  Xis[ty4+2]=jvecpoint[jty4+2];\
  Xis[ty4+3]=jvecpoint[jty4+3];\
  __syncthreads();\
  jvecpoint=(float *)jvec;\
  if(sameatypei==0){\
  /*if(by<10){*/\
  /*if(1){*/\
    gpu_kernel_core(ty4,nj,natj,pscal,atypei,Xis,jvecpoint,Bs,As,al2,size,r2max,rcut21,alpha2,api2,scaleqj_1,alpha3fac,xi,Csub,Csub2,Csubd); \
  }\
  else{\
    gpu_kernel_core2(ty4,nj,natj,pscal,atypei,Xis,jvecpoint,Bs,As,al2,size,r2max,rcut21,alpha2,api2,scaleqj_1,alpha3fac,xi,Csub,Csub2,Csubd); \
  }\
  WRITE_BACK_TO_FVEC;\
}
#endif // end of VG_MINIMUM_PARTICLE_BLOCK_J==64 && VG_JDIV==4 && MD_UNROLL_KERNEL==8

#else // else of MD_USE_QAUNION
#define MD_KERNEL(gpu_funcname,gpu_kernel_core) \
extern "C" __global__ \
void gpu_funcname(int nii, int nj, int natj,\
		  VG_IVEC *ivec, VG_JVEC *jvec,\
		  VG_FVEC *fvec, VG_PSCALER *pscal)\
{\
  int by = blockIdx.x;\
  int ty = threadIdx.x;\
  int aBegin = VG_MINIMUM_PARTICLE_BLOCK_I*by+ty;\
  int jEnd = (nj+VG_MINIMUM_PARTICLE_BLOCK_J-1)/VG_MINIMUM_PARTICLE_BLOCK_J;\
  int jty;\
  int itype,k;\
  \
  float l = d_scalers->volume[0];\
  float al = 1.e0 / d_scalers->volume[0];\
  float dn2, dn6, alpha2=(d_scalers->alpha)*(d_scalers->alpha);	\
  float dx, dy, dz, tmp, eps2=d_scalers->eps2;\
  float api2=(float)(2.e0/sqrt(M_PI));\
  \
  float Csub[3];\
  \
  __shared__ float As[3][VG_MINIMUM_PARTICLE_BLOCK_I];\
  __shared__ int Ds[1][VG_MINIMUM_PARTICLE_BLOCK_I];\
  __shared__ float Bs[3][VG_MINIMUM_PARTICLE_BLOCK_I];\
  __shared__ int Es[1][VG_MINIMUM_PARTICLE_BLOCK_I];\
  __shared__ float Qis[1][VG_MINIMUM_PARTICLE_BLOCK_I];\
  __shared__ float Qjs[1][VG_MINIMUM_PARTICLE_BLOCK_I];\
  \
  Csub[0] = 0.e0;\
  Csub[1] = 0.e0;\
  Csub[2] = 0.e0;\
  \
  for(k=0;k<3;k++) AS(k,ty)=ivec[aBegin].r[k];\
  DS(0,ty)=ivec[aBegin].atype;\
  QIS(0,ty)=ivec[aBegin].q;\
  \
  for (int j = 0; j < jEnd; j++){\
  \
    jty=j*VG_MINIMUM_PARTICLE_BLOCK_I+ty;\
    for(k=0;k<3;k++) BS(k,ty)=jvec[jty].r[k];\
    ES(0,ty)=jvec[jty].atype;\
    QJS(0,ty)=jvec[jty].q;\
    \
    __syncthreads();\
    for (int jj = 0; jj < VG_MINIMUM_PARTICLE_BLOCK_J; jj++){\
      dx = AS(0,ty) - BS(0,jj);\
      dy = AS(1,ty) - BS(1,jj);\
      dz = AS(2,ty) - BS(2,jj);\
      dx = dx - rintf(dx * al) * l;\
      dy = dy - rintf(dy * al) * l;\
      dz = dz - rintf(dz * al) * l;\
      itype = natj * DS(0,ty) + ES(0,jj);\
      gpu_kernel_core;\
    }\
    __syncthreads();\
  }\
  for(k=0;k<3;k++) fvec[aBegin].fi[k]=Csub[k];\
}
#endif // end of MD_USE_QAUNION


MD_KERNEL(lj_kernel_gpu,LJ_KERNEL_CORE);
MD_KERNEL(ljpot_kernel_gpu,LJPOT_KERNEL_CORE);
MD_KERNEL(grav_kernel_gpu,GRAV_KERNEL_CORE);
MD_KERNEL(gravpot_kernel_gpu,GRAVPOT_KERNEL_CORE);
MD_KERNEL(real_kernel_gpu,REAL_KERNEL_CORE);
MD_KERNEL(realpot_kernel_gpu,REALPOT_KERNEL_CORE);
MD_KERNEL(wave_kernel_gpu,WAVE_KERNEL_CORE);
MD_KERNEL(wavepot_kernel_gpu,WAVEPOT_KERNEL_CORE);
MD_KERNEL(r1_kernel_gpu,R1_KERNEL_CORE);
MD_KERNEL(rsqrt_kernel_gpu,RSQRT_KERNEL_CORE);
MD_KERNEL(mul_kernel_gpu,MUL_KERNEL_CORE);
#ifdef MD_NACL
MD_KERNEL(nacl_kernel_gpu,NACL_KERNEL_CORE);
#endif

// vdw force
//MD_KERNEL_FB_SUB(ljfb_kernel_gpu_sub,LJFB_KERNEL_CORE,LJ_LOAD_J,0);
MD_KERNEL_FB_SUB(ljfb_kernel_gpu_sub,LJFB_KERNEL_CORE,LJ_LOAD_J_ATYPEBIT,0);
//MD_KERNEL_FB_SUB(ljfb_kernel_gpu_sub_rounduploop,LJFB_KERNEL_CORE,LJ_LOAD_J,1);
MD_KERNEL_FB_SUB(ljfb_kernel_gpu_sub_rounduploop,LJFB_KERNEL_CORE,LJ_LOAD_J_ATYPEBIT,1);
#ifdef MD_SORT_ATYPEI
MD_KERNEL_FB_SUB(ljfb_sameatypei_kernel_gpu_sub,LJFB_KERNEL_CORE2,LJ_LOAD_J2,0);
MD_KERNEL_FB_SUB(ljfb_sameatypei_kernel_gpu_sub_rounduploop,LJFB_KERNEL_CORE2,LJ_LOAD_J2,1);
MD_KERNEL_FBLJ(ljfb_kernel_gpu,ljfb_kernel_gpu_sub,ljfb_sameatypei_kernel_gpu_sub);
MD_KERNEL_FBLJ(ljfb_kernel_gpu_rounduploop,ljfb_kernel_gpu_sub_rounduploop,ljfb_sameatypei_kernel_gpu_sub_rounduploop);
#else
MD_KERNEL_FBLJ(ljfb_kernel_gpu,ljfb_kernel_gpu_sub,ljfb_kernel_gpu_sub);
MD_KERNEL_FBLJ(ljfb_kernel_gpu_rounduploop,ljfb_kernel_gpu_sub_rounduploop,ljfb_kernel_gpu_sub);
#endif

// vdw pot
//MD_KERNEL_FB_SUB(ljpotfb_kernel_gpu_sub,LJPOTFB_KERNEL_CORE,LJ_LOAD_J,0);
MD_KERNEL_FB_SUB(ljpotfb_kernel_gpu_sub,LJPOTFB_KERNEL_CORE,LJ_LOAD_J_ATYPEBIT,0);
//MD_KERNEL_FB_SUB(ljpotfb_kernel_gpu_sub_rounduploop,LJPOTFB_KERNEL_CORE,LJ_LOAD_J,1);
MD_KERNEL_FB_SUB(ljpotfb_kernel_gpu_sub_rounduploop,LJPOTFB_KERNEL_CORE,LJ_LOAD_J_ATYPEBIT,1);
#ifdef MD_SORT_ATYPEI
MD_KERNEL_FB_SUB(ljpotfb_sameatypei_kernel_gpu_sub,LJPOTFB_KERNEL_CORE2,LJ_LOAD_J2,0);
MD_KERNEL_FB_SUB(ljpotfb_sameatypei_kernel_gpu_sub_rounduploop,LJPOTFB_KERNEL_CORE2,LJ_LOAD_J2,1);
MD_KERNEL_FBLJ(ljpotfb_kernel_gpu,ljpotfb_kernel_gpu_sub,ljpotfb_sameatypei_kernel_gpu_sub);
MD_KERNEL_FBLJ(ljpotfb_kernel_gpu_rounduploop,ljpotfb_kernel_gpu_sub_rounduploop,ljpotfb_sameatypei_kernel_gpu_sub_rounduploop);
#else
MD_KERNEL_FBLJ(ljpotfb_kernel_gpu,ljpotfb_kernel_gpu_sub,ljpotfb_kernel_gpu_sub);
MD_KERNEL_FBLJ(ljpotfb_kernel_gpu_rounduploop,ljpotfb_kernel_gpu_sub_rounduploop,ljpotfb_kernel_gpu_sub);
#endif

// coulomb force
MD_KERNEL_FB(gravfb_kernel_gpu,GRAVFB_KERNEL_CORE,0);
MD_KERNEL_FB(gravfb_kernel_gpu_rounduploop,GRAVFB_KERNEL_CORE,1);

// coulomb potential
MD_KERNEL_FB(gravpotfb_kernel_gpu,GRAVPOTFB_KERNEL_CORE,0);
MD_KERNEL_FB(gravpotfb_kernel_gpu_rounduploop,GRAVPOTFB_KERNEL_CORE,1);

// real force
MD_KERNEL_FB(realfb_kernel_gpu,REALFB_KERNEL_CORE,0);
MD_KERNEL_FB(realfb_kernel_gpu_rounduploop,REALFB_KERNEL_CORE,1);

// real potential
MD_KERNEL_FB(realpotfb_kernel_gpu,REALPOTFB_KERNEL_CORE,0);
MD_KERNEL_FB(realpotfb_kernel_gpu_rounduploop,REALPOTFB_KERNEL_CORE,1);


#ifdef MD_CELLINDEX
#ifdef MD_SORT_ATYPEI
#else
MD_KERNEL_FB_SUB_CI(ljfbci_kernel_gpu_sub,LJFB_KERNEL_CORE,LJ_LOAD_J);
//MD_KERNEL_FB_SUB_CI(ljfbci_kernel_gpu_sub,REALLJFB_KERNEL_CORE,LJ_LOAD_J);
MD_KERNEL_FB_SUB_CI(ljpotfbci_kernel_gpu_sub,LJPOTFB_KERNEL_CORE,LJ_LOAD_J);
MD_KERNEL_FBLJ_CI(ljfbci_kernel_gpu,ljfbci_kernel_gpu_sub,ljfbci_kernel_gpu_sub);
MD_KERNEL_FBLJ_CI(ljpotfbci_kernel_gpu,ljpotfbci_kernel_gpu_sub,ljpotfbci_kernel_gpu_sub);

#if VG_JDIV==1
MD_KERNEL_CVFB_SUB_CI(realljfbci_kernel_gpu_sub,REALLJFB_KERNEL_CORE,LJ_LOAD_J);
MD_KERNEL_CVFB_CI(realljfbci_kernel_gpu,realljfbci_kernel_gpu_sub,realljfbci_kernel_gpu_sub);
#else
MD_KERNEL_CVFB_SUB_CI_JDIV(realljfbci_kernel_gpu_sub,REALLJFB_KERNEL_CORE,LJ_LOAD_J);
MD_KERNEL_CVFB_CI_JDIV(realljfbci_kernel_gpu,realljfbci_kernel_gpu_sub,realljfbci_kernel_gpu_sub);
#endif

//MD_KERNEL_CI(realljfbci_kernel_gpu,REAL_KERNEL_CORE);

#if VG_JDIV==1
MD_KERNEL_CVFB_SUB_CI(realljpotfbci_kernel_gpu_sub,REALLJPOTFB_KERNEL_CORE,LJ_LOAD_J);
MD_KERNEL_CVFB_CI(realljpotfbci_kernel_gpu,realljpotfbci_kernel_gpu_sub,realljpotfbci_kernel_gpu_sub);
#else
MD_KERNEL_CVFB_SUB_CI_JDIV(realljpotfbci_kernel_gpu_sub,REALLJPOTFB_KERNEL_CORE,LJ_LOAD_J);
MD_KERNEL_CVFB_CI_JDIV(realljpotfbci_kernel_gpu,realljpotfbci_kernel_gpu_sub,realljpotfbci_kernel_gpu_sub);
#endif
MD_KERNEL_FB_SUB_CI(realljpotfbci_nojdiv_kernel_gpu_sub,REALLJPOTFB_NOJDIV_KERNEL_CORE,LJ_LOAD_J);
MD_KERNEL_FBLJ_CI(realljpotfbci_nojdiv_kernel_gpu,realljpotfbci_nojdiv_kernel_gpu_sub,realljpotfbci_nojdiv_kernel_gpu_sub);

#endif
MD_KERNEL_CI(realci_kernel_gpu,REAL_KERNEL_CORE);
MD_KERNEL_CI(realpotci_kernel_gpu,REALPOT_KERNEL_CORE);
MD_KERNEL_FB_CI(realfbci_kernel_gpu,REALFB_KERNEL_CORE);
MD_KERNEL_FB_CI(realpotfbci_kernel_gpu,REALPOTFB_KERNEL_CORE);
#endif



#else

#ifdef MD_USE_CONSTANT
extern "C"
__global__ 
void lj_kernel_gpu(int nii, int nj, int natj, 
		   VG_IVEC *ivec, VG_JVEC *jvec, 
		   VG_FVEC *fvec, VG_PSCALER *pscal)
{
  int by = blockIdx.x;
  int ty = threadIdx.x;
  int aBegin = VG_MINIMUM_PARTICLE_BLOCK_I*by+ty;
  int jEnd = (nj+VG_MINIMUM_PARTICLE_BLOCK_J-1)/VG_MINIMUM_PARTICLE_BLOCK_J;
  int ty3, j3, jty;
  int itype,k;

  float l = d_scalers->volume[0];
  float al = 1.e0 / d_scalers->volume[0];
  float dn2, dn6;
  float dx, dy, dz, tmp;

  float Csub[3];

  __shared__ float As[1][VG_MINIMUM_PARTICLE_BLOCK_I * 3];
  __shared__ int Ds[1][VG_MINIMUM_PARTICLE_BLOCK_I];
  __shared__ float Bs[1][VG_MINIMUM_PARTICLE_BLOCK_I * 3];
  __shared__ int Es[1][VG_MINIMUM_PARTICLE_BLOCK_I];

  ty3 = ty * 3;

  Csub[0] = 0.e0;
  Csub[1] = 0.e0;
  Csub[2] = 0.e0;

  for(k=0;k<3;k++) AS(0,ty3+k)=ivec[aBegin].r[k];
  DS(0,ty)=ivec[aBegin].atype;

  for (int j = 0; j < jEnd; j++){

    jty=j*VG_MINIMUM_PARTICLE_BLOCK_I+ty;
    for(k=0;k<3;k++) BS(0,ty3+k)=jvec[jty].r[k];
    ES(0,ty)=jvec[jty].atype;

    __syncthreads();
    for (int jj = 0; jj < VG_MINIMUM_PARTICLE_BLOCK_J; jj++){
      j3 = jj*3;
      dx = AS(0,ty3)   - BS(0,j3);
      dy = AS(0,ty3+1) - BS(0,j3+1);
      dz = AS(0,ty3+2) - BS(0,j3+2);
      dx = dx - rintf(dx * al) * l;
      dy = dy - rintf(dy * al) * l;
      dz = dz - rintf(dz * al) * l;
      itype = natj * DS(0,ty) + ES(0,jj);
      dn2 = (dx * dx + dy * dy + dz * dz) * d_matrix[itype].rscale;
      if(j*VG_MINIMUM_PARTICLE_BLOCK_I+jj<nj && dn2 >= MD_LJ_R2MIN && dn2<MD_LJ_R2MAX && dn2!=0.0f){
	dn6 = 1.e0 / (dn2 * dn2 * dn2);
	tmp = d_matrix[itype].gscale * dn6 / dn2 * (2.e0 * dn6 - 1.e0);
	Csub[0] += tmp * dx;
	Csub[1] += tmp * dy;
	Csub[2] += tmp * dz;
      }
    }
    __syncthreads();
  }
  for(k=0;k<3;k++) fvec[aBegin].fi[k]=Csub[k];
}
#else // else of MD_USE_CONSTANT
extern "C"
__global__ 
void lj_kernel_gpu(int nii, int nj, int natj, 
		   VG_IVEC *ivec, VG_JVEC *jvec, 
		   VG_MATRIX *matrix, VG_SCALER *scalers,
		   VG_FVEC *fvec, VG_PSCALER *pscal)
{
  int by = blockIdx.x;
  int ty = threadIdx.x;
  int aBegin = VG_MINIMUM_PARTICLE_BLOCK_I*by+ty;
  int jEnd = (nj+VG_MINIMUM_PARTICLE_BLOCK_J-1)/VG_MINIMUM_PARTICLE_BLOCK_J;
  int ty3, j3, jty;
  int itype,k;

  float l = scalers->volume[0];
  float al = 1.e0 / scalers->volume[0];
  float dn2, dn6;
  float dx, dy, dz, tmp;

  float Csub[3];

  __shared__ float As[1][VG_MINIMUM_PARTICLE_BLOCK_I * 3];
  __shared__ float Fs[1][VG_MINIMUM_PARTICLE_BLOCK_I];
  __shared__ float Gs[1][VG_MINIMUM_PARTICLE_BLOCK_I];
  __shared__ int Ds[1][VG_MINIMUM_PARTICLE_BLOCK_I];
  __shared__ float Bs[1][VG_MINIMUM_PARTICLE_BLOCK_I * 3];
  __shared__ int Es[1][VG_MINIMUM_PARTICLE_BLOCK_I];

  ty3 = ty * 3;

  Csub[0] = 0.e0;
  Csub[1] = 0.e0;
  Csub[2] = 0.e0;

  FS(0,ty)=matrix[ty].rscale;
  GS(0,ty)=matrix[ty].gscale;

  for(k=0;k<3;k++) AS(0,ty3+k)=ivec[aBegin].r[k];
  DS(0,ty)=ivec[aBegin].atype;

  for (int j = 0; j < jEnd; j++){

    jty=j*VG_MINIMUM_PARTICLE_BLOCK_I+ty;
    for(k=0;k<3;k++) BS(0,ty3+k)=jvec[jty].r[k];
    ES(0,ty)=jvec[jty].atype;

    __syncthreads();
    for (int jj = 0; jj < VG_MINIMUM_PARTICLE_BLOCK_J; jj++){
      j3 = jj*3;
      dx = AS(0,ty3)   - BS(0,j3);
      dy = AS(0,ty3+1) - BS(0,j3+1);
      dz = AS(0,ty3+2) - BS(0,j3+2);
      dx = dx - rintf(dx * al) * l;
      dy = dy - rintf(dy * al) * l;
      dz = dz - rintf(dz * al) * l;
      itype = natj * DS(0,ty) + ES(0,jj);
      dn2 = (dx * dx + dy * dy + dz * dz) * FS(0,itype);
      if(j*VG_MINIMUM_PARTICLE_BLOCK_I+jj<nj && dn2 >= MD_LJ_R2MIN && dn2<MD_LJ_R2MAX && dn2!=0.0f){
      //      if(dn2!=0.0f){
	dn6 = 1.e0 / (dn2 * dn2 * dn2);
	tmp = GS(0,itype) * dn6 / dn2 * (2.e0 * dn6 - 1.e0);
	Csub[0] += tmp * dx;
	Csub[1] += tmp * dy;
	Csub[2] += tmp * dz;
      }
    }
    __syncthreads();
  }
  for(k=0;k<3;k++) fvec[aBegin].fi[k]=Csub[k];
}
#endif // end of MD_USE_CONSTANT
#endif 
