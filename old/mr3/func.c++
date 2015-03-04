#include <pthread.h>
#include "md.h"

#ifdef GPUCODE
#define SUB_HEADER __device__
#else
#define SUB_HEADER 
#endif

template <int SHIFT_FUNC, class DF>
__inline__ SUB_HEADER
DF calc_vdw(DF r1, DF r2, DF r6, 
	    DF shift, DF rcut21, 
	    DF gs, DF factorv)
{
  DF dtmp;
  DF one=(DF)1.0,two=(DF)2.0,three=(DF)3.0;

  switch(SHIFT_FUNC){
  case 0:     // normal Lennard Jones
    dtmp=gs*r6*r1*(two*r6-one);
    break;
  case 100:
    dtmp=gs*r6*(r6-one);
    break;
  case 1:     // VSHIFT
    dtmp=r1;
#ifdef VDW_SHIFT_CORRECT
    dtmp*=r6*(two*r6-one)-shift*(three*shift-two);
#else
    r2*=rcut21;
    r2*=r2*r2;
    dtmp*=r6*(two*r6-one)-r2*shift*(two*shift-one);
#endif
    dtmp*=gs;
    break;
  case 101:
    dtmp=r6*(r6-one);
    r2*=rcut21;
    r2*=r2*r2;
    dtmp+=r2*shift*(two*shift-one)+shift*(two-three*shift);
    dtmp*=gs;
    break;
  case 2:     // VSWITCH
    break;
  case 202:
    break;
  }
  dtmp*=factorv;
  
  return dtmp;
}


#ifndef GPUCODE

extern "C"
float calc_dtmp_vdw_force(float r1, float r2, float r6, 
#ifdef VDW_SHIFT
			       float shift, float rcut21, 
#endif
			       float gs, float factorv)
{
  float dtmp;

#ifdef VDW_SHIFT
  dtmp=r1;
  r1=shift;
#ifdef VDW_SHIFT_CORRECT
  dtmp*=r6*(2.0f*r6-1.0f)-r1*(3.0f*r1-2.0f);
#else
  r2*=rcut21;
  r2*=r2*r2;
  dtmp*=r6*(2.0f*r6-1.0f)-r2*r1*(2.0f*r1-1.0f);
#endif
  dtmp*=gs;
#else
  dtmp=gs*r6*r1*(2.0f*r6-1.0f);
#endif
  dtmp*=factorv;
  
  return dtmp;
}


extern "C"
float calc_dtmp_vdw_force_notwork3(float r1, float r2, float r6, 
			  float shift, float rcut21, 
			  float gs, float factorv)
{
  float dtmp;

  dtmp=r1;
  r1=shift;
#ifdef VDW_SHIFT_CORRECT
  dtmp*=r6*(2.0f*r6-1.0f)-r1*(3.0f*r1-2.0f);
#else
  r2*=rcut21;
  r2*=r2*r2;
  dtmp*=r6*(2.0f*r6-1.0f)-r2*r1*(2.0f*r1-1.0f);
#endif
  dtmp*=gs;
  dtmp*=factorv;
  
  return dtmp;
}


extern "C"
float calc_dtmp_vdw_force_notwork2(float r1, float r2, float r6, 
			  float shift, float rcut21, 
			  float gs, float factorv)
{
  float dtmp;
  float one=(float)1.0,two=(float)2.0,three=(float)3.0;

  dtmp=r1;
#ifdef VDW_SHIFT_CORRECT
  dtmp*=r6*(two*r6-one)-shift*(three*shift-two);
#else
  r2*=rcut21;
  r2*=r2*r2;
  dtmp*=r6*(two*r6-one)-r2*shift*(two*shift-one);
#endif
  dtmp*=gs;
  dtmp*=factorv;
  
  return dtmp;
}


extern "C"
float calc_dtmp_vdw_force_notwork(float r1, float r2, float r6, 
#ifdef VDW_SHIFT
				     float shift, float rcut21, 
#endif
				     float gs, float factorv)
{
#ifdef VDW_SHIFT
  return calc_vdw<VDW_SHIFT,float>(r1,r2,r6,shift,rcut21,gs,factorv);
#else
  float shift, rcut21;
  return calc_vdw<0,float>(r1,r2,r6,shift,rcut21,gs,factorv);
#endif
}


extern "C"
float calc_dtmp_vdw_pot(float r1, float r2, float r6, 
#ifdef VDW_SHIFT
				     float shift, float rcut21, 
#endif
				     float gs, float factorv)
{
#ifdef VDW_SHIFT
  return calc_vdw<VDW_SHIFT+100,float>(r1,r2,r6,shift,rcut21,gs,factorv);
#else
  float shift, rcut21;
  return calc_vdw<100,float>(r1,r2,r6,shift,rcut21,gs,factorv);
#endif
}


#endif // end of not GPUCODE

