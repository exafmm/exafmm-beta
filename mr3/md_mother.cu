#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

#if 0 // use cuda util in SDK
#include <cutil.h> 
#else             
#define CUDA_SAFE_CALL(x) x
#define CUT_CHECK_ERROR(x) 
#endif                  

//#define CHECK_CUDA_ERROR // check cuda error after cuda calls
#define MD_USE_CONSTANT

#include "md.h"
#include "md_kernel.cu"


#ifdef CUDA_SDK_2
//extern int Argc;
//extern char **Argv;
//#define MY_CUT_DEVICE_INIT(argc,argv) CUT_DEVICE_INIT(argc,argv)
#define MY_CUT_DEVICE_INIT(argc,argv) 
#else
//static int Argc=0;
//static char **argv=NULL;
#define MY_CUT_DEVICE_INIT(argc,argv) CUT_DEVICE_INIT()
#endif

#ifdef CHECK_CUDA_ERROR // error check
#define CHECK_ERROR(x) CudaError(x,1);
#else
#define CHECK_ERROR(x) 
#endif


static int CudaError(char *message, int flag)
{
  /*
    return 0 -- no error
           1 -- error
    flag   0 -- do not display message
           1 -- display message
   */
  int ret=0;
  cudaError_t err;

#ifdef CHECK_CUDA_ERROR
  {
    static int ini=0;
    if(ini==0){
      printf("** CudaError is called because CHECK_ERROR is defined **\n");
      ini=1;
    }
  }
#endif
  err=cudaGetLastError();
  if(cudaSuccess!=err){
    if(flag) fprintf(stderr,"** CUDA error occurred : %s **\n",message);
    ret=1;
  }
  else{
    if(flag) fprintf(stderr,"%s\n",message);
  }
  return ret;
}

static void mycudaMalloc(void **p, long size, char *s)
{
#ifdef CHECK_CUDA_ERROR
  printf("before cudaMalloc : size=%d %s\n",size,s);
#endif
  CUDA_SAFE_CALL(cudaMalloc(p,size));
#ifdef CHECK_CUDA_ERROR
  CHECK_ERROR("After cudaMalloc");
  printf("Cuda Allocated %d memory at %016llx : %s\n",size,*p,s);
#endif
}


static __inline__ void mycudaFree(void **p, char *s)
{
  CUDA_SAFE_CALL(cudaFree(*p));
#ifdef CHECK_CUDA_ERROR
  CHECK_ERROR("After cudaFree");
  printf("Cuda Freeed memory at %016llx : %s\n",*p,s);
  *p=NULL;
#endif
}


void mySetDevice(int deviceid)
{
  static int ini=0;

  if(ini==0){
#ifdef MD_PRINT_WARN
    printf("mySetDevice is called (deviceid=%d)\n",deviceid);fflush(stdout);
#endif
    CUDA_SAFE_CALL(cudaSetDevice(deviceid));
    CudaError("after cudaSetDevice",0);
    cudaDeviceProp deviceProp;
    cudaGetDeviceProperties(&deviceProp, deviceid);
    printf("Device %d: %s\n", deviceid, deviceProp.name);
    ini=1;
  }
}


#define MD_MOTHER(mother_funcname,gpu_funcname,gpu_funcname2) \
extern "C" \
void mother_funcname(void *unit_org)\
{\
  VG_UNIT *unit=(VG_UNIT *)unit_org;\
  int ni=unit->ni;\
  int nj=unit->nj;\
  int nf=unit->nf;\
  int natj=unit->natj;\
  VG_IVEC *ivec=(VG_IVEC *)unit->ivectors;\
  VG_JVEC *jvec=(VG_JVEC *)unit->jvectors;\
  VG_MATRIX *matrix=(VG_MATRIX *)unit->matrices;\
  VG_SCALER *scalers=(VG_SCALER *)unit->scalers;\
  VG_FVEC *fvec=(VG_FVEC *)unit->fvectors;\
  VG_PSCALER *pscal=(VG_PSCALER *)unit->pscalers;\
  unsigned int mem_size_i = sizeof(VG_IVEC)*NI_ROUNDUP(ni);\
  unsigned int mem_size_j = sizeof(VG_JVEC)*NJ_ROUNDUP(nj);\
  unsigned int mem_size_matrix = sizeof(VG_MATRIX)*VG_MINIMUM_ATYPE_BLOCK;\
  unsigned int mem_size_scalers = sizeof(VG_SCALER);\
  unsigned int mem_size_f = sizeof(VG_FVEC)*NF_ROUNDUP(nf);\
  unsigned int mem_size_p = sizeof(VG_PSCALER);\
  VG_IVEC *d_ivec;\
  VG_JVEC *d_jvec;\
  VG_FVEC *d_fvec;\
  VG_PSCALER *d_pscal;\
  int s_gpuCount=0;\
  \
  MY_CUT_DEVICE_INIT(Argc,Argv);				\
  CHECK_ERROR("After MY_CUT_DEVICE_INIT");\
  CUDA_SAFE_CALL(cudaGetDeviceCount(&s_gpuCount));\
  CHECK_ERROR("After cudaGetDeviceCount");\
  if(unit->deviceid<s_gpuCount){\
    mySetDevice(unit->deviceid);\
  }\
  else{\
    fprintf(stderr,"** error : no such device id=%d **\n",unit->deviceid);\
    exit(1);\
  }\
  CHECK_ERROR("After cudaSetDevice");\
  mycudaMalloc((void**)&d_ivec,mem_size_i,"d_ivec MD_MOTHER");			\
  mycudaMalloc((void**)&d_jvec,mem_size_j,"d_jvec MD_MOTHER");				\
  mycudaMalloc((void**)&d_fvec,mem_size_f,"d_fvec MD_MOTHER");				\
  mycudaMalloc((void**)&d_pscal,mem_size_p,"d_pscal MD_MOTHER");				\
  \
  CUDA_SAFE_CALL(cudaMemcpy(d_ivec,ivec,mem_size_i,cudaMemcpyHostToDevice));\
  CHECK_ERROR("After cudaMemcpy d_ivec");\
  CUDA_SAFE_CALL(cudaMemcpy(d_jvec,jvec,mem_size_j,cudaMemcpyHostToDevice));\
  CHECK_ERROR("After cudaMemcpy d_jvec");\
  if(matrix!=NULL) CUDA_SAFE_CALL(cudaMemcpyToSymbol(d_matrix,matrix,mem_size_matrix));\
  CHECK_ERROR("After cudaMemcpy d_matrix");\
  CUDA_SAFE_CALL(cudaMemcpyToSymbol(d_scalers,scalers,mem_size_scalers));\
  CHECK_ERROR("After cudaMemcpy d_scalers");\
  CUDA_SAFE_CALL(cudaMemcpy(d_fvec,fvec,mem_size_f,cudaMemcpyHostToDevice));\
  CHECK_ERROR("After cudaMemcpy d_fvec");\
  CUDA_SAFE_CALL(cudaMemcpy(d_pscal,pscal,mem_size_p,cudaMemcpyHostToDevice));\
  CHECK_ERROR("After cudaMemcpy d_pscal");\
  \
  dim3 threads(VG_MINIMUM_PARTICLE_BLOCK_I);\
  dim3 grid((ni+VG_MINIMUM_PARTICLE_BLOCK_I-1)/VG_MINIMUM_PARTICLE_BLOCK_I);\
  if((nj % VG_MINIMUM_PARTICLE_BLOCK_J)!=0){\
    gpu_funcname<<< grid, threads >>>(ni,nj,natj,d_ivec,d_jvec,		\
				    d_fvec,d_pscal);			\
    CHECK_ERROR("After kernel gpu_funcname");\
  }\
  else{\
    gpu_funcname2<<< grid, threads >>>(ni,nj,natj,d_ivec,d_jvec,		\
				    d_fvec,d_pscal);			\
    CHECK_ERROR("After kernel gpu_funcname2");\
  }\
  CUT_CHECK_ERROR("Kernel execution failed");\
  \
  CUDA_SAFE_CALL(cudaMemcpy(fvec,d_fvec,mem_size_f,cudaMemcpyDeviceToHost));\
  CUDA_SAFE_CALL(cudaMemcpy(pscal,d_pscal,mem_size_p,cudaMemcpyDeviceToHost));\
  CHECK_ERROR("After cudaMemcpy fvec, pscal");\
  \
  mycudaFree((void **)&d_ivec,"d_ivec in MD_MOTHER");		\
  mycudaFree((void **)&d_jvec,"d_jvec in MD_MOTHER");				\
  mycudaFree((void **)&d_fvec,"d_fvec in MD_MOTHER");					\
  mycudaFree((void **)&d_pscal,"d_pscal in MD_MOTHER");					\
}


#define MD_MOTHER_CI(mother_funcname,gpu_funcname) \
extern "C" \
void mother_funcname(void *unit_org)\
{\
  VG_UNIT *unit=(VG_UNIT *)unit_org;\
  int ni=unit->ni;\
  int nj=unit->nj;\
  int nf=unit->nf;\
  int natj=unit->natj;\
  VG_IVEC *ivec=(VG_IVEC *)unit->ivectors;\
  VG_JVEC *jvec=(VG_JVEC *)unit->jvectors;\
  VG_MATRIX *matrix=(VG_MATRIX *)unit->matrices;\
  VG_SCALER *scalers=(VG_SCALER *)unit->scalers;\
  VG_FVEC *fvec=(VG_FVEC *)unit->fvectors;\
  VG_PSCALER *pscal=(VG_PSCALER *)unit->pscalers;\
  unsigned int mem_size_i = sizeof(VG_IVEC)*NI_ROUNDUP(ni);\
  unsigned int mem_size_j = sizeof(VG_JVEC)*NJ_ROUNDUP(nj);\
  unsigned int mem_size_matrix = sizeof(VG_MATRIX)*VG_MINIMUM_ATYPE_BLOCK;\
  unsigned int mem_size_scalers = sizeof(VG_SCALER);\
  unsigned int mem_size_f = sizeof(VG_FVEC)*NF_ROUNDUP(nf);\
  unsigned int mem_size_p = sizeof(VG_PSCALER);\
  static VG_IVEC *d_ivec;\
  static VG_JVEC *d_jvec;\
  static VG_FVEC *d_fvec;\
  static VG_PSCALER *d_pscal;\
  int s_gpuCount=0;\
  int gpuoverlapflag=unit->gpuoverlapflag;\
  \
  if(gpuoverlapflag==0 || gpuoverlapflag==1){\
    MY_CUT_DEVICE_INIT(Argc,Argv);				\
    CHECK_ERROR("After MY_CUT_DEVICE_INIT");\
    CUDA_SAFE_CALL(cudaGetDeviceCount(&s_gpuCount));\
    CHECK_ERROR("After cudaGetDeviceCount");\
    if(unit->deviceid<s_gpuCount){\
      mySetDevice(unit->deviceid);\
    }\
    else{\
      fprintf(stderr,"** error : no such device id=%d **\n",unit->deviceid);\
      exit(1);\
    }\
    CHECK_ERROR("After cudaSetDevice");\
    mycudaMalloc((void**)&d_ivec,mem_size_i,"d_ivec MD_MOTHER_CI");	\
    mycudaMalloc((void**)&d_jvec,mem_size_j,"d_jvec MD_MOTHER_CI");				\
    mycudaMalloc((void**)&d_fvec,mem_size_f,"d_fvec MD_MOTHER_CI");\
    mycudaMalloc((void**)&d_pscal,mem_size_p,"d_pscal MD_MOTHER_CI");			   \
    \
    CUDA_SAFE_CALL(cudaMemcpy(d_ivec,ivec,mem_size_i,cudaMemcpyHostToDevice));\
    CHECK_ERROR("After cudaMemcpy d_ivec");\
    CUDA_SAFE_CALL(cudaMemcpy(d_jvec,jvec,mem_size_j,cudaMemcpyHostToDevice));\
    CHECK_ERROR("After cudaMemcpy d_jvec");\
    if(matrix!=NULL) CUDA_SAFE_CALL(cudaMemcpyToSymbol(d_matrix,matrix,mem_size_matrix));\
    CHECK_ERROR("After cudaMemcpy d_matrix");\
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(d_scalers,scalers,mem_size_scalers));\
    CHECK_ERROR("After cudaMemcpy d_scalers");\
    CUDA_SAFE_CALL(cudaMemcpy(d_fvec,fvec,mem_size_f,cudaMemcpyHostToDevice));\
    CHECK_ERROR("After cudaMemcpy d_fvec");\
    CUDA_SAFE_CALL(cudaMemcpy(d_pscal,pscal,mem_size_p,cudaMemcpyHostToDevice));\
    CHECK_ERROR("After cudaMemcpy d_pscal");\
    \
    dim3 threads(VG_MINIMUM_PARTICLE_BLOCK_I);\
    dim3 grid((ni+VG_MINIMUM_PARTICLE_BLOCK_I-1)/VG_MINIMUM_PARTICLE_BLOCK_I);\
    if(1){\
    /*if((nj % VG_MINIMUM_PARTICLE_BLOCK_J)!=0){*/\
      gpu_funcname<<< grid, threads >>>(ni,nj,natj,d_ivec,d_jvec,		\
	  			        d_fvec,d_pscal);			\
      CHECK_ERROR("After kernel gpu_funcname");\
    }\
    else{\
      printf("nj=%d is not multiple of %d, ni=%d nf=%d\n",nj,VG_MINIMUM_PARTICLE_BLOCK_J,ni,nf);\
    }\
    CUT_CHECK_ERROR("Kernel execution failed");\
  }\
  \
  if(gpuoverlapflag==0 || gpuoverlapflag==2){\
    CUDA_SAFE_CALL(cudaMemcpy(fvec,d_fvec,mem_size_f,cudaMemcpyDeviceToHost));\
    CUDA_SAFE_CALL(cudaMemcpy(pscal,d_pscal,mem_size_p,cudaMemcpyDeviceToHost));\
    CHECK_ERROR("After cudaMemcpy fvec, pscal");\
    \
    mycudaFree((void **)&d_ivec,"d_ivec in MD_MOTHER_CI");		\
    mycudaFree((void **)&d_jvec,"d_jvec in MD_MOTHER_CI");		\
    mycudaFree((void **)&d_fvec,"d_fvec in MD_MOTHER_CI");		\
    mycudaFree((void **)&d_pscal,"d_pscal in MD_MOTHER_CI");		\
  }\
}


#if 1
MD_MOTHER(lj_mother,lj_kernel_gpu,lj_kernel_gpu);
MD_MOTHER(ljpot_mother,ljpot_kernel_gpu,ljpot_kernel_gpu);
MD_MOTHER(grav_mother,grav_kernel_gpu,grav_kernel_gpu);
MD_MOTHER(gravpot_mother,gravpot_kernel_gpu,gravpot_kernel_gpu);
MD_MOTHER(real_mother,real_kernel_gpu,real_kernel_gpu);
MD_MOTHER(realpot_mother,realpot_kernel_gpu,realpot_kernel_gpu);
MD_MOTHER(wave_mother,wave_kernel_gpu,wave_kernel_gpu);
MD_MOTHER(wavepot_mother,wavepot_kernel_gpu,wavepot_kernel_gpu);
MD_MOTHER(r1_mother,r1_kernel_gpu,r1_kernel_gpu);
MD_MOTHER(rsqrt_mother,rsqrt_kernel_gpu,rsqrt_kernel_gpu);
MD_MOTHER(mul_mother,mul_kernel_gpu,mul_kernel_gpu);
#ifdef MD_NACL
MD_MOTHER(nacl_mother,nacl_kernel_gpu,nacl_kernel_gpu);
#endif
MD_MOTHER(gravfb_mother,gravfb_kernel_gpu,gravfb_kernel_gpu_rounduploop);
MD_MOTHER(gravpotfb_mother,gravpotfb_kernel_gpu,gravpotfb_kernel_gpu_rounduploop);
MD_MOTHER(ljfb_mother,ljfb_kernel_gpu,ljfb_kernel_gpu_rounduploop);
MD_MOTHER(ljpotfb_mother,ljpotfb_kernel_gpu,ljpotfb_kernel_gpu_rounduploop);
MD_MOTHER(realfb_mother,realfb_kernel_gpu,realfb_kernel_gpu_rounduploop);
MD_MOTHER(realpotfb_mother,realpotfb_kernel_gpu,realpotfb_kernel_gpu_rounduploop);

#ifdef MD_CELLINDEX
#if 1
//MD_MOTHER(ljfbci_mother,ljfb_kernel_gpu,ljfb_kernel_gpu_rounduploop);
//MD_MOTHER(ljpotfbci_mother,ljpotfb_kernel_gpu,ljpotfb_kernel_gpu_rounduploop);
//MD_MOTHER(ljfbci_mother,ljfbci_kernel_gpu,ljfbci_kernel_gpu);
//MD_MOTHER(ljpotfbci_mother,ljpotfbci_kernel_gpu,ljpotfbci_kernel_gpu);
MD_MOTHER_CI(ljfbci_mother,ljfbci_kernel_gpu);
MD_MOTHER_CI(ljpotfbci_mother,ljpotfbci_kernel_gpu);
#endif
MD_MOTHER_CI(realci_mother,realci_kernel_gpu);
MD_MOTHER_CI(realpotci_mother,realpotci_kernel_gpu);
MD_MOTHER_CI(realfbci_mother,realfbci_kernel_gpu);
MD_MOTHER_CI(realpotfbci_mother,realpotfbci_kernel_gpu);
//MD_MOTHER_CI(realljfbci_mother,realljfbci_kernel_gpu);
#if defined(MD_QAUNION_ATYPEBIT) && 1
MD_MOTHER_CI(realljfbci_mother,realljfbci_kernel_gpu);
MD_MOTHER_CI(realljfbci2_mother,realljfbci_kernel_gpu);
MD_MOTHER_CI(realljpotfbci_nojdiv_mother,realljpotfbci_nojdiv_kernel_gpu);
MD_MOTHER_CI(realljpotfbci_mother,realljpotfbci_kernel_gpu);
#endif
#endif

extern "C"
void debug_mother(void *unit_org)
{
#if 0
  VG_UNIT *unit=(VG_UNIT *)unit_org;
  int nf=unit->nf;
  VG_FVEC *fvec=(VG_FVEC *)unit->fvectors;
  unsigned int mem_size_f = sizeof(VG_FVEC)*NF_ROUNDUP(nf);
  VG_FVEC *d_fvec;

  CUDA_SAFE_CALL(cudaMalloc((void**)&d_fvec,mem_size_f));
  CUDA_SAFE_CALL(cudaFree(d_fvec));
#endif

#if 0 // print the size of struct
  printf("in debug_mother: size of VG_UNIT=%d, M3_CELL=%d, VG_JVEC=%d VG_IVEC=%d, VG_MATRIX=%d, VG_SCALER=%d, VG_FVEC=%d, VG_PSCALER=%d, VG_UNION_FI=%d\n",
	 sizeof(VG_UNIT),sizeof(M3_CELL),sizeof(VG_JVEC),sizeof(VG_IVEC),
	 sizeof(VG_MATRIX),sizeof(VG_SCALER),sizeof(VG_FVEC),
	 sizeof(VG_PSCALER),sizeof(VG_UNION_FI));
#endif  

}

#else
extern "C"
void lj_mother(void *unit_org)
{
  VG_UNIT *unit=(VG_UNIT *)unit_org;
  int ni=unit->ni;
  int nj=unit->nj;
  int natj=unit->natj;
  VG_IVEC *ivec=(VG_IVEC *)unit->ivectors;
  VG_JVEC *jvec=(VG_JVEC *)unit->jvectors; 
  VG_MATRIX *matrix=(VG_MATRIX *)unit->matrices;
  VG_SCALER *scalers=(VG_SCALER *)unit->scalers;
  VG_FVEC *fvec=(VG_FVEC *)unit->fvectors;
  VG_PSCALER *pscal=(VG_PSCALER *)unit->pscalers;
  unsigned int mem_size_i = sizeof(VG_IVEC)*NI_ROUNDUP(ni);
  unsigned int mem_size_j = sizeof(VG_JVEC)*NJ_ROUNDUP(nj);
  unsigned int mem_size_matrix = sizeof(VG_MATRIX)*VG_MINIMUM_ATYPE_BLOCK;
  unsigned int mem_size_scalers = sizeof(VG_SCALER);
  unsigned int mem_size_f = sizeof(VG_FVEC)*NF_ROUNDUP(ni);
  unsigned int mem_size_p = sizeof(VG_PSCALER);
  VG_IVEC *d_ivec;
  VG_JVEC *d_jvec;
#ifndef MD_USE_CONSTANT
  VG_MATRIX *d_matrix;
  VG_SCALER *d_scalers;
#endif
  VG_FVEC *d_fvec;
  VG_PSCALER *d_pscal;

#if VG_MINIMUM_PARTICLE_BLOCK_I!=VG_MINIMUM_PARTICLE_BLOCK_J
  fprintf(stderr,"** VG_MINIMUM_PARTICLE_BLOCK_I and VG_MINIMUM_PARTICLE_BLOCK_J must be the same **\n");
  exit(1);
#endif
  MY_CUT_DEVICE_INIT(Argc,Argv);
    
  CUDA_SAFE_CALL(cudaMalloc((void**)&d_ivec,mem_size_i));
  CUDA_SAFE_CALL(cudaMalloc((void**)&d_jvec,mem_size_j));
#ifndef MD_USE_CONSTANT
  CUDA_SAFE_CALL(cudaMalloc((void**)&d_matrix,mem_size_matrix));
  CUDA_SAFE_CALL(cudaMalloc((void**)&d_scalers,mem_size_scalers));
#endif
  CUDA_SAFE_CALL(cudaMalloc((void**)&d_fvec,mem_size_f));
  CUDA_SAFE_CALL(cudaMalloc((void**)&d_pscal,mem_size_p));

  CUDA_SAFE_CALL(cudaMemcpy(d_ivec,ivec,mem_size_i,cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(d_jvec,jvec,mem_size_j,cudaMemcpyHostToDevice));
#ifdef MD_USE_CONSTANT
  CUDA_SAFE_CALL(cudaMemcpyToSymbol(d_matrix,matrix,mem_size_matrix));
  CUDA_SAFE_CALL(cudaMemcpyToSymbol(d_scalers,scalers,mem_size_scalers));
#else
  CUDA_SAFE_CALL(cudaMemcpy(d_matrix,matrix,mem_size_matrix,cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(d_scalers,scalers,mem_size_scalers,cudaMemcpyHostToDevice));
#endif
  CUDA_SAFE_CALL(cudaMemcpy(d_fvec,fvec,mem_size_f,cudaMemcpyHostToDevice));
  CUDA_SAFE_CALL(cudaMemcpy(d_pscal,pscal,mem_size_p,cudaMemcpyHostToDevice));
  
  dim3 threads(VG_MINIMUM_PARTICLE_BLOCK_I);
  dim3 grid((ni+VG_MINIMUM_PARTICLE_BLOCK_I-1)/VG_MINIMUM_PARTICLE_BLOCK_I);
  //  printf("ni=%d nj=%d natj=%d\n",ni,nj,natj);
  lj_kernel_gpu<<< grid, threads >>>(ni,nj,natj,d_ivec,d_jvec,
#ifndef MD_USE_CONSTANT
				     d_matrix,d_scalers,
#endif
				     d_fvec,d_pscal);
  CUT_CHECK_ERROR("Kernel execution failed");
  
  CUDA_SAFE_CALL(cudaMemcpy(fvec,d_fvec,mem_size_f,cudaMemcpyDeviceToHost));
  CUDA_SAFE_CALL(cudaMemcpy(pscal,d_pscal,mem_size_p,cudaMemcpyDeviceToHost));

  //  for(int i=0;i<10;i++) printf("fvec[%d].fi=%e %e %e\n",i,fvec[i].fi[0],fvec[i].fi[1],fvec[i].fi[2]);
			   
  
  CUDA_SAFE_CALL(cudaFree(d_ivec));
  CUDA_SAFE_CALL(cudaFree(d_jvec));
#ifndef MD_USE_CONSTANT
  CUDA_SAFE_CALL(cudaFree(d_matrix));
  CUDA_SAFE_CALL(cudaFree(d_scalers));
#endif
  CUDA_SAFE_CALL(cudaFree(d_fvec));
  CUDA_SAFE_CALL(cudaFree(d_pscal));
}
#endif
