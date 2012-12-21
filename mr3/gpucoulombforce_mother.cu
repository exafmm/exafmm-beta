#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

#include <cutil.h>

  //#define THREAD_CONTIGUOUS  // not fast
//#define HALF_THREADS_PER_BLOCK  // not fast

//extern int Deviceid;
static int Deviceid=0;

#include <gpucoulombforce_kernel.cu>
////////////////////////////////////////////////////////////////////////////////
extern "C"
void gpucoulombforce_ (double*, int, double*, double, int, double, int, int, double*, double );
void gpucoulombforce__ (double*, int, double*, double, int, double, int, int, double*, double *);
extern "C"
void gpucoulombforce_ij_(double* xi, int ni, double *qj, double *xj, int nj,
			 double, int, double, int, int, double*, double );

extern "C"
void coulombforceGold_f( float*, const float*, unsigned int, float, float);
////////////////////////////////////////////////////////////////////////////////

//extern "C"
void
gpucoulombforce_ (double* x, int n, double* q, double rscale, int tblno, 
		  double xmax, int periodicflag, int natchangeflag, 
		  double* force, double eps2d)
{
#ifndef CUDA_SDK_2
   CUT_DEVICE_INIT();
#endif
//   CUT_CHECK_DEVICE();

   unsigned int size_A = ((n+THD-1)/THD*THD) * 4;
   unsigned int mem_size_A = sizeof(float) * size_A;
   float* x_float = (float*) malloc(mem_size_A);
   float eps2=(float)eps2d;

   //double stime,ltime;

   for (int i = 0; i < size_A/4; i++){
#ifdef THREAD_CONTIGUOUS
     int offset,ii,i2;
     ii=i/THD;
     offset=i % THD;
     i2=ii*4*THD+offset;
     if(i<n){
       x_float[i2]       = (float)x[i*3];
       x_float[i2+THD]   = (float)x[i*3+1];
       x_float[i2+THD*2] = (float)x[i*3+2];
       x_float[i2+THD*3] = (float)q[i];
     }
     else{
       x_float[i2]       = 0.0f;
       x_float[i2+THD]   = 0.0f;
       x_float[i2+THD*2] = 0.0f;
       x_float[i2+THD*3] = 0.0f;
     }
#else
     if(i<n){
       x_float[i*4]   = (float)x[i*3];
       x_float[i*4+1] = (float)x[i*3+1];
       x_float[i*4+2] = (float)x[i*3+2];
       x_float[i*4+3] = (float)q[i];
     }
     else{
       x_float[i*4]   = 0.0f;
       x_float[i*4+1] = 0.0f;
       x_float[i*4+2] = 0.0f;
       x_float[i*4+3] = 0.0f;
     }
#endif
   }
   
   float xmax_float = (float)xmax;
   float alpha_float = (float)rscale;

   float* d_A;
   CUDA_SAFE_CALL(cudaMalloc((void**) &d_A, mem_size_A));
   CUDA_SAFE_CALL(cudaMemcpy(d_A, x_float, mem_size_A,cudaMemcpyHostToDevice) );

   unsigned int size_C = ((n+THD-1)/THD*THD) * 3;
   unsigned int mem_size_C = sizeof(float) * size_C;
   float* d_C;
   CUDA_SAFE_CALL(cudaMalloc((void**) &d_C, mem_size_C));
   float* f_float = (float*) malloc(mem_size_C);

   //get_cputime(&ltime,&stime);

#ifdef HALF_THREADS_PER_BLOCK   
   dim3 threads(THD/2);
#else
   dim3 threads(THD);
#endif
   dim3 grid((n+THD-1) / THD);
   coulombforce_kernel<<< grid, threads >>>(d_C, d_A, n, xmax_float,
					    alpha_float,eps2);
   CUT_CHECK_ERROR("Kernel execution failed");
   CUDA_SAFE_CALL(cudaMemcpy(f_float, d_C, mem_size_C,cudaMemcpyDeviceToHost) );

   //get_cputime(&ltime,&stime);
   //printf("GPU  Processing time: %10.3f (sec)\n", stime);

   //float* reference = (float*) malloc(mem_size_C);
   /*
   for (int i = 0; i < n; ++i){
     reference[i*3]   = 0.e0;
     reference[i*3+1] = 0.e0;
     reference[i*3+2] = 0.e0;
     }*/

   //get_cputime(&ltime,&stime);
   
   //ewreforceGold_f(reference, x_float, n, xmax_float, alpha_float);

   /*
   unsigned int size_E = n * 3;
   unsigned int mem_size_E = size_E * sizeof(double);
   double* force_double = (double*) malloc(mem_size_E);
   for (int i = 0; i < n; ++i){
     force_double[i*3] = 0.e0;
     force_double[i*3+1] = 0.e0;
     force_double[i*3+2] =0.e0;
   }
   */
   //ewrepotGold_d(force_double, x, x, n, xmax);
   
   //get_cputime(&ltime,&stime);
   //printf("HOST Processing time: %10.3f (sec)\n", stime);
   
   //double sum_gpu = 0.e0;
   //double sum_host = 0.e0;
   for (int i = 0; i < n; ++i){
     //sum_gpu += (double)h_C[i*3];
     //sum_host += reference[i*3+1];
     //     printf("%16.6f %d \n",f_float[i*3],i);
     //     printf("%16.6f %d \n",f_float[i*3+1],i);
     //     printf("%16.6f %d \n",f_float[i*3+2],i);
     force[i*3]   += (double)f_float[i*3];
     force[i*3+1] += (double)f_float[i*3+1];
     force[i*3+2] += (double)f_float[i*3+2];
     /*
     force[i*3]   += force_double[i*3];
     force[i*3+1] += force_double[i*3+1];
     force[i*3+2] += force_double[i*3+2];
     */
   }
   //printf("GPU : %20.8f  \n",sum_gpu);
   
   //printf("GPU : %20.8f  \n",f_float[6]);
   //printf("HOST: %20.8f  \n",reference[6]);

   free(x_float);
   free(f_float);
   //free(reference);
   CUDA_SAFE_CALL(cudaFree(d_A));
   //CUDA_SAFE_CALL(cudaFree(d_B));
   CUDA_SAFE_CALL(cudaFree(d_C));
}

extern "C"
void
gpucoulombforce__ (double* x, int *n, double* q, double *rscale, int *tblno, 
		   double *xmax, int *periodicflag, int *natchangeflag, 
		   double* force, double *eps2d)
{
  gpucoulombforce_ (x,*n,q,*rscale,*tblno,*xmax,
		    *periodicflag,*natchangeflag,force,*eps2d);
  }

extern "C" 
void
gpucoulombforce_ij_(double* xi, int ni, double* qj, double *xj, int nj,
		    double rscale, int tblno, 
		    double xmax, int periodicflag, int natchangeflag, 
		    double* force, double eps2d)
{
   int s_gpuCount=0;
#ifndef CUDA_SDK_2
   CUT_DEVICE_INIT();
#endif
   CUDA_SAFE_CALL(cudaGetDeviceCount(&s_gpuCount));
   if(Deviceid<s_gpuCount){
     CUDA_SAFE_CALL(cudaSetDevice(Deviceid));
   }
   else{
     fprintf(stderr,"** error : no such device id=%d **\n",Deviceid);
     exit(1);
   }

   unsigned int size_A = ((nj+THD-1)/THD*THD) * 4;
   unsigned int mem_size_A = sizeof(float) * size_A;
   float* xj_float;
   unsigned int size_B = ((ni+THD-1)/THD*THD) * 4;
   unsigned int mem_size_B = sizeof(float) * size_B;
   float* xi_float;
   float eps2=(float)eps2d;
   unsigned int size_C = ((ni+THD-1)/THD*THD) * 3;
   unsigned int mem_size_C = sizeof(float) * size_C;
   float* f_float;

#ifdef NO_SHARED_FOR_I
#if VMP!=1 || !defined(UNROLL)
   fprintf(stderr,"** error : VMP must be 1, and UNROLL must be defined when NO_SHARED_FOR_I is defined **\n");
   exit(1);
#endif
#endif

   xj_float = (float*) malloc(mem_size_A);
   xi_float = (float*) malloc(mem_size_B);
   f_float = (float*) malloc(mem_size_C);
   
   for (int i = 0; i < size_A/4; i++){
     if(i<nj){
       xj_float[i*4]   = (float)xj[i*3];
       xj_float[i*4+1] = (float)xj[i*3+1];
       xj_float[i*4+2] = (float)xj[i*3+2];
       xj_float[i*4+3] = (float)qj[i];
     }
     else{
       xj_float[i*4]   = 0.0f;
       xj_float[i*4+1] = 0.0f;
       xj_float[i*4+2] = 0.0f;
       xj_float[i*4+3] = 0.0f;
     }
   }
   for (int i = 0; i < size_B/4; i++){
     if(i<ni){
       xi_float[i*4]   = (float)xi[i*3];
       xi_float[i*4+1] = (float)xi[i*3+1];
       xi_float[i*4+2] = (float)xi[i*3+2];
       xi_float[i*4+3] = 0.0f;
     }
     else{
       xi_float[i*4]   = 0.0f;
       xi_float[i*4+1] = 0.0f;
       xi_float[i*4+2] = 0.0f;
       xi_float[i*4+3] = 0.0f;
     }
   }
   
   float xmax_float = (float)xmax;
   float alpha_float = (float)rscale;
   static float* d_A=NULL,*d_B=NULL;
   if(d_A==NULL) CUDA_SAFE_CALL(cudaMalloc((void**) &d_A, mem_size_A));
   if(d_B==NULL) CUDA_SAFE_CALL(cudaMalloc((void**) &d_B, mem_size_B));
   CUDA_SAFE_CALL(cudaMemcpy(d_A, xj_float, mem_size_A,cudaMemcpyHostToDevice) );
   CUDA_SAFE_CALL(cudaMemcpy(d_B, xi_float, mem_size_B,cudaMemcpyHostToDevice) );

   static float* d_C=NULL;
   if(d_C==NULL) CUDA_SAFE_CALL(cudaMalloc((void**) &d_C, mem_size_C));

   dim3 threads(THD);
   dim3 grid((ni+THD-1) / THD);
   coulombforce_ij_kernel<<< grid, threads >>>(d_C, d_B, d_A, nj, xmax_float,
					    alpha_float,eps2);
   CUT_CHECK_ERROR("Kernel execution failed");
   CUDA_SAFE_CALL(cudaMemcpy(f_float, d_C, mem_size_C,cudaMemcpyDeviceToHost) );
   if(eps2d!=0.0){
     for (int i = 0; i < ni; ++i){
       force[i*3]   = -(double)f_float[i*3];
       force[i*3+1] = -(double)f_float[i*3+1];
       force[i*3+2] = -(double)f_float[i*3+2];
     }
   }
   else{
     for (int i = 0; i < ni; ++i){
       force[i*3]   += (double)f_float[i*3];
       force[i*3+1] += (double)f_float[i*3+1];
       force[i*3+2] += (double)f_float[i*3+2];
     }
   }

   free(xi_float);
   free(xj_float);
   free(f_float);
}


void
coulombforceGold_f(float* C, const float* A, unsigned int num_a, float xmax, float alpha)
{

  float alpha2 = alpha  * alpha;
  float alpha3 = alpha2 * alpha;
  float dx;
  float dy;
  float dz;
  float l2 = xmax * 0.5;
  float sqdn, ar2, tmp;
  //float exclude_radius2 = 0.01e0;
  //float cutoff_radius2 = 9.e0;
  float sum[3];
  float api2 = 2.e0 / sqrt(3.141592653e0);

  for (unsigned int i = 0; i < num_a; ++i){
    sum[0] = 0.e0;
    sum[1] = 0.e0;
    sum[2] = 0.e0;
    for (unsigned int j = 0; j < num_a; ++j) {

      dx = (A[i*4]   - A[j*4]  );
      dy = (A[i*4+1] - A[j*4+1]);
      dz = (A[i*4+2] - A[j*4+2]);

      if (!(dx < l2 && dx > -l2))
	if (dx > l2){
	  dx = dx - xmax;
	}else{
	  dx = dx + xmax;
	}
      if (!(dy < l2 && dy > -l2))
	if (dy > l2){
	  dy = dy - xmax;
	}else{
	  dy = dy + xmax;
	}
      if (!(dz < l2 && dz > -l2))
	if (dz > l2){
	  dz = dz - xmax;
	}else{
	  dz = dz + xmax;
	}

      ar2 = (dx * dx + dy * dy + dz * dz) * alpha2;
      //if ((i != j) && dn2 < cutoff_radius2){
      if ((i != j)){
	sqdn = sqrt(ar2);
	tmp = A[j*4+3] * (api2 * exp(-ar2) + erfc(sqdn) / sqdn) / ar2;
	sum[0] += tmp * dx;
	sum[1] += tmp * dy;
	sum[2] += tmp * dz;
      }
    }
    tmp = A[i*4+3] * alpha3;
    C[i*3]   = sum[0] * tmp;
    C[i*3+1] = sum[1] * tmp;
    C[i*3+2] = sum[2] * tmp;
  }
}
