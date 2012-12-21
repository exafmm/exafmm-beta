#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

#include <cutil.h>

#include <gpuewwvforce_kernel.cu>
////////////////////////////////////////////////////////////////////////////////
extern "C"
void gpuewwvforce_  (double*, int, double*, int, double*, double *);
void gpuewwvforce__ (double*, int, double*, int, double*, double *);

extern "C"
void computeGold_d(double*, const double*, const double*, int , int);
////////////////////////////////////////////////////////////////////////////////
//! 
////////////////////////////////////////////////////////////////////////////////
#if 0
extern "C"
void get_cputime(double *laptime, double *sprittime)
{
  struct timeval tv;
  struct timezone tz;
  double sec,microsec;

  gettimeofday(&tv, &tz);
  sec=tv.tv_sec;
  microsec=tv.tv_usec;

  *sprittime = sec + microsec * 1e-6 - *laptime;
  *laptime = sec + microsec * 1e-6;
}
#endif


extern "C"
void gpuewwvforce_ (double* xq, int num_atm, double* g_vec, int num_k, 
		    double* force, double *tpot)
{
#ifndef CUDA_SDK_2
   CUT_DEVICE_INIT();
#endif
//   CUT_CHECK_DEVICE();

//   double tpot[1];
   unsigned int size_A = ((num_atm+THD-1)/THD*THD) * 4;
   unsigned int mem_size_A = sizeof(float) * size_A;
   float* xq_float = (float*) malloc(mem_size_A);

   unsigned int size_B = ((num_k+THD-1)/THD*THD) * 4;
   unsigned int mem_size_B = sizeof(float) * size_B;
   float* gv_float = (float*) malloc(mem_size_B);

   unsigned int size_C = ((num_k+THD-1)/THD*THD) * 3;
   unsigned int mem_size_C = sizeof(float) * size_C;
   float* pot_float = (float*) malloc(mem_size_C);

   unsigned int size_D = ((num_atm+THD-1)/THD*THD) * 3;
   unsigned int mem_size_D = sizeof(float) * size_D;
   float* f_float = (float*) malloc(mem_size_D);
   //   unsigned int mem_size_E = sizeof(double) * size_D;
   //   double* force_double = (double*) malloc(mem_size_E);

   //double stime,ltime;

   for (int i = 0; i < size_A ; i++){
     if(i<num_atm*4) xq_float[i] = (float)xq[i];
     else            xq_float[i] = 0.0f;
   }
   
   for (int i = 0; i < size_B; i++){
     if(i<num_k*4) gv_float[i] = (float)g_vec[i];
     else          gv_float[i] = 0.0f;
     //printf("%16.6f %d \n",gv_float[i],i);
   }
   /*   for (int i = 0; i < num_k; i++)
     {
       printf("%16.6f %d \n",gv_float[i*4+3],i);
       }*/

   for (int i = 0; i < size_D; i++)
     {
       f_float[i] = 0.e0;
       //       force_double[i] = 0.e0;
     }

   //get_cputime(&ltime,&stime);

   float* d_A;
   CUDA_SAFE_CALL(cudaMalloc((void**) &d_A, mem_size_A));
   CUDA_SAFE_CALL(cudaMemcpy(d_A, xq_float, mem_size_A,cudaMemcpyHostToDevice) );
   float* d_B;
   CUDA_SAFE_CALL(cudaMalloc((void**) &d_B, mem_size_B));
   CUDA_SAFE_CALL(cudaMemcpy(d_B, gv_float, mem_size_B,cudaMemcpyHostToDevice) );
   float* d_C;
   CUDA_SAFE_CALL(cudaMalloc((void**) &d_C, mem_size_C));
   float* d_D;
   CUDA_SAFE_CALL(cudaMalloc((void**) &d_D, mem_size_D));

   dim3 threads(THD);
   dim3 grid((num_k+THD-1) / THD);
   gpuewwvforce_kernel<<< grid, threads >>>(d_C, d_A, d_B, num_atm);
   CUT_CHECK_ERROR("Kernel execution failed");
#if 1
   CUDA_SAFE_CALL(cudaMemcpy(pot_float, d_C, mem_size_C,cudaMemcpyDeviceToHost) );
   *tpot=0.0;
   for(int i=0;i<num_k;i++) *tpot+=pot_float[i*3+2];
   (*tpot)*=0.5;
   //   printf("tpot=%e 0:%f\n",*tpot);
#endif

   dim3 grid2((num_atm+THD-1) / THD);
   gpuewwvforce_kernel2<<< grid2, threads >>>(d_C, d_B, d_A, d_D, num_k);

   CUT_CHECK_ERROR("Kernel execution failed");
   CUDA_SAFE_CALL(cudaMemcpy(f_float, d_D, mem_size_D,cudaMemcpyDeviceToHost) );

   //get_cputime(&ltime,&stime);

   //printf("GPU  Processing time: %10.3f (sec)\n", stime);

   // host cumputation ////////////////////////////////////////

   //   get_cputime(&ltime,&stime);
   
   //   computeGold_d(force_double, g_vec, xq, num_atm, num_k);
   
   //   get_cputime(&ltime,&stime);
   
   //   printf("HOST Processing time: %10.3f (sec)\n", stime);
   // post preccess  ////////////////////////////////////////////
   
   //double err;
   for (int i = 0; i < num_atm; ++i){
     //printf("%16.6f %16.6f %d \n",f_float[i*3],force_double[i*3],i);
     
     force[i*3]   = (double)f_float[i*3];
     force[i*3+1] = (double)f_float[i*3+1];
     force[i*3+2] = (double)f_float[i*3+2];
     /*
     err += fabs((force_double[i*3] - (double)f_float[i*3]) 
		 / force_double[i*3]);
     err += fabs((force_double[i*3+1] - (double)f_float[i*3+1]) 
		 / force_double[i*3+1]);
     err += fabs((force_double[i*3+2] - (double)f_float[i*3+2]) 
     / force_double[i*3+2]);*/
     /*
     force[i*3]   = force_double[i*3];
     force[i*3+1] = force_double[i*3+1];
     force[i*3+2] = force_double[i*3+2];
     */
   }
   //err = err / (3.e0 * (double)num_atm) * 100.e0;
   //printf("err : %20.8f  \n",err);

   //printf("GPU : %20.8f  \n",sum_gpu);
   //printf("GPU : %20.8f \n",force[0]);
   //printf("HOST: %20.8f \n",force_double[0]);

   CUDA_SAFE_CALL(cudaFree(d_A));
   CUDA_SAFE_CALL(cudaFree(d_B));
   CUDA_SAFE_CALL(cudaFree(d_C));
   CUDA_SAFE_CALL(cudaFree(d_D));

   free(xq_float);
   free(gv_float);
   free(f_float);
   free(pot_float);
   //free(force_double);
}

extern "C"
void
gpuewwvforce__ (double* xq, int* num_atm, double* g_vec, int* num_k, 
		double* force, double *tpot)
{
  gpuewwvforce_ (xq,*num_atm,g_vec,*num_k,force,tpot);
}

extern "C"
void
computeGold_d(double* C, const double* A, const double* B, int num_atm, int num_k)
{
  //A : g_vec
  //B : xq

  double kr, qsin, qcos, tmp;

  unsigned int size_sc = num_atm;
  unsigned int mem_size_sc = sizeof(double) * size_sc;
  double* sin_theta = (double*) malloc(mem_size_sc);
  double* cos_theta = (double*) malloc(mem_size_sc);

  for (unsigned int i = 0; i < num_atm; ++i){
    C[i*3]   = 0.e0;
    C[i*3+1] = 0.e0;
    C[i*3+2] = 0.e0;
  }

  for (unsigned int i = 0; i < num_k; ++i){
    qsin = 0.e0;
    qcos = 0.e0;
    for (unsigned int j = 0; j < num_atm; ++j){
      kr = A[i*4]   * B[j*4]
	 + A[i*4+1] * B[j*4+1]
	 + A[i*4+2] * B[j*4+2];
      sin_theta[j] = sin(kr);
      cos_theta[j] = cos(kr);
      qsin += B[j*4+3] * sin_theta[j];
      qcos += B[j*4+3] * cos_theta[j];
    }
    //C[i*3]   = qsin;
    //C[i*3+1] = qcos;
    //C[i*3+2] = A[i*4+3] * ( qsin * qsin + qcos * qcos );
    qsin *= A[i*4+3];
    qcos *= A[i*4+3];
    for (unsigned int j = 0; j < num_atm; ++j){
      tmp = qcos * sin_theta[j] - qsin * cos_theta[j];
      C[j*3]   += tmp * A[i*4];
      C[j*3+1] += tmp * A[i*4+1];
      C[j*3+2] += tmp * A[i*4+2];
    }
  }

  for (unsigned int i = 0; i < num_atm; ++i){
    C[i*3]   *= B[i*4+3];
    C[i*3+1] *= B[i*4+3];
    C[i*3+2] *= B[i*4+3];
  }

  free(sin_theta);
  free(cos_theta);

}
