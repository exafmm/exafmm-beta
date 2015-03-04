#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

#include <cutil.h>

#include <gpucoulombpot_kernel.cu>
////////////////////////////////////////////////////////////////////////////////
extern "C"
void gpucoulombpot_ (double*, int, double*, double, int, double, int, int, double*);
void gpucoulombpot__ (double*, int, double*, double, int, double, int, int, double*);

extern "C"
void coulombpotGold_f( float*, const float*, unsigned int, float, float);
////////////////////////////////////////////////////////////////////////////////

//extern "C"
void
gpucoulombpot_ (double* x, int n, double* q, double rscale, int tblno, double xmax, int periodicflag, int natchangeflag, double* force)
{
#ifndef CUDA_SDK_2
   CUT_DEVICE_INIT();
#endif
//   CUT_CHECK_DEVICE();

   unsigned int size_A = ((n+THD-1)/THD*THD) * 4;
   unsigned int mem_size_A = sizeof(float) * size_A;
   float* x_float = (float*) malloc(mem_size_A);

   //double stime,ltime;

   for (int i = 0; i < size_A/4; i++){
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

   dim3 threads(THD);
   dim3 grid((n+THD-1) / THD);
   coulombpot_kernel<<< grid, threads >>>(d_C, d_A, n, xmax_float,alpha_float);
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
   
   //get_cputime(&ltime,&stime);
   //printf("HOST Processing time: %10.3f (sec)\n", stime);
   
   //double sum_gpu = 0.e0;
   //double sum_host = 0.e0;
   for (int i = 0; i < n; ++i){
     //sum_gpu += (double)h_C[i*3];
     //sum_host += reference[i*3+1];
     //printf("%16.6f %16.6f %d \n",f_float[i*3],reference[i*3],i);
     force[i*3]   += (double)f_float[i*3];
     //force[i*3+1] += (double)f_float[i*3+1];
     //force[i*3+2] += (double)f_float[i*3+2];
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

/*
extern "C"
void
gpuvdwpot_ (double* x, int *n, int* atype, int *nat, double* gscale, double* rscale, int *tblno, double *xmax, int *periodicflag, int *natchangeflag, double* force)
{
  gpuvdwpot_ (x,*n,atype,*nat,gscale,rscale,*tblno,*xmax,*periodicflag,*natchangeflag,force);
}
*/

extern "C"
void
gpucoulombpot__ (double* x, int *n, double* q, double *rscale, int *tblno, double *xmax, int *periodicflag, int *natchangeflag, double* force)
{
  gpucoulombpot_ (x,*n,q,*rscale,*tblno,*xmax,*periodicflag,*natchangeflag,force);
  }

/*
extern "C"
void printDiff(float *data1, float *data2, int width, int height)
{
  int i,j,k;
  int error_count=0;
  for (j=0; j<height; j++) {
    for (i=0; i<width; i++) {
      k = j*width+i;
      if (data1[k] != data2[k]) {
	printf("diff(%d,%d) CPU=%4.4f, GPU=%4.4f\n", i,j, data1[k], data2[k]);
	error_count++;
      }
    }
  }
  printf("\nTotal Errors = %d\n", error_count);
}
*/
void
coulombpotGold_f(float* C, const float* A, unsigned int num_a, float xmax, float alpha)
{

  float sum;
  float dn2;
  float dx;
  float dy;
  float dz;
  float l2 = xmax * 0.5;
  float sqdn;
  //float exclude_radius2 = 0.01e0;
  //float cutoff_radius2 = 9.e0;

  for (unsigned int i = 0; i < num_a; ++i){
    sum = 0.e0;
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

      dn2 = dx * dx + dy * dy + dz * dz;
      //if ((i != j) && dn2 < cutoff_radius2){
      if ((i != j)){
	sqdn = sqrt(dn2) * alpha;
	sum += erfc(sqdn) / sqdn * A[j*4+3];
      }
    }
    C[i*3] = sum * A[i*4+3] * alpha;
  }
}
