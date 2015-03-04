#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

#include <cutil.h>

#include <gpuewreforce_kernel.cu>
////////////////////////////////////////////////////////////////////////////////
extern "C"
void gpuewreforce_ (double*, int, double*, double, int, double, int, int, double*);
void gpuewreforce__ (double*, int, double*, double, int, double, int, int, double*);

extern "C"
//void ewreforceGold( float*, const double*, const double*, unsigned int, double);
//void ewreforceGold_d( double*, const double*, const double*, unsigned int, double);
void ewreforceGold_f( float*, const float*, unsigned int, float, float);
////////////////////////////////////////////////////////////////////////////////

//extern "C"
void
gpuewreforce_ (double* x, int n, double* q, double rscale, int tblno, double xmax, int periodicflag, int natchangeflag, double* force)
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
   ewreforce_kernel<<< grid, threads >>>(d_C, d_A, n, xmax_float,alpha_float);
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
     //printf("%16.6f %16.6f %d \n",f_float[i*3],reference[i*3],i);
     //printf("%16.6f %16.6f %d \n",f_float[i*3+1],reference[i*3+1],i);
     //printf("%16.6f %16.6f %d \n",f_float[i*3+2],reference[i*3+2],i);
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
gpuewreforce__ (double* x, int *n, double* q, double *rscale, int *tblno, double *xmax, int *periodicflag, int *natchangeflag, double* force)
{
  gpuewreforce_ (x,*n,q,*rscale,*tblno,*xmax,*periodicflag,*natchangeflag,force);
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
 /*
void
ewrepotGold(float* C, const double* A, const double* B, unsigned int num_a, double xmax)
{

  double sum;
  double dn2;
  double dn6;
  double dx;
  double dy;
  double dz;
  double l2 = xmax * 0.5;
  double exclude_radius2 = 0.1;

  for (unsigned int i = 0; i < num_a; ++i){
    sum = 0;
    for (unsigned int j = 0; j < num_a; ++j) {

      dx = (A[i*3]   - A[j*3]  );
      dy = (A[i*3+1] - A[j*3+1]);
      dz = (A[i*3+2] - A[j*3+2]);

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

      dn2 = (dx * dx + dy * dy + dz * dz) * 1;
      if (dn2 > exclude_radius2){
	dn6 = 1.e0 / (dn2 * dn2 * dn2);
	sum += 4 * dn6 * (dn6 - 1.e0);
      }
    }
    C[i*3] = (float)sum;
  }
}
void
ewrepotGold_d(double* C, const double* A, const double* B, unsigned int num_a, double xmax)
{

  double sum;
  double dn2;
  double dn6;
  double dx;
  double dy;
  double dz;
  double l2 = xmax * 0.5;
  double exclude_radius2 = 0.01e0;
  double cutoff_radius2 = 9.e0;

  for (unsigned int i = 0; i < num_a; ++i){
    sum = 0.e0;
    for (unsigned int j = 0; j < num_a; ++j) {

      dx = (A[i*3]   - A[j*3]  );
      dy = (A[i*3+1] - A[j*3+1]);
      dz = (A[i*3+2] - A[j*3+2]);

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
      if ((i != j) && dn2 < cutoff_radius2){
	dn6 = 1.e0 / dn2 / dn2 / dn2 / dn2 / dn2 / dn2
	  - 1.e0 / dn2 / dn2 / dn2;
	sum += dn6 * 4.e0;
      }
      
      dn2 = (dx * dx + dy * dy + dz * dz) * 1;
      if (dn2 > exclude_radius2 && dn2 < cutoff_radius2){
	dn6 = 1.e0 / (dn2 * dn2 * dn2);
	sum += 4 * dn6 * (dn6 - 1.e0);
	}
    }
    C[i*3] = sum;
  }
}
*/
void
ewreforceGold_f(float* C, const float* A, unsigned int num_a, float xmax, float alpha)
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
