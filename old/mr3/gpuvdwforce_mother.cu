#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

#include <cutil.h>

#include <gpuvdwforce_kernel.cu>
////////////////////////////////////////////////////////////////////////////////
extern "C"
void gpuvdwforce_ (double*, int, int*, int, double*, double*, int, double, int, int, double*,
		   double, double);
void gpuvdwforce__ (double*, int, int*, int, double*, double*, int, double, int, int, double*,
		    double *, double *);
//void printDiff(float*, float*, int, int);

extern "C"
//void computeGold2( float*, const double*, const double*, unsigned int, double);
//void computeGold2_d(double*, const double*, const double*, unsigned int, double);
//void computeGold2_d2(double*, const double*, const double*, unsigned int, double);
//void computeGold2_f2(float*, const float*, const float*, unsigned int, float);
////////////////////////////////////////////////////////////////////////////////
//! One dimensional dimensionless Lennard-Jones Simulation
////////////////////////////////////////////////////////////////////////////////
/*
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
*/

extern "C"
void
gpuvdwforce_ (double* x, int n, int* atype, int nat_org, double* gscale, double* rscale, int tblno, double xmax, int periodicflag, int natchangeflag, double* force,
	      double r2mind, double r2maxd)
{
#ifndef CUDA_SDK_2
   CUT_DEVICE_INIT();
#endif
//   CUT_CHECK_DEVICE();

   unsigned int size_A = ((n+THD-1)/THD*THD) * 3;
   unsigned int mem_size_A = sizeof(float) * size_A;
   float* x_float = (float*) malloc(mem_size_A);
   float* f_float = (float*) malloc(mem_size_A);
   unsigned int size_B = THD;
   unsigned int mem_size_B = sizeof(float) * size_B;
   float* gr_float = (float*) malloc(mem_size_B);

   //double stime,ltime;
   float r2min=(float)r2mind;
   float r2max=(float)r2maxd;
   int nat=nat_org+1;

   //   printf("n=%d n*3=%d size_A=%d\n",n,n*3,size_A);
   for (int i = 0; i < size_A; i++){
     if(i<n*3) x_float[i] = (float)x[i];
     else      x_float[i]=0.0;
   }
   for(int i=0;i<nat;i++){
     for(int j=0;j<nat;j++){
       if(i<nat_org && j<nat_org){
	 gr_float[(i*nat+j)*2]   = (float)gscale[i*nat_org+j];
	 gr_float[(i*nat+j)*2+1] = (float)rscale[i*nat_org+j];
       }
       else{
	 gr_float[(i*nat+j)*2]   = 0.0f;
	 gr_float[(i*nat+j)*2+1] = 1.0f;
       }
     }
   }
   
   float* d_B;
   CUDA_SAFE_CALL(cudaMalloc((void**) &d_B, mem_size_B));
   CUDA_SAFE_CALL(cudaMemcpy(d_B, gr_float, mem_size_B,cudaMemcpyHostToDevice) );
   float* d_A;
   CUDA_SAFE_CALL(cudaMalloc((void**) &d_A, mem_size_A));
   CUDA_SAFE_CALL(cudaMemcpy(d_A, x_float, mem_size_A,cudaMemcpyHostToDevice) );
   unsigned int size_D = (n+THD-1)/THD*THD;
   unsigned int mem_size_D = sizeof(int) * size_D;
   int* d_D;
   int *atype2 = (int *)malloc(mem_size_D);

   for(int i=0;i<size_D;i++){
     if(i<n) atype2[i]=atype[i];
     else    atype2[i]=nat_org;
   }
   CUDA_SAFE_CALL(cudaMalloc((void**) &d_D, mem_size_D));
   CUDA_SAFE_CALL(cudaMemcpy(d_D, atype2, mem_size_D,cudaMemcpyHostToDevice) );

   unsigned int size_C = size_A;
   unsigned int mem_size_C = sizeof(float) * size_C;
   float* d_C;
   CUDA_SAFE_CALL(cudaMalloc((void**) &d_C, mem_size_C));

   //get_cputime(&ltime,&stime);

   float xmax_float = (float)xmax;

   dim3 threads(THD);
   dim3 grid((n+THD-1)/THD);
   gpuvdwforce_kernel<<< grid, threads >>>(d_C, d_A, d_B, d_D, n, xmax_float, nat,
					   r2min,r2max);

   CUT_CHECK_ERROR("Kernel execution failed");

   CUDA_SAFE_CALL(cudaMemcpy(f_float, d_C, mem_size_C,cudaMemcpyDeviceToHost) );
   //get_cputime(&ltime,&stime);
   //printf("GPU  Processing time: %10.3f (sec)\n", stime);
   /*
   float* reference = (float*) malloc(mem_size_A);
   for (int i = 0; i < n; ++i){
     reference[i*3] = 0.e0;
     reference[i*3+1] = 0.e0;
     reference[i*3+2] = 0.e0;
   }
   */
   //get_cputime(&ltime,&stime);
   //computeGold2_f2(reference, x_float, x_float, n, xmax_float);
   
   //computeGold(reference, x, x, n, xmax);
   /*
   unsigned int size_E = n * 3;
   unsigned int mem_size_E = size_E * sizeof(double);
   double* force_double = (double*) malloc(mem_size_E);

   for (int i = 0; i < n; ++i){
     force_double[i*3] = 0.e0;
     force_double[i*3+1] = 0.e0;
     force_double[i*3+2] = 0.e0;
     }*/
   //computeGold2_d2(force_double, x, x, n, xmax);
   
   //get_cputime(&ltime,&stime);
   //printf("HOST Processing time: %10.3f (sec)\n", stime);

   //double sum_gpu = 0.e0;
   //double sum_host = 0.e0;
   
   for (int i = 0; i < n; ++i){
     //sum_gpu += (double)f_float[i*3+1];
     //sum_host += force_double[i*3+1];
     //sum_host += (double)reference[i*3+1];
     //printf("%16.6f %16.6f %d \n",h_C[i*3],reference[i*3],i);

     //printf("%16.6f %16.6f %d \n",f_float[i*3],f_float[i*3+1],i);
     
     force[i*3]   += (double)f_float[i*3];
     force[i*3+1] += (double)f_float[i*3+1];
     force[i*3+2] += (double)f_float[i*3+2];
     
     /*
     force[i*3]   += force_double[i*3];
     force[i*3+1] += force_double[i*3+1];
     force[i*3+2] += force_double[i*3+2];
     */
   }
   
   //printf("HOST: %20.8f  \n",force_double[5]);
   //printf("HOST: %20.8f  \n",reference[5]);
   //printf("GPU : %20.8f  \n",f_float[5]);
   

   free(x_float);
   free(f_float);
   free(gr_float);
   //free(force_double);
   //free(reference);
   CUDA_SAFE_CALL(cudaFree(d_A));
   CUDA_SAFE_CALL(cudaFree(d_B));
   CUDA_SAFE_CALL(cudaFree(d_C));
   CUDA_SAFE_CALL(cudaFree(d_D));
}

/*
extern "C"
void
gpuvdwforce_ (double* x, int *n, int* atype, int *nat, double* gscale, double* rscale, int *tblno, double *xmax, int *periodicflag, int *natchangeflag, double* force)
{
  gpuvdwforce_ (x,*n,atype,*nat,gscale,rscale,*tblno,*xmax,*periodicflag,*natchangeflag,force);
}
*/

extern "C"
void
gpuvdwforce__ (double* x, int *n, int* atype, int *nat, double* gscale, double* rscale, int *tblno, double *xmax, int *periodicflag, int *natchangeflag, double* force,
	       double *r2mind, double *r2maxd)
{
  gpuvdwforce_ (x,*n,atype,*nat,gscale,rscale,*tblno,*xmax,*periodicflag,*natchangeflag,force,
		*r2mind,*r2maxd);
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
computeGold2(float* C, const double* A, const double* B, unsigned int num_a, double xmax)
{

  double dn2, adn2, tmp2;
  double dx, dy, dz;
  double f, fx, fy, fz;
  double l2 = 0.5e0 * xmax;
  double exclude_radius2 = 0.1;

  for (unsigned int i = 0; i < num_a; ++i){
    fx = 0.e0;
    fy = 0.e0;
    fz = 0.e0;
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
	adn2 = 1.e0 / dn2;
	tmp2 = adn2 * adn2 * adn2;
	f = 24.e0 * adn2 * tmp2 * (2.e0 * tmp2 - 1.e0);
	fx += f * dx;
	fy += f * dy;
	fz += f * dz;
	
      }
    }
    C[i*3]   = (float)fx;
    C[i*3+1] = (float)fy;
    C[i*3+2] = (float)fz;
  }
}

void
computeGold2_d(double* C, const double* A, const double* B, unsigned int num_a, double xmax)
{

  double dn2, tmp2;
  double dx, dy, dz;
  double f, fx, fy, fz;
  double l2 = 0.5e0 * xmax;
  //double exclude_radius2 = 0.01e0;
  double cutoff_radius2 = 9.e0;

  for (unsigned int i = 0; i < num_a; ++i){
    fx = 0.e0;
    fy = 0.e0;
    fz = 0.e0;
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
	tmp2 = 2.e0 / dn2 / dn2 / dn2 / dn2 / dn2 / dn2 / dn2
	  - 1.e0 / dn2 / dn2 / dn2 / dn2;
	f = tmp2 * 24.e0;
	fx += f * dx;
	fy += f * dy;
	fz += f * dz;

	
      }
      /*
      dn2 = (dx * dx + dy * dy + dz * dz) * 1;
      if (dn2 > exclude_radius2 && dn2 < cutoff_radius2){
	adn2 = 1.e0 / dn2;
	tmp2 = adn2 * adn2 * adn2;
	f = 24.e0 * adn2 * tmp2 * (2.e0 * tmp2 - 1.e0);
	fx += f * dx;
	fy += f * dy;
	fz += f * dz;
	
	}*/
    }
    C[i*3]   = fx;
    C[i*3+1] = fy;
    C[i*3+2] = fz;
	if (fz > 100.e0)
	  printf("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
  }
}
void
computeGold2_d2(double* C, const double* A, const double* B, unsigned int num_a, double xmax)
{

  double dn2, tmp2;
  double dx, dy, dz;
  double f;
  double l2 = 0.5e0 * xmax;
  //double exclude_radius2 = 0.01e0;
  double cutoff_radius2 = 9.e0;

  for (unsigned int i = 0; i < num_a-1; ++i){
    for (unsigned int j = i+1; j < num_a; ++j) {

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
      if (dn2 < cutoff_radius2){
	tmp2 = 2.e0 / dn2 / dn2 / dn2 / dn2 / dn2 / dn2 / dn2
	  - 1.e0 / dn2 / dn2 / dn2 / dn2;
	f = tmp2 * 24.e0;
	C[i*3]   += f * dx;
	C[i*3+1] += f * dy;
	C[i*3+2] += f * dz;
	C[j*3]   -= f * dx;
	C[j*3+1] -= f * dy;
	C[j*3+2] -= f * dz;
	
      }
      /*
      dn2 = (dx * dx + dy * dy + dz * dz) * 1;
      if (dn2 > exclude_radius2 && dn2 < cutoff_radius2){
	adn2 = 1.e0 / dn2;
	tmp2 = adn2 * adn2 * adn2;
	f = 24.e0 * adn2 * tmp2 * (2.e0 * tmp2 - 1.e0);
	fx += f * dx;
	fy += f * dy;
	fz += f * dz;
	
	}*/
    }
  }
}
void
computeGold2_f2(float* C, const float* A, const float* B, unsigned int num_a, float xmax)
{

  float dn2, tmp2;
  float dx, dy, dz;
  float f;
  float l2 = 0.5e0 * xmax;
  //double exclude_radius2 = 0.01e0;
  float cutoff_radius2 = 9.e0;

  for (unsigned int i = 0; i < num_a-1; ++i){
    for (unsigned int j = i+1; j < num_a; ++j) {

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
      if (dn2 < cutoff_radius2){
	tmp2 = 2.e0 / dn2 / dn2 / dn2 / dn2 / dn2 / dn2 / dn2
	  - 1.e0 / dn2 / dn2 / dn2 / dn2;
	f = tmp2 * 24.e0;
	C[i*3]   += f * dx;
	C[i*3+1] += f * dy;
	C[i*3+2] += f * dz;
	C[j*3]   -= f * dx;
	C[j*3+1] -= f * dy;
	C[j*3+2] -= f * dz;
	
      }
      /*
      dn2 = (dx * dx + dy * dy + dz * dz) * 1;
      if (dn2 > exclude_radius2 && dn2 < cutoff_radius2){
	adn2 = 1.e0 / dn2;
	tmp2 = adn2 * adn2 * adn2;
	f = 24.e0 * adn2 * tmp2 * (2.e0 * tmp2 - 1.e0);
	fx += f * dx;
	fy += f * dy;
	fz += f * dz;
	
	}*/
    }
  }
}
