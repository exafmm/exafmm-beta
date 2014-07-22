/*

A simple program that demonstrate building of an octree using cilk.
This programm test the transpose version

author: Nikos Sismanis
date: Jul 2014

*/

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "sys/time.h"
#include "utils.h"

#define DIM 3
#define LDIM 12

int main(int argc, char** argv){

  struct timeval startwtime, endwtime;

#ifdef LARRAY

  struct Body *data, *Rdata;

#endif

  int N = atoi(argv[1]); // Number of points
  int dist = atoi(argv[2]); // Distribution identifier
  int repeat = atoi(argv[3]);

  printf("Running for %d points...\n", N);
  printf("With distribution: %d\n", dist);

#ifndef LARRAY

  float *X = (float*)malloc(DIM * N * sizeof(float)); // Matrix to hold the initial data
  float *Y = (float*)malloc(DIM * N * sizeof(float)); // Matrix to hold the rearranged data

#else

  data = (Body*)malloc(N*sizeof(Body));
  Rdata = (Body*)malloc(N*sizeof(Body));
  float *X = (float*)data;
  float *Y = (float*)Rdata;
#endif


  uint* codes = (uint*)malloc(DIM * N * sizeof(uint));
  unsigned long int* mcodes = (unsigned long int*)malloc(N*sizeof(unsigned long int));
  uint long* scodes = (uint long*)malloc(N*sizeof(uint long)); // Buffer for the zcodes
  uint* pointIds = (uint*)malloc(N*sizeof(uint)); // Array to hold the final indexing of pointa to bins
  uint* index = (uint*)malloc(N*sizeof(uint)); //Array to hold the initial index of points 
  uint long* bins = (uint long*)malloc(N*sizeof(uint long*)); // Not needed must be removed
  int* levels = (int*)malloc(N*sizeof(int)); // Not needed must be removed

  // initialize the index
  for(int i=0; i<N; i++){
    index[i] = i;
  }

  /* Generate a 3-dimensional data distribution */

  printf("Generating the dataset...\n");

  gettimeofday (&startwtime, NULL);

#ifndef LARRAY
  create_dataset_T(X, N, dist);

#else

  create_dataset_TL(X, N, dist);

#endif

  gettimeofday (&endwtime, NULL);

  double data_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
  printf("Time to create the data: %f\n", data_time);

  /* Compute the quantization codes */
  int maxlev = 10; // Maximum level of the tree
  int maxheight = 8; // Maximun height of the tree (extra control for debugin)
  int nbins = (1 << maxlev); // number of maximum bins at the leaf level

  double quant1[repeat+1], mortn1[repeat+1], radix1[repeat+1], permu1[repeat+1];
  double quant2[repeat+1], mortn2[repeat+1], radix2[repeat+1], permu2[repeat+1];
  for (int it=0; it<repeat+1; it++) {
    gettimeofday (&startwtime, NULL);
    printf("%d\n",it);

#ifndef LARRAY

    compute_quantization_codes_T(codes, X, N, nbins);

#else

    compute_quantization_codes_TL(codes, X, N, nbins);

#endif

    gettimeofday (&endwtime, NULL);

    double quant_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
    quant1[it] = quant_time;
    printf("Quant: %f\n", quant_time);
  
    /* Compute the morton codes */
  
    gettimeofday (&startwtime, NULL);

    morton_encoding_T(mcodes, codes, N);

    gettimeofday (&endwtime, NULL);

    double mcodes_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
    mortn1[it] = mcodes_time;
    printf("Mortn: %f\n", mcodes_time);
  
    /* Build the Octree */
  
    gettimeofday (&startwtime, NULL);

    bin_sort_radix6(mcodes, scodes, pointIds, index, bins, levels, N, 3*(maxlev-2), 0, 0, 3*(maxlev-maxheight));
    //bin_sort_serial_radix6(mcodes, scodes, pointIds, index, bins, levels, N, 3*(maxlev-2), 0, 0, 3*(maxlev-maxheight));

    gettimeofday (&endwtime, NULL);

    double tbuild_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
    radix1[it] = tbuild_time;
    printf("Radix: %f\n", tbuild_time);
  
    /* Rearrange the data into memory */
  
    gettimeofday (&startwtime, NULL);

#ifndef LARRAY

    rearrange_dataT(Y, X, pointIds, N);

#else

    rearrange_dataTL(Y, X, pointIds, N);

#endif

    gettimeofday (&endwtime, NULL);

    double datare_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
    permu1[it] = datare_time;
    printf("Permu: %f\n", datare_time);

#ifndef LARRAY

    compute_quantization_codes_T(codes, X, N, nbins);

#else

    compute_quantization_codes_TL(codes, X, N, nbins);

#endif

    gettimeofday (&endwtime, NULL);

    quant_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
    quant2[it] = quant_time;
    printf("Quant: %f\n", quant_time);
  
    /* Compute the morton codes */
  
    gettimeofday (&startwtime, NULL);

    morton_encoding_T(mcodes, codes, N);

    gettimeofday (&endwtime, NULL);

    mcodes_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
    mortn2[it] = mcodes_time;
    printf("Mortn: %f\n", mcodes_time);
  
    /* Build the Octree */
  
    gettimeofday (&startwtime, NULL);

    bin_sort_radix6(mcodes, scodes, pointIds, index, bins, levels, N, 3*(maxlev-2), 0, 0, 3*(maxlev-maxheight));
    //bin_sort_serial_radix6(mcodes, scodes, pointIds, index, bins, levels, N, 3*(maxlev-2), 0, 0, 3*(maxlev-maxheight));

    gettimeofday (&endwtime, NULL);

    tbuild_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
    radix2[it] = tbuild_time;
    printf("Radix: %f\n", tbuild_time);
  
    /* Rearrange the data into memory */
  
    gettimeofday (&startwtime, NULL);

#ifndef LARRAY

    rearrange_dataT(Y, X, pointIds, N);

#else

    rearrange_dataTL(Y, X, pointIds, N);

#endif

    gettimeofday (&endwtime, NULL);

    datare_time = (double)((endwtime.tv_usec - startwtime.tv_usec)/1.0e6 + endwtime.tv_sec - startwtime.tv_sec);
    permu2[it] = datare_time;
    printf("Permu: %f\n", datare_time);

  }

  double quant1ave = 0, mortn1ave = 0, radix1ave = 0, permu1ave = 0;
  double quant2ave = 0, mortn2ave = 0, radix2ave = 0, permu2ave = 0;
  for (int it=0; it<repeat; it++) {
    quant1ave += quant1[it+1];
    mortn1ave += mortn1[it+1];
    radix1ave += radix1[it+1];
    permu1ave += permu1[it+1];
    quant2ave += quant2[it+1];
    mortn2ave += mortn2[it+1];
    radix2ave += radix2[it+1];
    permu2ave += permu2[it+1];
  }
  quant1ave /= repeat;
  mortn1ave /= repeat;
  radix1ave /= repeat;
  permu1ave /= repeat;
  quant2ave /= repeat;
  mortn2ave /= repeat;
  radix2ave /= repeat;
  permu2ave /= repeat;
  double quant1std = 0, mortn1std = 0, radix1std = 0, permu1std = 0;
  double quant2std = 0, mortn2std = 0, radix2std = 0, permu2std = 0;
  for (int it=0; it<repeat; it++) {
    quant1std = (quant1[it+1] - quant1ave) * (quant1[it+1] - quant1ave);
    mortn1std = (mortn1[it+1] - mortn1ave) * (mortn1[it+1] - mortn1ave);
    radix1std = (radix1[it+1] - radix1ave) * (radix1[it+1] - radix1ave);
    permu1std = (permu1[it+1] - permu1ave) * (permu1[it+1] - permu1ave);
    quant2std = (quant2[it+1] - quant2ave) * (quant2[it+1] - quant2ave);
    mortn2std = (mortn2[it+1] - mortn2ave) * (mortn2[it+1] - mortn2ave);
    radix2std = (radix2[it+1] - radix2ave) * (radix2[it+1] - radix2ave);
    permu2std = (permu2[it+1] - permu2ave) * (permu2[it+1] - permu2ave);
  }
  quant1std /= repeat;
  mortn1std /= repeat;
  radix1std /= repeat;
  permu1std /= repeat;
  quant2std /= repeat;
  mortn2std /= repeat;
  radix2std /= repeat;
  permu2std /= repeat;
  quant1std = sqrt(quant1std);
  mortn1std = sqrt(mortn1std);
  radix1std = sqrt(radix1std);
  permu1std = sqrt(permu1std);
  quant2std = sqrt(quant2std);
  mortn2std = sqrt(mortn2std);
  radix2std = sqrt(radix2std);
  permu2std = sqrt(permu2std);

  printf("Quant1: %lf+-%lf Mortn1: %lf+-%lf Radix1: %lf+-%lf Permu1: %lf+-%lf\n",
	 quant1ave,quant1std,mortn1ave,mortn1std,radix1ave,radix1std,permu1ave,permu1std);
  printf("Quant2: %lf+-%lf Mortn2: %lf+-%lf Radix2: %lf+-%lf Permu2: %lf+-%lf\n",
	 quant2ave,quant2std,mortn2ave,mortn2std,radix2ave,radix2std,permu2ave,permu2std);
  
#ifndef LARRAY
  free(X);
  free(Y);
#else
  free(data);
  free(Rdata);
#endif
  free(codes);
  free(mcodes);
  free(scodes);
  free(pointIds);
  free(index);
  free(bins);
  free(levels);
}
