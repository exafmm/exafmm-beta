/*
Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
#define KERNEL
#include "kernel.h"
#undef KERNEL
#  define CUDA_SAFE_CALL( call) {                                    \
    cudaError err = call;                                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error in file '%s' in line %i : %s.\n",        \
                __FILE__, __LINE__, cudaGetErrorString( err) );              \
        exit(EXIT_FAILURE);                                                  \
    } }

#ifdef _DEBUG
#  define CUT_CHECK_ERROR(errorMessage) {                                    \
    cudaError_t err = cudaGetLastError();                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
                errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
        exit(EXIT_FAILURE);                                                  \
    }                                                                        \
    err = CUT_DEVICE_SYNCHRONIZE();                                           \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
                errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
        exit(EXIT_FAILURE);                                                  \
    }                                                                        \
    }
#else
#  define CUT_CHECK_ERROR(errorMessage) {                                    \
    cudaError_t err = cudaGetLastError();                                    \
    if( cudaSuccess != err) {                                                \
        fprintf(stderr, "Cuda error: %s in file '%s' in line %i : %s.\n",    \
                errorMessage, __FILE__, __LINE__, cudaGetErrorString( err) );\
        exit(EXIT_FAILURE);                                                  \
    }                                                                        \
    }
#endif

__device__ __constant__ gpureal constDevc[514];                 // Constants on device

void Kernel<Stokes>::initialize()
{
    startTimer("Init GPU     ");                                  // Start timer
    cudaThreadExit();                                             // Exit GPU thread
    cudaSetDevice(DEVICE);                                        // Set GPU device
    cudaThreadSynchronize();                                      // Sync GPU threads
#ifdef CUDA_4_1
    cudaSetDeviceFlags(cudaDeviceMapHost);
#endif
    stopTimer("Init GPU     ", MPIRANK == 0);                     // Stop timer & print
    eraseTimer("Init GPU     ");                                  // Erase timer
}

void Kernel<Stokes>::finalize()
{
}

void Kernel<Stokes>::allocate()
{
    cudaThreadSynchronize();
    startTimer("cudaMalloc   ");
#ifdef CUDA_4_1
    if ( keysHost.size() > keysDevcSize )
    {
        CUDA_SAFE_CALL(cudaHostRegister(&keysHost[0], sizeof(int) * keysHost.size(), cudaHostRegisterMapped));
        CUDA_SAFE_CALL(cudaHostGetDevicePointer((void **) &keysDevc, (void *)&keysHost[0], 0));
        keysDevcSize = keysHost.size();
    }
    if ( rangeHost.size() > rangeDevcSize )
    {
        CUDA_SAFE_CALL(cudaHostRegister(&rangeHost[0], sizeof(int) * rangeHost.size(), cudaHostRegisterMapped));
        CUDA_SAFE_CALL(cudaHostGetDevicePointer((void **) &rangeDevc, (void *)&rangeHost[0], 0));
        rangeDevcSize = rangeHost.size();
    }
    if ( sourceHost.size() > sourceDevcSize )
    {
        CUDA_SAFE_CALL(cudaHostRegister(&sourceDevc[0], sizeof(gpureal) * sourceHost.size(), cudaHostRegisterMapped));
        CUDA_SAFE_CALL(cudaHostGetDevicePointer((void **) &sourceDevc, (void *)&sourceDevc[0], 0));
        sourceDevcSize = sourceHost.size();
    }
    if ( targetHost.size() > targetDevcSize )
    {
        CUDA_SAFE_CALL(cudaHostRegister(&targetHost[0], sizeof(gpureal) * targetHost.size(), cudaHostRegisterMapped));
        CUDA_SAFE_CALL(cudaHostGetDevicePointer((void **) &targetDevc, (void *)&targetHost[0], 0));
        targetDevcSize = targetHost.size();
    }
#else
    if ( keysHost.size() > keysDevcSize )
    {
        if ( keysDevcSize != 0 ) CUDA_SAFE_CALL(cudaFree(keysDevc));
        CUDA_SAFE_CALL(cudaMalloc( (void**) &keysDevc, keysHost.size()*sizeof(int) ));
        keysDevcSize = keysHost.size();
    }
    if ( rangeHost.size() > rangeDevcSize )
    {
        if ( rangeDevcSize != 0 ) CUDA_SAFE_CALL(cudaFree(rangeDevc));
        CUDA_SAFE_CALL(cudaMalloc( (void**) &rangeDevc, rangeHost.size()*sizeof(int) ));
        rangeDevcSize = rangeHost.size();
    }
    if ( sourceHost.size() > sourceDevcSize )
    {
        if ( sourceDevcSize != 0 ) CUDA_SAFE_CALL(cudaFree(sourceDevc));
        CUDA_SAFE_CALL(cudaMalloc( (void**) &sourceDevc, sourceHost.size()*sizeof(gpureal) ));
        sourceDevcSize = sourceHost.size();
    }
    if ( targetHost.size() > targetDevcSize )
    {
        if ( targetDevcSize != 0 ) CUDA_SAFE_CALL(cudaFree(targetDevc));
        CUDA_SAFE_CALL(cudaMalloc( (void**) &targetDevc, targetHost.size()*sizeof(gpureal) ));
        targetDevcSize = targetHost.size();
    }
#endif
    cudaThreadSynchronize();
    stopTimer("cudaMalloc   ");
}


void Kernel<Stokes>::hostToDevice()
{
    cudaThreadSynchronize();
    startTimer("cudaMemcpy   ");
#ifndef CUDA_4_1
    CUDA_SAFE_CALL(cudaMemcpy(keysDevc,  &keysHost[0],  keysHost.size()*sizeof(int), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(rangeDevc, &rangeHost[0], rangeHost.size()*sizeof(int), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(sourceDevc, &sourceHost[0], sourceHost.size()*sizeof(gpureal), cudaMemcpyHostToDevice));
    CUDA_SAFE_CALL(cudaMemcpy(targetDevc, &targetHost[0], targetHost.size()*sizeof(gpureal), cudaMemcpyHostToDevice));
#endif
    CUDA_SAFE_CALL(cudaMemcpyToSymbol(constDevc, &constHost[0], constHost.size()*sizeof(gpureal)));
    cudaThreadSynchronize();
    stopTimer("cudaMemcpy   ");
}

void Kernel<Stokes>::deviceToHost()
{
    cudaThreadSynchronize();
    startTimer("cudaMemcpy   ");
#ifdef CUDA_4_1
    CUDA_SAFE_CALL(cudaHostUnregister(&keysHost[0]));
    CUDA_SAFE_CALL(cudaHostUnregister(&rangeHost[0]));
    CUDA_SAFE_CALL(cudaHostUnregister(&sourceHost[0]));
    CUDA_SAFE_CALL(cudaHostUnregister(&targetHost[0]));
#else
    CUDA_SAFE_CALL(cudaMemcpy(&targetHost[0], targetDevc, targetHost.size()*sizeof(gpureal), cudaMemcpyDeviceToHost));
#endif
    cudaThreadSynchronize();
    stopTimer("cudaMemcpy   ");
}

__device__ void StokesP2P_core(gpureal *target, gpureal *targetX, gpureal *sourceShrd, float3 d, int i, float delta)
{
    d.x += targetX[0];
    d.x -= sourceShrd[6*i+0];
    d.y += targetX[1];
    d.y -= sourceShrd[6*i+1];
    d.z += targetX[2];
    d.z -= sourceShrd[6*i+2];

    float3 force = {sourceShrd[6*i+3], sourceShrd[6*i+4], sourceShrd[6*i+5]};

    float r2 = d.x * d.x + d.y * d.y + d.z * d.z;
    float d2 = delta * delta;
    float R1 = r2 + d2;
    float R2 = R1 + d2;
    float invR = 1.0 / R1;
    float H = sqrt(invR) * invR;

    float fdx =  force.x * d.x + force.y * d.y + force.z * d.z;

    target[0] += H * (force.x * R2 + fdx * d.x);
    target[1] += H * (force.y * R2 + fdx * d.y);
    target[2] += H * (force.z * R2 + fdx * d.z);

}

__global__ void StokesP2P_GPU(int *keysGlob, int *rangeGlob, gpureal *targetGlob, gpureal *sourceGlob, float delta)
{
    int keys = keysGlob[blockIdx.x];
    int numList = rangeGlob[keys];
    gpureal D0 = -constDevc[0];
    gpureal targetX[3];
    gpureal target[4] = {0, 0, 0, 0};
    __shared__ gpureal sourceShrd[6*THREADS];
    int itarget = blockIdx.x * THREADS + threadIdx.x;
    targetX[0] = targetGlob[4*itarget+0];
    targetX[1] = targetGlob[4*itarget+1];
    targetX[2] = targetGlob[4*itarget+2];
    for ( int ilist = 0; ilist < numList; ++ilist )
    {
        int begin     = rangeGlob[keys+3*ilist+1];
        int size      = rangeGlob[keys+3*ilist+2];
        int Iperiodic = rangeGlob[keys+3*ilist+3];
        for ( int iblok = 0; iblok < (size - 1) / THREADS; ++iblok ){
            int isource = begin + iblok * THREADS + threadIdx.x;
            __syncthreads();
            sourceShrd[6*threadIdx.x+0] = sourceGlob[6*isource+0];
            sourceShrd[6*threadIdx.x+1] = sourceGlob[6*isource+1];
            sourceShrd[6*threadIdx.x+2] = sourceGlob[6*isource+2];
            sourceShrd[6*threadIdx.x+3] = sourceGlob[6*isource+3];
            sourceShrd[6*threadIdx.x+4] = sourceGlob[6*isource+4];
            sourceShrd[6*threadIdx.x+5] = sourceGlob[6*isource+5];
            __syncthreads();
            int I = 0;
            for ( int ix = -1; ix <= 1; ++ix ){
                for ( int iy = -1; iy <= 1; ++iy ){
                    for ( int iz = -1; iz <= 1; ++iz, ++I ){
                        if ( Iperiodic & (1 << I) )
                        {
                            float3 d;
                            d.x = ix * D0;
                            d.y = iy * D0;
                            d.z = iz * D0;
#pragma unroll 64
                            for ( int i = 0; i < THREADS; ++i ){
                                StokesP2P_core(target, targetX, sourceShrd, d, i, delta);
                            }
                        }
                    }
                }
            }
        }
        int iblok = (size - 1) / THREADS;
        int isource = begin + iblok * THREADS + threadIdx.x;
        __syncthreads();
        if ( threadIdx.x < size - iblok * THREADS )
        {
            sourceShrd[6*threadIdx.x+0] = sourceGlob[6*isource+0];
            sourceShrd[6*threadIdx.x+1] = sourceGlob[6*isource+1];
            sourceShrd[6*threadIdx.x+2] = sourceGlob[6*isource+2];
            sourceShrd[6*threadIdx.x+3] = sourceGlob[6*isource+3];
            sourceShrd[6*threadIdx.x+4] = sourceGlob[6*isource+4];
            sourceShrd[6*threadIdx.x+5] = sourceGlob[6*isource+5];
        }
        __syncthreads();
        int I = 0;
        int icounter = 0;
        for ( int ix = -1; ix <= 1; ++ix ){
            for ( int iy = -1; iy <= 1; ++iy ){
                for ( int iz = -1; iz <= 1; ++iz, ++I ){
                    if ( Iperiodic & (1 << I) )
                    {
                        icounter++;
                        float3 d;
                        d.x = ix * D0;
                        d.y = iy * D0;
                        d.z = iz * D0;
                        for ( int i = 0; i < size - iblok*THREADS; ++i )
                        {
                            StokesP2P_core(target, targetX, sourceShrd, d, i, delta);
                        }
                    }
                }
            }
        }
    }
    targetGlob[4*itarget+0] = target[0];
    targetGlob[4*itarget+1] = target[1];
    targetGlob[4*itarget+2] = target[2];
    targetGlob[4*itarget+3] = target[3];
}

void Kernel<Stokes>::P2P()
{
    cudaThreadSynchronize();
    startTimer("P2P GPUkernel");
    int numBlocks = keysHost.size();
    if ( numBlocks != 0 )
    {
        StokesP2P_GPU <<< numBlocks, THREADS >>>(keysDevc, rangeDevc, targetDevc, sourceDevc, delta);
    }
    CUT_CHECK_ERROR("Kernel execution failed");
    cudaThreadSynchronize();
    stopTimer("P2P GPUkernel");
}


