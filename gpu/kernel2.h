#pragma once

namespace {
  template<int nx, int ny, int nz>
  struct Index {
    static const int      I = Index<nx,ny+1,nz-1>::I + 1;
    static const uint64_t F = Index<nx,ny,nz-1>::F * nz;
    static __host__ __device__ __forceinline__
    float power(const fvec3 &dX) {
      return Index<nx,ny,nz-1>::power(dX) * dX[2] / nz;
    }
  };

  template<int nx, int ny>
  struct Index<nx,ny,0> {
    static const int      I = Index<nx+1,0,ny-1>::I + 1;
    static const uint64_t F = Index<nx,ny-1,0>::F * ny;
    static __host__ __device__ __forceinline__
    float power(const fvec3 &dX) {
      return Index<nx,ny-1,0>::power(dX) * dX[1] / ny;
    }
  };

  template<int nx>
  struct Index<nx,0,0> {
    static const int      I = Index<0,0,nx-1>::I + 1;
    static const uint64_t F = Index<nx-1,0,0>::F * nx;
    static __host__ __device__ __forceinline__
    float power(const fvec3 &dX) {
      return Index<nx-1,0,0>::power(dX) * dX[0] / nx;
    }
  };

  template<>
  struct Index<0,0,0> {
    static const int      I = 0;
    static const uint64_t F = 1;
    static __host__ __device__ __forceinline__
    float power(const fvec3&) { return 1.0f; }
  };


  template<int n, int kx, int ky , int kz, int d>
  struct DerivativeTerm {
    static const int coef = 1 - 2 * n;
    static __host__ __device__ __forceinline__
    float kernel(const fvecP &C, const fvec3 &dX) {
      return coef * dX[d] * C[Index<kx,ky,kz>::I];
    }
  };

  template<int n, int kx, int ky , int kz>
  struct DerivativeTerm<n,kx,ky,kz,-1> {
    static const int coef = 1 - n;
    static __host__ __device__ __forceinline__
      float kernel(const fvecP &C, const fvec3&) {
      return coef * C[Index<kx,ky,kz>::I];
    }
  };


  template<int nx, int ny, int nz, int kx=nx, int ky=ny, int kz=nz, int flag=5>
  struct DerivativeSum {
    static const int nextflag = 5 - (kz < nz || kz == 1);
    static const int dim = kz == (nz-1) ? -1 : 2;
    static const int n = nx + ny + nz;
    static __host__ __device__ __forceinline__
    float loop(const fvecP &C, const fvec3 &dX) {
      return DerivativeSum<nx,ny,nz,nx,ny,kz-1,nextflag>::loop(C,dX)
	+ DerivativeTerm<n,nx,ny,kz-1,dim>::kernel(C,dX);
    }
  };

  template<int nx, int ny, int nz, int kx, int ky, int kz>
  struct DerivativeSum<nx,ny,nz,kx,ky,kz,4> {
    static const int nextflag = 3 - (ny == 0);
    static __host__ __device__ __forceinline__
    float loop(const fvecP &C, const fvec3 &dX) {
      return DerivativeSum<nx,ny,nz,nx,ny,nz,nextflag>::loop(C,dX);
    }
  };

  template<int nx, int ny, int nz, int kx, int ky, int kz>
  struct DerivativeSum<nx,ny,nz,kx,ky,kz,3> {
    static const int nextflag = 3 - (ky < ny || ky == 1);
    static const int dim = ky == (ny-1) ? -1 : 1;
    static const int n = nx + ny + nz;
    static __host__ __device__ __forceinline__
    float loop(const fvecP &C, const fvec3 &dX) {
      return DerivativeSum<nx,ny,nz,nx,ky-1,nz,nextflag>::loop(C,dX)
	+ DerivativeTerm<n,nx,ky-1,nz,dim>::kernel(C,dX);
    }
  };

  template<int nx, int ny, int nz, int kx, int ky, int kz>
  struct DerivativeSum<nx,ny,nz,kx,ky,kz,2> {
    static const int nextflag = 1 - (nx == 0);
    static __host__ __device__ __forceinline__
    float loop(const fvecP &C, const fvec3 &dX) {
      return DerivativeSum<nx,ny,nz,nx,ny,nz,nextflag>::loop(C,dX);
    }
  };

  template<int nx, int ny, int nz, int kx, int ky, int kz>
  struct DerivativeSum<nx,ny,nz,kx,ky,kz,1> {
    static const int nextflag = 1 - (kx < nx || kx == 1);
    static const int dim = kx == (nx-1) ? -1 : 0;
    static const int n = nx + ny + nz;
    static __host__ __device__ __forceinline__
    float loop(const fvecP &C, const fvec3 &dX) {
      return DerivativeSum<nx,ny,nz,kx-1,ny,nz,nextflag>::loop(C,dX)
	+ DerivativeTerm<n,kx-1,ny,nz,dim>::kernel(C,dX);
    }
  };

  template<int nx, int ny, int nz, int kx, int ky, int kz>
  struct DerivativeSum<nx,ny,nz,kx,ky,kz,0> {
    static __host__ __device__ __forceinline__
    float loop(const fvecP&, const fvec3&) {
      return 0;
    }
  };

  template<int nx, int ny, int nz, int kx, int ky>
  struct DerivativeSum<nx,ny,nz,kx,ky,0,5> {
    static __host__ __device__ __forceinline__
    float loop(const fvecP &C, const fvec3 &dX) {
      return DerivativeSum<nx,ny,nz,nx,ny,0,4>::loop(C,dX);
    }
  };


  template<int nx, int ny, int nz, int kx=nx, int ky=ny, int kz=nz>
  struct MultipoleSum {
    static __host__ __device__ __forceinline__
    float kernel(const fvec3 &dX, const fvecP &M) {
      return MultipoleSum<nx,ny,nz,kx,ky,kz-1>::kernel(dX,M)
	+ Index<nx-kx,ny-ky,nz-kz>::power(dX)*M[Index<kx,ky,kz>::I];
    }
  };

  template<int nx, int ny, int nz, int kx, int ky>
  struct MultipoleSum<nx,ny,nz,kx,ky,0> {
    static __host__ __device__ __forceinline__
    float kernel(const fvec3 &dX, const fvecP &M) {
      return MultipoleSum<nx,ny,nz,kx,ky-1,nz>::kernel(dX,M)
	+ Index<nx-kx,ny-ky,nz>::power(dX)*M[Index<kx,ky,0>::I];
    }
  };

  template<int nx, int ny, int nz, int kx>
  struct MultipoleSum<nx,ny,nz,kx,0,0> {
    static __host__ __device__ __forceinline__
    float kernel(const fvec3 &dX, const fvecP &M) {
      return MultipoleSum<nx,ny,nz,kx-1,ny,nz>::kernel(dX,M)
	+ Index<nx-kx,ny,nz>::power(dX)*M[Index<kx,0,0>::I];
    }
  };

  template<int nx, int ny, int nz>
  struct MultipoleSum<nx,ny,nz,0,0,0> {
    static __host__ __device__ __forceinline__
    float kernel(const fvec3 &dX, const fvecP &M) {
      return Index<nx,ny,nz>::power(dX)*M[Index<0,0,0>::I];
    }
  };


  template<int nx, int ny, int nz>
  struct Kernels {
    static const int n = nx + ny + nz;
    static const int x = nx > 0;
    static const int y = ny > 0;
    static const int z = nz > 0;
    static __host__ __device__ __forceinline__
    void P2M(fvecP &M, const fvec3 &dX) {
      Kernels<nx,ny+1,nz-1>::P2M(M,dX);
      M[Index<nx,ny,nz>::I] = Index<nx,ny,nz>::power(dX) * M[0];
    }
    static __host__ __device__ __forceinline__
    void M2M(fvecP &MI, const fvec3 &dX, const fvecP &MJ) {
      Kernels<nx,ny+1,nz-1>::M2M(MI,dX,MJ);
      MI[Index<nx,ny,nz>::I] += MultipoleSum<nx,ny,nz>::kernel(dX,MJ);
    }
    static __host__ __device__ __forceinline__
    void M2P(fvec4 &TRG, fvecP &C, const fvec3 &dX, const float &invR2, const fvecP &M) {
      Kernels<nx,ny+1,nz-1>::M2P(TRG,C,dX,invR2,M);
      C[Index<nx,ny,nz>::I] = DerivativeSum<nx,ny,nz>::loop(C,dX) / n * invR2;
      TRG[0] -= M[Index<nx,ny,nz>::I] * C[Index<nx,ny,nz>::I] * Index<nx,ny,nz>::F;
      TRG[1] += M[Index<(nx-1)*x,ny,nz>::I] * C[Index<nx,ny,nz>::I] * Index<nx,ny,nz>::F * x;
      TRG[2] += M[Index<nx,(ny-1)*y,nz>::I] * C[Index<nx,ny,nz>::I] * Index<nx,ny,nz>::F * y;
      TRG[3] += M[Index<nx,ny,(nz-1)*z>::I] * C[Index<nx,ny,nz>::I] * Index<nx,ny,nz>::F * z;
    }
  };

  template<int nx, int ny>
  struct Kernels<nx,ny,0> {
    static const int n = nx + ny;
    static const int x = nx > 0;
    static const int y = ny > 0;
    static __host__ __device__ __forceinline__
    void P2M(fvecP &M, const fvec3 &dX) {
      Kernels<nx+1,0,ny-1>::P2M(M,dX);
      M[Index<nx,ny,0>::I] = Index<nx,ny,0>::power(dX) * M[0];
    }
    static __host__ __device__ __forceinline__
    void M2M(fvecP &MI, const fvec3 &dX, const fvecP &MJ) {
      Kernels<nx+1,0,ny-1>::M2M(MI,dX,MJ);
      MI[Index<nx,ny,0>::I] += MultipoleSum<nx,ny,0>::kernel(dX,MJ);
    }
    static __host__ __device__ __forceinline__
    void M2P(fvec4 &TRG, fvecP &C, const fvec3 &dX, const float &invR2, const fvecP &M) {
      Kernels<nx+1,0,ny-1>::M2P(TRG,C,dX,invR2,M);
      C[Index<nx,ny,0>::I] = DerivativeSum<nx,ny,0>::loop(C,dX) / n * invR2;
      TRG[0] -= M[Index<nx,ny,0>::I] * C[Index<nx,ny,0>::I] * Index<nx,ny,0>::F;
      TRG[1] += M[Index<(nx-1)*x,ny,0>::I] * C[Index<nx,ny,0>::I] * Index<nx,ny,0>::F * x;
      TRG[2] += M[Index<nx,(ny-1)*y,0>::I] * C[Index<nx,ny,0>::I] * Index<nx,ny,0>::F * y;
    }
  };

  template<int nx>
  struct Kernels<nx,0,0> {
    static const int n = nx;
    static __host__ __device__ __forceinline__
    void P2M(fvecP &M, const fvec3 &dX) {
      Kernels<0,0,nx-1>::P2M(M,dX);
      M[Index<nx,0,0>::I] = Index<nx,0,0>::power(dX) * M[0];
    }
    static __host__ __device__ __forceinline__
    void M2M(fvecP &MI, const fvec3 &dX, const fvecP &MJ) {
      Kernels<0,0,nx-1>::M2M(MI,dX,MJ);
      MI[Index<nx,0,0>::I] += MultipoleSum<nx,0,0>::kernel(dX,MJ);
    }
    static __host__ __device__ __forceinline__
    void M2P(fvec4 &TRG, fvecP &C, const fvec3 &dX, const float &invR2, const fvecP &M) {
      Kernels<0,0,nx-1>::M2P(TRG,C,dX,invR2,M);
      C[Index<nx,0,0>::I] = DerivativeSum<nx,0,0>::loop(C,dX) / n * invR2;
      TRG[0] -= M[Index<nx,0,0>::I] * C[Index<nx,0,0>::I] * Index<nx,0,0>::F;
      TRG[1] += M[Index<nx-1,0,0>::I] * C[Index<nx,0,0>::I] * Index<nx,0,0>::F;
    }
  };

  template<>
  struct Kernels<0,0,0> {
    static __host__ __device__ __forceinline__
    void P2M(fvecP &M, const fvec3 &dX) {}
    static __host__ __device__ __forceinline__
    void M2M(fvecP &MI, const fvec3 &dX, const fvecP &MJ) {
      MI[Index<0,0,0>::I] += MultipoleSum<0,0,0>::kernel(dX,MJ);
    }
    static __host__ __device__ __forceinline__
    void M2P(fvec4 &TRG, fvecP &C, const fvec3&, const float&, const fvecP&) {
      TRG[0] -= C[0];
    }
  };


  __device__ __forceinline__
  void P2M(const int begin,
	   const int end,
	   const fvec4 center,
	   fvecP & Mi) {
    for (int i=begin; i<end; i++) {
      fvec4 body = tex1Dfetch(texBody,i);
      fvec3 dX = make_fvec3(center - body);
      fvecP M;
      M[0] = body[3];
      Kernels<0,0,P-1>::P2M(M,dX);
      Mi += M;
    }
  }

  __device__ __forceinline__
  void M2M(const int begin,
	   const int end,
	   const fvec4 Xi,
	   fvec4 * sourceCenter,
	   fvec4 * Multipole,
	   fvecP & Mi) {
    for (int i=begin; i<end; i++) {
      fvecP Mj = *(fvecP*)&Multipole[NVEC4*i];
      fvec4 Xj = sourceCenter[i];
      fvec3 dX = make_fvec3(Xi - Xj);
      Kernels<0,0,P-1>::M2M(Mi,dX,Mj);
    }
  }

  __device__ __forceinline__
  fvec4 P2P(fvec4 acc,
	    const fvec3 pos_i,
	    const fvec3 pos_j,
	    const float q_j,
	    const float EPS2) {
    fvec3 dX = pos_j - pos_i;
    const float R2 = norm(dX) + EPS2;
    const float invR = rsqrtf(R2);
    const float invR2 = invR * invR;
    const float invR1 = q_j * invR;
    dX *= invR1 * invR2;
    acc[0] -= invR1;
    acc[1] += dX[0];
    acc[2] += dX[1];
    acc[3] += dX[2];
    return acc;
  }

  __device__ __forceinline__
  fvec4 M2P(fvec4 acc,
	    const fvec3 & pos_i,
	    const fvec3 & pos_j,
	    fvecP & M,
	    float EPS2) {
    const fvec3 dX = pos_i - pos_j;
    const float R2 = norm(dX) + EPS2;
    const float invR = rsqrtf(R2);
#if 0
    const float invR2 = -invR * invR;
    const float invR1 = M[0] * invR;
    const float invR3 = invR2 * invR1;
    const float invR5 = 3 * invR2 * invR3;
    const float invR7 = 5 * invR2 * invR5;
    const float q11 = M[4];
    const float q12 = 0.5f * M[5];
    const float q13 = 0.5f * M[6];
    const float q22 = M[7];
    const float q23 = 0.5f * M[8];
    const float q33 = M[9];
    const float q = q11 + q22 + q33;
    fvec3 qR;
    qR[0] = q11 * dX[0] + q12 * dX[1] + q13 * dX[2];
    qR[1] = q12 * dX[0] + q22 * dX[1] + q23 * dX[2];
    qR[2] = q13 * dX[0] + q23 * dX[1] + q33 * dX[2];
    const float qRR = qR[0] * dX[0] + qR[1] * dX[1] + qR[2] * dX[2];
    acc[0] -= invR1 + invR3 * q + invR5 * qRR;
    const float C = invR3 + invR5 * q + invR7 * qRR;
    acc[1] += C * dX[0] + 2 * invR5 * qR[0];
    acc[2] += C * dX[1] + 2 * invR5 * qR[1];
    acc[3] += C * dX[2] + 2 * invR5 * qR[2];
#else
    const float M0 = M[0];
    const float invR1 = M0 * invR;
    float invR2 = invR * invR;
    fvecP C;
    C[0] = invR1;
    M[0] = 1;
    Kernels<0,0,P-1>::M2P(acc,C,dX,invR2,M);
    M[0] = M0;
#endif
    return acc;
  }
}
