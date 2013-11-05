#pragma once

namespace {
  template<int nx, int ny, int nz>
  struct Index {
    static const int                I = Index<nx,ny+1,nz-1>::I + 1;
    static const unsigned long long F = Index<nx,ny,nz-1>::F * nz;
  };

  template<int nx, int ny>
  struct Index<nx,ny,0> {
    static const int                I = Index<nx+1,0,ny-1>::I + 1;
    static const unsigned long long F = Index<nx,ny-1,0>::F * ny;
  };

  template<int nx>
  struct Index<nx,0,0> {
    static const int                I = Index<0,0,nx-1>::I + 1;
    static const unsigned long long F = Index<nx-1,0,0>::F * nx;
  };

  template<>
  struct Index<0,0,0> {
    static const int                I = 0;
    static const unsigned long long F = 1;
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
    float kernel(const fvecP &C, const fvecP &M) {
      return MultipoleSum<nx,ny,nz,kx,ky,kz-1>::kernel(C,M)
	+ C[Index<nx-kx,ny-ky,nz-kz>::I]*M[Index<kx,ky,kz>::I];
    }
  };

  template<int nx, int ny, int nz, int kx, int ky>
  struct MultipoleSum<nx,ny,nz,kx,ky,0> {
    static __host__ __device__ __forceinline__
    float kernel(const fvecP &C, const fvecP &M) {
      return MultipoleSum<nx,ny,nz,kx,ky-1,nz>::kernel(C,M)
	+ C[Index<nx-kx,ny-ky,nz>::I]*M[Index<kx,ky,0>::I];
    }
  };

  template<int nx, int ny, int nz, int kx>
  struct MultipoleSum<nx,ny,nz,kx,0,0> {
    static __host__ __device__ __forceinline__
    float kernel(const fvecP &C, const fvecP &M) {
      return MultipoleSum<nx,ny,nz,kx-1,ny,nz>::kernel(C,M)
	+ C[Index<nx-kx,ny,nz>::I]*M[Index<kx,0,0>::I];
    }
  };

  template<int nx, int ny, int nz>
  struct MultipoleSum<nx,ny,nz,0,0,0> {
    static __host__ __device__ __forceinline__
    float kernel(const fvecP&, const fvecP&) { return 0; }
  };


  template<int nx, int ny, int nz, int kx=0, int ky=0, int kz=P-1-nx-ny-nz>
    struct LocalSum {
      static __host__ __device__ __forceinline__
      float kernel(const fvecP &C, const fvecP &M) {
	return LocalSum<nx,ny,nz,kx,ky+1,kz-1>::kernel(C,M)
	  + M[Index<kx,ky,kz>::I] * C[Index<nx+kx,ny+ky,nz+kz>::I];
      }
    };

  template<int nx, int ny, int nz, int kx, int ky>
  struct LocalSum<nx,ny,nz,kx,ky,0> {
    static __host__ __device__ __forceinline__
    float kernel(const fvecP &C, const fvecP &M) {
      return LocalSum<nx,ny,nz,kx+1,0,ky-1>::kernel(C,M)
	+ M[Index<kx,ky,0>::I] * C[Index<nx+kx,ny+ky,nz>::I];
    }
  };

  template<int nx, int ny, int nz, int kx>
  struct LocalSum<nx,ny,nz,kx,0,0> {
    static __host__ __device__ __forceinline__
    float kernel(const fvecP &C, const fvecP &M) {
      return LocalSum<nx,ny,nz,0,0,kx-1>::kernel(C,M)
	+ M[Index<kx,0,0>::I] * C[Index<nx+kx,ny,nz>::I];
    }
  };

  template<int nx, int ny, int nz>
  struct LocalSum<nx,ny,nz,0,0,0> {
    static __host__ __device__ __forceinline__
    float kernel(const fvecP&, const fvecP&) { return 0; }
  };


  template<int nx, int ny, int nz>
  struct Kernels {
    static const int n = nx + ny + nz;
    static __host__ __device__ __forceinline__
    void power(fvecP &C, const fvec3 &dX) {
      Kernels<nx,ny+1,nz-1>::power(C,dX);
      C[Index<nx,ny,nz>::I] = C[Index<nx,ny,nz-1>::I] * dX[2] / nz;
    }
    static __host__ __device__ __forceinline__
    void derivative(fvecP &C, const fvec3 &dX, const float &invR2) {
      Kernels<nx,ny+1,nz-1>::derivative(C,dX,invR2);
      C[Index<nx,ny,nz>::I] = DerivativeSum<nx,ny,nz>::loop(C,dX) / n * invR2;
    }
    static __host__ __device__ __forceinline__
    void scale(fvecP &C) {
      Kernels<nx,ny+1,nz-1>::scale(C);
      C[Index<nx,ny,nz>::I] *= Index<nx,ny,nz>::F;
    }
    static __host__ __device__ __forceinline__
    void M2M(fvecP &MI, const fvecP &C, const fvecP &MJ) {
      Kernels<nx,ny+1,nz-1>::M2M(MI,C,MJ);
      MI[Index<nx,ny,nz>::I] += MultipoleSum<nx,ny,nz>::kernel(C,MJ);
    }
    static __host__ __device__ __forceinline__
    void M2P(fvec4 &TRG, const fvecP &C, const fvecP &M) {
      Kernels<nx,ny+1,nz-1>::M2P(TRG,C,M);
      TRG[Index<nx,ny,nz>::I] += LocalSum<nx,ny,nz>::kernel(C,M);
    }
  };

  template<int nx, int ny>
  struct Kernels<nx,ny,0> {
    static const int n = nx + ny;
    static __host__ __device__ __forceinline__
    void power(fvecP &C, const fvec3 &dX) {
      Kernels<nx+1,0,ny-1>::power(C,dX);
      C[Index<nx,ny,0>::I] = C[Index<nx,ny-1,0>::I] * dX[1] / ny;
    }
    static __host__ __device__ __forceinline__
    void derivative(fvecP &C, const fvec3 &dX, const float &invR2) {
      Kernels<nx+1,0,ny-1>::derivative(C,dX,invR2);
      C[Index<nx,ny,0>::I] = DerivativeSum<nx,ny,0>::loop(C,dX) / n * invR2;
    }
    static __host__ __device__ __forceinline__
    void scale(fvecP &C) {
      Kernels<nx+1,0,ny-1>::scale(C);
      C[Index<nx,ny,0>::I] *= Index<nx,ny,0>::F;
    }
    static __host__ __device__ __forceinline__
    void M2M(fvecP &MI, const fvecP &C, const fvecP &MJ) {
      Kernels<nx+1,0,ny-1>::M2M(MI,C,MJ);
      MI[Index<nx,ny,0>::I] += MultipoleSum<nx,ny,0>::kernel(C,MJ);
    }
    static __host__ __device__ __forceinline__
    void M2P(fvec4 &TRG, const fvecP &C, const fvecP &M) {
      Kernels<nx+1,0,ny-1>::M2P(TRG,C,M);
      TRG[Index<nx,ny,0>::I] += LocalSum<nx,ny,0>::kernel(C,M);
    }
  };

  template<int nx>
  struct Kernels<nx,0,0> {
    static const int n = nx;
    static __host__ __device__ __forceinline__
    void power(fvecP &C, const fvec3 &dX) {
      Kernels<0,0,nx-1>::power(C,dX);
      C[Index<nx,0,0>::I] = C[Index<nx-1,0,0>::I] * dX[0] / nx;
    }
    static __host__ __device__ __forceinline__
    void derivative(fvecP &C, const fvec3 &dX, const float &invR2) {
      Kernels<0,0,nx-1>::derivative(C,dX,invR2);
      C[Index<nx,0,0>::I] = DerivativeSum<nx,0,0>::loop(C,dX) / n * invR2;
    }
    static __host__ __device__ __forceinline__
    void scale(fvecP &C) {
      Kernels<0,0,nx-1>::scale(C);
      C[Index<nx,0,0>::I] *= Index<nx,0,0>::F;
    }

    static __host__ __device__ __forceinline__
    void M2M(fvecP &MI, const fvecP &C, const fvecP &MJ) {
      Kernels<0,0,nx-1>::M2M(MI,C,MJ);
      MI[Index<nx,0,0>::I] += MultipoleSum<nx,0,0>::kernel(C,MJ);
    }
    static __host__ __device__ __forceinline__
    void M2P(fvec4 &TRG, const fvecP &C, const fvecP &M) {
      Kernels<0,0,nx-1>::M2P(TRG,C,M);
      TRG[Index<nx,0,0>::I] += LocalSum<nx,0,0>::kernel(C,M);
    }
  };

  template<>
  struct Kernels<0,0,0> {
    static __host__ __device__ __forceinline__
    void power(fvecP&, const fvec3&) {}
    static __host__ __device__ __forceinline__
    void derivative(fvecP&, const fvec3&, const float&) {}
    static __host__ __device__ __forceinline__
    void scale(fvecP&) {}
    static __host__ __device__ __forceinline__
    void M2M(fvecP&, const fvecP&, const fvecP&) {}
    static __host__ __device__ __forceinline__
    void M2P(fvec4&, const fvecP&, const fvecP&) {}
  };


  template<int PP>
  __host__ __device__ __forceinline__
  void getCoef(fvecP &C, const fvec3 &dX, float &invR2, const float &invR) {
    C[0] = invR;
    Kernels<0,0,PP>::derivative(C,dX,invR2);
    Kernels<0,0,PP>::scale(C);
  }

  template<>
  __host__ __device__ __forceinline__
  void getCoef<1>(fvecP &C, const fvec3 &dX, float &invR2, const float &invR) {
    C[0] = invR;
    invR2 = -invR2;
    float x = dX[0], y = dX[1], z = dX[2];
    float invR3 = invR * invR2;
    C[1] = x * invR3;
    C[2] = y * invR3;
    C[3] = z * invR3;
  }

  template<>
  __host__ __device__ __forceinline__
  void getCoef<2>(fvecP &C, const fvec3 &dX, float &invR2, const float &invR) {
    getCoef<1>(C,dX,invR2,invR);
    float x = dX[0], y = dX[1], z = dX[2];
    float invR3 = invR * invR2;
    float invR5 = 3 * invR3 * invR2;
    float t = x * invR5;
    C[4] = x * t + invR3;
    C[5] = y * t;
    C[6] = z * t;
    t = y * invR5;
    C[7] = y * t + invR3;
    C[8] = z * t;
    C[9] = z * z * invR5 + invR3;
  }

  template<>
  __host__ __device__ __forceinline__
  void getCoef<3>(fvecP &C, const fvec3 &dX, float &invR2, const float &invR) {
    getCoef<2>(C,dX,invR2,invR);
    float x = dX[0], y = dX[1], z = dX[2];
    float invR3 = invR * invR2;
    float invR5 = 3 * invR3 * invR2;
    float invR7 = 5 * invR5 * invR2;
    float t = x * x * invR7;
    C[10] = x * (t + 3 * invR5);
    C[11] = y * (t +     invR5);
    C[12] = z * (t +     invR5);
    t = y * y * invR7;
    C[13] = x * (t +     invR5);
    C[16] = y * (t + 3 * invR5);
    C[17] = z * (t +     invR5);
    t = z * z * invR7;
    C[15] = x * (t +     invR5);
    C[18] = y * (t +     invR5);
    C[19] = z * (t + 3 * invR5);
    C[14] = x * y * z * invR7;
  }

  template<>
  __host__ __device__ __forceinline__
  void getCoef<4>(fvecP &C, const fvec3 &dX, float &invR2, const float &invR) {
    getCoef<3>(C,dX,invR2,invR);
    float x = dX[0], y = dX[1], z = dX[2];
    float invR3 = invR * invR2;
    float invR5 = 3 * invR3 * invR2;
    float invR7 = 5 * invR5 * invR2;
    float invR9 = 7 * invR7 * invR2;
    float t = x * x * invR9;
    C[20] = x * x * (t + 6 * invR7) + 3 * invR5;
    C[21] = x * y * (t + 3 * invR7);
    C[22] = x * z * (t + 3 * invR7);
    C[23] = y * y * (t +     invR7) + x * x * invR7 + invR5;
    C[24] = y * z * (t +     invR7);
    C[25] = z * z * (t +     invR7) + x * x * invR7 + invR5;
    t = y * y * invR9;
    C[26] = x * y * (t + 3 * invR7);
    C[27] = x * z * (t +     invR7);
    C[30] = y * y * (t + 6 * invR7) + 3 * invR5;
    C[31] = y * z * (t + 3 * invR7);
    C[32] = z * z * (t +     invR7) + y * y * invR7 + invR5;
    t = z * z * invR9;
    C[28] = x * y * (t +     invR7);
    C[29] = x * z * (t + 3 * invR7);
    C[33] = y * z * (t + 3 * invR7);
    C[34] = z * z * (t + 6 * invR7) + 3 * invR5;
  }

  template<>
  __host__ __device__ __forceinline__
  void getCoef<5>(fvecP &C, const fvec3 &dX, float &invR2, const float &invR) {
    getCoef<4>(C,dX,invR2,invR);
    float x = dX[0], y = dX[1], z = dX[2];
    float invR3 = invR * invR2;
    float invR5 = 3 * invR3 * invR2;
    float invR7 = 5 * invR5 * invR2;
    float invR9 = 7 * invR7 * invR2;
    float invR11 = 9 * invR9 * invR2;
    float t = x * x * invR11;
    C[35] = x * x * x * (t + 10 * invR9) + 15 * x * invR7;
    C[36] = x * x * y * (t +  6 * invR9) +  3 * y * invR7;
    C[37] = x * x * z * (t +  6 * invR9) +  3 * z * invR7;
    C[38] = x * y * y * (t +  3 * invR9) + x * x * x * invR9 + 3 * x * invR7;
    C[39] = x * y * z * (t +  3 * invR9);
    C[40] = x * z * z * (t +  3 * invR9) + x * x * x * invR9 + 3 * x * invR7;
    C[41] = y * y * y * (t +      invR9) + 3 * x * x * y * invR9 + 3 * y * invR7;
    C[42] = y * y * z * (t +      invR9) + x * x * z * invR9 + z * invR7;
    C[43] = y * z * z * (t +      invR9) + x * x * y * invR9 + y * invR7;
    C[44] = z * z * z * (t +      invR9) + 3 * x * x * z * invR9 + 3 * z * invR7;
    t = y * y * invR11;
    C[45] = x * y * y * (t +  6 * invR9) +  3 * x * invR7;
    C[46] = x * y * z * (t +  3 * invR9);
    C[47] = x * z * z * (t +      invR9) + x * y * y * invR9 + x * invR7;
    C[50] = y * y * y * (t + 10 * invR9) + 15 * y * invR7;
    C[51] = y * y * z * (t +  6 * invR9) + 3 * z * invR7;
    C[52] = y * z * z * (t +  3 * invR9) + y * y * y * invR9 + 3 * y * invR7;
    C[53] = z * z * z * (t +      invR9) + 3 * y * y * z * invR9 + 3 * z * invR7;
    t = z * z * invR11;
    C[48] = x * y * z * (t +  3 * invR9);
    C[49] = x * z * z * (t +  6 * invR9) +  3 * x * invR7;
    C[54] = y * z * z * (t +  6 * invR9) +  3 * y * invR7;
    C[55] = z * z * z * (t + 10 * invR9) + 15 * z * invR7;
  }

  template<>
  __host__ __device__ __forceinline__
  void getCoef<6>(fvecP &C, const fvec3 &dX, float &invR2, const float &invR) {
    getCoef<5>(C,dX,invR2,invR);
    float x = dX[0], y = dX[1], z = dX[2];
    float invR3 = invR * invR2;
    float invR5 = 3 * invR3 * invR2;
    float invR7 = 5 * invR5 * invR2;
    float invR9 = 7 * invR7 * invR2;
    float invR11 = 9 * invR9 * invR2;
    float invR13 = 11 * invR11 * invR2;
    float t = x * x * invR13;
    C[56] = x * x * x * x * (t + 15 * invR11) + 45 * x * x * invR9 + 15 * invR7;
    C[57] = x * x * x * y * (t + 10 * invR11) + 15 * x * y * invR9;
    C[58] = x * x * x * z * (t + 10 * invR11) + 15 * x * z * invR9;
    C[59] = x * x * y * y * (t +  6 * invR11) + x * x * x * x * invR11 + (6 * x * x + 3 * y * y) * invR9 + 3 * invR7;
    C[60] = x * x * y * z * (t +  6 * invR11) + 3 * y * z * invR9;
    C[61] = x * x * z * z * (t +  6 * invR11) + x * x * x * x * invR11 + (6 * x * x + 3 * z * z) * invR9 + 3 * invR7;
    C[62] = x * y * y * y * (t +  3 * invR11) + 3 * x * x * x * y * invR11 + 9 * x * y * invR9;
    C[63] = x * y * y * z * (t +  3 * invR11) + x * x * x * z * invR11 + 3 * x * z * invR9;
    C[64] = x * y * z * z * (t +  3 * invR11) + x * x * x * y * invR11 + 3 * x * y * invR9;
    C[65] = x * z * z * z * (t +  3 * invR11) + 3 * x * x * x * z * invR11 + 9 * x * z * invR9;
    C[66] = y * y * y * y * (t +      invR11) + 6 * x * x * y * y * invR11 + (3 * x * x + 6 * y * y) * invR9 + 3 * invR7;
    C[67] = y * y * y * z * (t +      invR11) + 3 * x * x * y * z * invR11 + 3 * y * z * invR9;
    C[68] = y * y * z * z * (t +      invR11) + (x * x * y * y + x * x * z * z) * invR11 + (x * x + y * y + z * z) * invR9 + invR7;
    C[69] = y * z * z * z * (t +      invR11) + 3 * x * x * y * z * invR11 + 3 * y * z * invR9;
    C[70] = z * z * z * z * (t +      invR11) + 6 * x * x * z * z * invR11 + (3 * x * x + 6 * z * z) * invR9 + 3 * invR7;
    t = y * y * invR13;
    C[71] = x * y * y * y * (t + 10 * invR11) + 15 * x * y * invR9;
    C[72] = x * y * y * z * (t +  6 * invR11) + 3 * x * z * invR9;
    C[73] = x * y * z * z * (t +  3 * invR11) + x * y * y * y * invR11 + 3 * x * y * invR9;
    C[74] = x * z * z * z * (t +      invR11) + 3 * x * y * y * z * invR11 + 3 * x * z * invR9;
    C[77] = y * y * y * y * (t + 15 * invR11) + 45 * y * y * invR9 + 15 * invR7;
    C[78] = y * y * y * z * (t + 10 * invR11) + 15 * y * z * invR9;
    C[79] = y * y * z * z * (t +  6 * invR11) + y * y * y * y * invR11 + (6 * y * y + 3 * z * z) * invR9 + 3 * invR7;
    C[80] = y * z * z * z * (t +  3 * invR11) + 3 * y * y * y * z * invR11 + 9 * y * z * invR9;
    C[81] = z * z * z * z * (t +      invR11) + 6 * y * y * z * z * invR11 + (3 * y * y + 6 * z * z) * invR9 + 3 * invR7;
    t = z * z * invR13;
    C[75] = x * y * z * z * (t +  6 * invR11) + 3 * x * y * invR9;
    C[76] = x * z * z * z * (t + 10 * invR11) + 15 * x * z * invR9;
    C[82] = y * z * z * z * (t + 10 * invR11) + 15 * y * z * invR9;
    C[83] = z * z * z * z * (t + 15 * invR11) + 45 * z * z * invR9 + 15 * invR7;
  }


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
      Kernels<0,0,P-1>::power(M,dX);
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
      fvecP C;
      C[0] = 1;
      Kernels<0,0,P-1>::power(C,dX);
      for (int j=0; j<NTERM; j++) Mi[j] += C[j] * Mj[0];
      Kernels<0,0,P-1>::M2M(Mi,C,Mj);
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
	    const fvecP & __restrict__ M,
	    float EPS2) {
    const fvec3 dX = pos_i - pos_j;
    const float R2 = norm(dX) + EPS2;
    const float invR = rsqrtf(R2);
#if 1
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
    const float invR1 = M[0] * invR;
    float invR2 = invR * invR;
    fvecP C;
    getCoef<P-1>(C,dX,invR2,invR1);
    acc[0] -= C[0];
    for (int i=1; i<NTERM; i++) acc[0] -= M[i] * C[i];
    for (int i=1; i<4; i++) acc[i] += C[i];
    Kernels<0,0,1>::M2P(acc,C,M);
#endif
    return acc;
  }
}
