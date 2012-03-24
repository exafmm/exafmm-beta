/*
 * Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
#define KERNEL
#include "kernel.h"
#undef KERNEL

namespace
{
    //! Get r,theta,phi from x,y,z
    void cart2sph(real& r, real& theta, real& phi, vect dist)
    {
        r = sqrt(norm(dist)) + EPS;                                 // r = sqrt(x^2 + y^2 + z^2) + eps
        theta = acos(dist[2] / r);                                  // theta = acos(z / r)
        if (fabs(dist[0]) + fabs(dist[1]) < EPS)                    // If |x| < eps & |y| < eps
        {
            phi = 0;                                                  //  phi can be anything so we set it to 0
        }
        else if (fabs(dist[0]) < EPS)                             // If |x| < eps
        {
            phi = dist[1] / fabs(dist[1]) * M_PI * 0.5;               //  phi = sign(y) * pi / 2
        }
        else if (dist[0] > 0)                                     // If x > 0
        {
            phi = atan(dist[1] / dist[0]);                            //  phi = atan(y / x)
        }
        else                                                      // If x < 0
        {
            phi = atan(dist[1] / dist[0]) + M_PI;                     //  phi = atan(y / x) + pi
        }                                                           // End if for x,y cases
    }

    //! Spherical to cartesian coordinates
    template<typename T>
    void sph2cart(real r, real theta, real phi, const T &spherical, T &cartesian)
    {
        cartesian[0] = sin(theta) * cos(phi) * spherical[0]         // x component (not x itself)
                       + cos(theta) * cos(phi) / r * spherical[1]
                       - sin(phi) / r / sin(theta) * spherical[2];
        cartesian[1] = sin(theta) * sin(phi) * spherical[0]         // y component (not y itself)
                       + cos(theta) * sin(phi) / r * spherical[1]
                       + cos(phi) / r / sin(theta) * spherical[2];
        cartesian[2] = cos(theta) * spherical[0]                    // z component (not z itself)
                       - sin(theta) / r * spherical[1];
    }

}


template<int n, int m>
struct Index
{
    static const int npm = Index < n, m - 1 >::npm + 1;
    static const int nmm = Index < n, m - 1 >::nmm - 1;
    static const int nms = Index < n, m - 1 >::nms + 1;
};

template<int n>
struct Index<n, 0>
{
    static const int npm = Index < n - 1, n - 1 >::npm + n + 1;
    static const int nmm = Index < n - 1, n - 1 >::npm + n + 1;
    static const int nms = Index < n - 1, n - 1 >::nms + 1;
};

template<>
struct Index<0, 0>
{
    static const int npm = 0;
    static const int nmm = 0;
    static const int nms = 0;
};

template < int p, int n = p - 1, int m = p - 1 >
struct Expansion
{
    static inline void getYnm(real &x, real &y, real &rho,
                              real &rhom, real &rhon, real &beta, complex &eim,
                              real &Pn, real &P0, real &P1, real &P2,
                              real *prefactor, complex *Ynm)
    {
        Expansion < p, n - 1, m >::getYnm(x, y, rho, rhom, rhon, beta, eim, Pn, P0, P1, P2, prefactor, Ynm);
        Ynm[Index<n,m>::npm] = rhon * P0 * prefactor[Index<n,m>::npm] * eim;
        Ynm[Index<n,m>::nmm] = std::conj(Ynm[Index<n,m>::npm]);
        P2 = P1;
        P1 = P0;
        P0 = (x * (2 * n + 1) * P1 - (n + m) * P2) / (n - m + 1);
        rhon *= rho;
    }
    static inline void getYnmTheta(real &x, real &y, real &rho,
                                   real &rhom, real &rhon, real &beta, complex &eim,
                                   real &Pn, real &P0, real &P1, real &P2,
                                   real *prefactor, complex *Ynm, complex *YnmTheta)
    {
        Expansion < p, n - 1, m >::getYnmTheta(x, y, rho, rhom, rhon, beta, eim, Pn, P0, P1, P2, prefactor, Ynm, YnmTheta);
        Ynm[Index<n,m>::npm] = rhon * P0 * prefactor[Index<n,m>::npm] * eim;
        Ynm[Index<n,m>::nmm] = std::conj(Ynm[Index<n,m>::npm]);
        P2 = P1;
        P1 = P0;
        P0 = (x * (2 * n + 1) * P1 - (n + m) * P2) / (n - m + 1);
        YnmTheta[Index<n,m>::npm] = rhon * ((n - m + 1) * P0 - (n + 1) * x * P1) / y * prefactor[Index<n,m>::npm] * eim;
        rhon *= rho;
    }
};

template<int p, int m>
struct Expansion<p, m, m>
{
    static inline void getYnm(real &x, real &y, real &rho,
                              real &rhom, real &rhon, real &beta, complex &eim,
                              real &Pn, real &P0, real &P1, real &P2,
                              real *prefactor, complex *Ynm)
    {
        Expansion < p, p - 1, m - 1 >::getYnm(x, y, rho, rhom, rhon, beta, eim, Pn, P0, P1, P2, prefactor, Ynm);
        const complex I(0., 1.);
        eim = std::exp(I * real(m * beta));
        Pn = -Pn * (2 * m - 1) * y;
        P0 = Pn;
        Ynm[Index<m,m>::npm] = rhom * P0 * prefactor[Index<m,m>::npm] * eim;
        Ynm[Index<m,m>::nmm] = std::conj(Ynm[Index<m,m>::npm]);
        P1 = P0;
        P0 = x * (2 * m + 1) * P0;
        rhom *= rho;
        rhon = rhom;
    }
    static inline void getYnmTheta(real &x, real &y, real &rho,
                                   real &rhom, real &rhon, real &beta, complex &eim,
                                   real &Pn, real &P0, real &P1, real &P2,
                                   real *prefactor, complex *Ynm, complex *YnmTheta)
    {
        Expansion < p, p - 1, m - 1 >::getYnmTheta(x, y, rho, rhom, rhon, beta, eim, Pn, P0, P1, P2, prefactor, Ynm, YnmTheta);
        const complex I(0., 1.);
        eim = std::exp(I * real(m * beta));
        Pn = -Pn * (2 * m - 1) * y;
        P0 = Pn;
        Ynm[Index<m,m>::npm] = rhom * P0 * prefactor[Index<m,m>::npm] * eim;
        Ynm[Index<m,m>::nmm] = std::conj(Ynm[Index<m,m>::npm]);
        P1 = P0;
        P0 = x * (2 * m + 1) * P0;
        YnmTheta[Index<m,m>::npm] = rhom * (P0 - (m + 1) * x * P1) / y * prefactor[Index<m,m>::npm] * eim;
        rhom *= rho;
        rhon = rhom;
    }
};

template<int p>
struct Expansion<p, 0, 0>
{
    static inline void getYnm(real &, real &, real &rho,
                              real &rhom, real &rhon, real&, complex&,
                              real&, real &, real &, real&,
                              real*, complex *Ynm)
    {
        Ynm[0] = rhom;
        rhom *= rho;
        rhon = rhom;
    }
    static inline void getYnmTheta(real &, real &, real &rho,
                                   real &rhom, real &rhon, real&, complex&,
                                   real&, real &, real &, real&,
                                   real*, complex *Ynm, complex *YnmTheta)
    {
        Ynm[0] = rhom;
        YnmTheta[0] = 0;
        rhom *= rho;
        rhon = rhom;
    }
};

template<int n, int m>
struct Terms
{
    static inline void P2M(Mset &M, const real C[4], complex *Ynm)
    {
        Terms < n, m - 1 >::P2M(M, C, Ynm);
        M[Index<n,m>::nms] += C[0] * Ynm[Index<n,m>::npm];
        M[Index<n,m>::nms+NTERM] += C[1] * Ynm[Index<n,m>::npm];
        M[Index<n,m>::nms+NTERM+NTERM] += C[2] * Ynm[Index<n,m>::npm];
        M[Index<n,m>::nms+NTERM+NTERM+NTERM] += C[3] * Ynm[Index<n,m>::npm];
    }
};

template<int n>
struct Terms<n, 0>
{
    static inline void P2M(Mset &M, const real C[4], complex *Ynm)
    {
        Terms < n - 1, n - 1 >::P2M(M, C, Ynm);
        M[Index<n,0>::nms] += C[0] * Ynm[Index<n,0>::npm];
        M[Index<n,0>::nms+NTERM] += C[1] * Ynm[Index<n,0>::npm];
        M[Index<n,0>::nms+NTERM+NTERM] += C[2] * Ynm[Index<n,0>::npm];
        M[Index<n,0>::nms+NTERM+NTERM+NTERM] += C[3] * Ynm[Index<n,0>::npm];
    }
};

template<>
struct Terms<0, 0>
{
    static inline void P2M(Mset &M, const real C[4], complex *Ynm)
    {
        M[Index<0,0>::nms] += C[0] * Ynm[Index<0,0>::npm];
        M[Index<0,0>::nms+NTERM] += C[1] * Ynm[Index<0,0>::npm];
        M[Index<0,0>::nms+NTERM+NTERM] += C[2] * Ynm[Index<0,0>::npm];
        M[Index<0,0>::nms+NTERM+NTERM+NTERM] += C[3] * Ynm[Index<0,0>::npm];
    }
};

template < int p, int j = p - 1, int k = p - 1, int n = p - 1, int m = p - 1 >
struct M2Ltemplate
{
    static inline void loop(Mset &M, Lset &L, complex *Cnm, complex *Ynm)
    {
        M2Ltemplate < p, j, k, n, m - 1 >::loop(M, L, Cnm, Ynm);
        int nm   = n * n + n + m;
        int nms  = n * (n + 1) / 2 + m;
        int jk = j * j + j + k;
        int jks = j * (j + 1) / 2 + k;
        int jknm = jk * P * P + nm;
        int jnkm = (j + n) * (j + n) + j + n + m - k;
        L[jks] += M[nms] * Cnm[jknm] * Ynm[jnkm];
        nm   = n * n + n - m;
        jknm = jk * P * P + nm;
        jnkm = (j + n) * (j + n) + j + n - m - k;
        L[jks] += std::conj(M[nms]) * Cnm[jknm] * Ynm[jnkm];
    }
};

template<int p, int j, int k, int n>
struct M2Ltemplate<p, j, k, n, 1>
{
    static inline void loop(Mset &M, Lset &L, complex *Cnm, complex *Ynm)
    {
        M2Ltemplate < p, j, k, n - 1, n - 1 >::loop(M, L, Cnm, Ynm);
        int nm   = n * n + n;
        int nms  = n * (n + 1) / 2;
        int jk = j * j + j + k;
        int jks = j * (j + 1) / 2 + k;
        int jknm = jk * P * P + nm;
        int jnkm = (j + n) * (j + n) + j + n - k;
        L[jks] += M[nms] * Cnm[jknm] * Ynm[jnkm];
        nm   = n * n + n + 1;
        nms  = n * (n + 1) / 2 + 1;
        jknm = jk * P * P + nm;
        jnkm = (j + n) * (j + n) + j + n + 1 - k;
        L[jks] += M[nms] * Cnm[jknm] * Ynm[jnkm];
        nm   = n * n + n - 1;
        jknm = jk * P * P + nm;
        jnkm = (j + n) * (j + n) + j + n - 1 - k;
        L[jks] += std::conj(M[nms]) * Cnm[jknm] * Ynm[jnkm];
    }
};

template<int p, int j, int k>
struct M2Ltemplate<p, j, k, 2, 1>
{
    static inline void loop(Mset &M, Lset &L, complex *Cnm, complex *Ynm)
    {
        M2Ltemplate < p, j, k - 1, p - 1, p - 1 >::loop(M, L, Cnm, Ynm);
        int jk = j * j + j + k;
        int jks = j * (j + 1) / 2 + k;
        L[jks] += M[0] * Cnm[jk*P*P] * Ynm[j*j+j-k];
        int jknm = jk * P * P + 6;
        int jnkm = (j + 2) * (j + 2) + j + 2 - k;
        L[jks] += M[3] * Cnm[jknm] * Ynm[jnkm];
        jknm = jk * P * P + 7;
        jnkm = (j + 2) * (j + 2) + j + 3 - k;
        L[jks] += M[4] * Cnm[jknm] * Ynm[jnkm];
        jknm = jk * P * P + 5;
        jnkm = (j + 2) * (j + 2) + j + 1 - k;
        L[jks] += std::conj(M[4]) * Cnm[jknm] * Ynm[jnkm];
    }
};

template<int p, int j>
struct M2Ltemplate<p, j, 0, 2, 1>
{
    static inline void loop(Mset &M, Lset &L, complex *Cnm, complex *Ynm)
    {
        M2Ltemplate < p, j - 1, j - 1, p - 1, p - 1 >::loop(M, L, Cnm, Ynm);
        int jk = j * j + j;
        int jks = j * (j + 1) / 2;
        L[jks] += M[0] * Cnm[jk*P*P] * Ynm[j*j+j];
        int jknm = jk * P * P + 6;
        int jnkm = (j + 2) * (j + 2) + j + 2;
        L[jks] += M[3] * Cnm[jknm] * Ynm[jnkm];
        jknm = jk * P * P + 7;
        jnkm = (j + 2) * (j + 2) + j + 3;
        L[jks] += M[4] * Cnm[jknm] * Ynm[jnkm];
        jknm = jk * P * P + 5;
        jnkm = (j + 2) * (j + 2) + j + 1;
        L[jks] += std::conj(M[4]) * Cnm[jknm] * Ynm[jnkm];
    }
};

template<int p>
struct M2Ltemplate<p, 0, 0, 2, 1>
{
    static inline void loop(Mset &M, Lset &L, complex *Cnm, complex *Ynm)
    {
        L[0] += M[0] * Cnm[0] * Ynm[0];
        L[0] += M[3] * Cnm[6] * Ynm[6];
        L[0] += M[4] * Cnm[7] * Ynm[7];
        L[0] += std::conj(M[4]) * Cnm[5] * Ynm[5];
    }
};


void Kernel<Stokes>::evalMultipole(real rho, real alpha, real beta, complex *Ynm) const
{
    real x = std::cos(alpha);
    real y = std::sin(alpha);
    real rhom = 1, rhon = rhom, P0 = x, P1 = 1, P2 = 1, Pn = 1;
    complex eim = 1;
    Expansion<P>::getYnm(x, y, rho, rhom, rhon, beta, eim, Pn, P0, P1, P2, prefactor, Ynm);
}

void Kernel<Stokes>::evalMultipoleTheta(real rho, real alpha, real beta, complex *Ynm, complex *YnmTheta) const
{
    real x = std::cos(alpha);
    real y = std::sin(alpha);
    real rhom = 1, rhon = rhom, P0 = x, P1 = 1, P2 = 1, Pn = 1;
    complex eim = 1;
    Expansion<P>::getYnmTheta(x, y, rho, rhom, rhon, beta, eim, Pn, P0, P1, P2, prefactor, Ynm, YnmTheta);
}

void Kernel<Stokes>::evalLocal(real rho, real alpha, real beta, complex *Ynm) const
{
    real x = std::cos(alpha);
    real y = std::sin(alpha);
    real invR = 1 / rho;
    real rhom = invR, rhon = rhom, P0 = x, P1 = 1, P2 = 1, Pn = 1;
    complex eim = 1;
    Expansion<2*P>::getYnm(x, y, invR, rhom, rhon, beta, eim, Pn, P0, P1, P2, prefactor, Ynm);
}

void Kernel<Stokes>::evalLocalTheta(real rho, real alpha, real beta, complex *Ynm, complex *YnmTheta) const
{
    real x = std::cos(alpha);
    real y = std::sin(alpha);
    real invR = 1 / rho;
    real rhom = invR, rhon = rhom, P0 = x, P1 = 1, P2 = 1, Pn = 1;
    complex eim = 1;
    Expansion<2*P>::getYnmTheta(x, y, invR, rhom, rhon, beta, eim, Pn, P0, P1, P2, prefactor, Ynm, YnmTheta);
}

#if STOKES
void Kernel<Stokes>::P2M(C_iter Ci) const
{
    real Rmax = 0;
    complex Ynm[4*P*P];
    for (B_iter B = Ci->LEAF; B != Ci->LEAF + Ci->NDLEAF; ++B)
    {
        vect dist = B->X - Ci->X;
        real R = std::sqrt(norm(dist));
        if (R > Rmax) Rmax = R;
        real rho, alpha, beta;
        cart2sph(rho, alpha, beta, dist);
        evalMultipole(rho, alpha, -beta, Ynm);
        real SRC[4] = {B->FORCE[0],
                       B->FORCE[1],
                       B->FORCE[2],
                       B->X[0] * B->FORCE[0] + B->X[1] * B->FORCE[1] + B->X[2] * B->FORCE[2]
                      };
        Terms < P - 1, P - 1 >::P2M(Ci->M, SRC, Ynm);
    }
    Ci->RMAX = Rmax;
    Ci->RCRIT = std::min(Ci->R, Rmax);
}
void Kernel<Stokes>::D2M(C_iter Ci) const
{
    const complex I(0., 1.);
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    for (B_iter B = Ci->LEAF; B != Ci->LEAF + Ci->NDLEAF; ++B)
    {
        vect dist = B->X - Ci->X;
        real rho, alpha, beta, factor;
        vect d = 0, gradient = 0;
        cart2sph(rho, alpha, beta, dist);
        evalMultipoleTheta(rho, alpha, -beta, Ynm, YnmTheta);
        d[0] = B->FORCE[0] * B->X[2];
        d[1] = B->FORCE[1] * B->X[2];
        d[2] = -B->FORCE[2] * B->X[2];
        real u = -B->FORCE[2];
        for (int n = 0; n != P; ++n)
        {
            int nm  = n * n + n;
            int nms = n * (n + 1) / 2;

            factor = 1.0 / rho * n;
            gradient[0] += std::real(Ynm[nm]) * factor;
            gradient[1] += std::real(YnmTheta[nm]);
            Ci->M[nms] += u * (gradient[0] * d[0] + gradient[1] * d[1]);

            for (int m = 1; m <= n; ++m)
            {
                const int nm  = n * (n + 1) + m;
                const int nms = n * (n + 1) / 2 + m;
                gradient[0] += 2 * std::real(Ynm[nm]) * factor;
                gradient[1] += 2 * std::real(YnmTheta[nm]);
                gradient[2] += 2 * std::real(Ynm[nm] * I) * m;
                Ci->M[nms] += u * (gradient[0] * d[0] + gradient[1] * d[1] + gradient[2] * d[2]);
            }
        }
    }
}
#endif
void Kernel<Stokes>::M2M(C_iter Ci) const
{
    real Rmax = Ci->RMAX;
    const complex I(0., 1.);                                      // Imaginary unit
    complex Ynm[4*P*P];
    for (C_iter Cj = Cj0 + Ci->CHILD; Cj != Cj0 + Ci->CHILD + Ci->NCHILD; ++Cj)
    {
        vect dist = Ci->X - Cj->X;
        real R = std::sqrt(norm(dist)) + Cj->RCRIT;
        if( R > Rmax ) Rmax = R;
        real rho, alpha, beta;
        cart2sph(rho, alpha, beta, dist);
        evalMultipole(rho, alpha, -beta, Ynm);
        int offsets[3] = {NTERM, 2*NTERM, 3*NTERM};
        for (int j = 0; j != P; ++j) {
            for (int k = 0; k <= j; ++k) {
                const int jk = j * j + j + k;
                const int jks = j * (j + 1) / 2 + k;
                complex M[4] = {0};
                for (int n = 0; n <= j; ++n) {
                    for (int m = -n; m <= std::min(k - 1, n); ++m) {
                        if (j - n >= k - m) {
                            const int jnkm  = (j - n) * (j - n) + j - n + k - m;
                            const int jnkms = (j - n) * (j - n + 1) / 2 + k - m;
                            const int nm    = n * n + n + m;
                            complex factor = std::pow(I, real(m - abs(m))) * Ynm[nm] * real(ODDEVEN(n) * Anm[nm] * Anm[jnkm] / Anm[jk]);
                            M[0] += Cj->M[jnkms] * factor;
                            M[1] += Cj->M[jnkms+offsets[0]] * factor;
                            M[2] += Cj->M[jnkms+offsets[1]] * factor;
                            M[3] += Cj->M[jnkms+offsets[2]] * factor;
                        }
                    }
                    for (int m = k; m <= n; ++m)
                    {
                        if (j - n >= m - k)
                        {
                            const int jnkm  = (j - n) * (j - n) + j - n + k - m;
                            const int jnkms = (j - n) * (j - n + 1) / 2 - k + m;
                            const int nm    = n * n + n + m;
                            complex factor = Ynm[nm] * real(ODDEVEN(k + n + m) * Anm[nm] * Anm[jnkm] / Anm[jk]);
                            M[0] += std::conj(Cj->M[jnkms]) * factor;
                            M[1] += Cj->M[jnkms+offsets[0]] * factor;
                            M[2] += Cj->M[jnkms+offsets[1]] * factor;
                            M[3] += Cj->M[jnkms+offsets[2]] * factor;
                        }
                    }
                }
                Ci->M[jks] += M[0] * EPS;
                Ci->M[jks+offsets[0]] += M[1] * EPS;
                Ci->M[jks+offsets[1]] += M[2] * EPS;
                Ci->M[jks+offsets[2]] += M[3] * EPS;
            }
        }
    }
    Ci->RCRIT = std::min(Ci->R,Rmax);
}

void Kernel<Stokes>::M2L(C_iter Ci, C_iter Cj) const
{
    complex Ynm[4*P*P];
    vect dist = Ci->X - Cj->X - Xperiodic;
    real rho, alpha, beta;
    cart2sph(rho, alpha, beta, dist);
    evalLocal(rho, alpha, beta, Ynm);
    int offsets[3] = {NTERM, 2*NTERM, 3*NTERM};
    for (int j = 0; j != P; ++j)
    {
        for (int k = 0; k <= j; ++k)
        {
            const int jk = j * j + j + k;
            const int jks = j * (j + 1) / 2 + k;
            complex L[4] = {0};
            for (int n = 0; n != P; ++n)
            {
                for (int m = -n; m < 0; ++m)
                {
                    const int nm   = n * n + n + m;
                    const int nms  = n * (n + 1) / 2 - m;
                    const int jknm = jk * P2 + nm;
                    const int jnkm = (j + n) * (j + n) + j + n + m - k;
                    L[0] += std::conj(Cj->M[nms]) * Cnm[jknm] * Ynm[jnkm];
                    L[1] += std::conj(Cj->M[nms+offsets[0]]) * Cnm[jknm] * Ynm[jnkm];
                    L[2] += std::conj(Cj->M[nms+offsets[1]]) * Cnm[jknm] * Ynm[jnkm];
                    L[3] += std::conj(Cj->M[nms+offsets[2]]) * Cnm[jknm] * Ynm[jnkm];
                }
                for (int m = 0; m <= n; ++m)
                {
                    const int nm   = n * n + n + m;
                    const int nms  = n * (n + 1) / 2 + m;
                    const int jknm = jk * P2 + nm;
                    const int jnkm = (j + n) * (j + n) + j + n + m - k;
                    L[0] += Cj->M[nms] * Cnm[jknm] * Ynm[jnkm];
                    L[1] += Cj->M[nms+offsets[0]] * Cnm[jknm] * Ynm[jnkm];
                    L[2] += Cj->M[nms+offsets[1]] * Cnm[jknm] * Ynm[jnkm];
                    L[3] += Cj->M[nms+offsets[2]] * Cnm[jknm] * Ynm[jnkm];
                }
            }
            Ci->L[jks] += L[0];
            Ci->L[jks+offsets[0]] += L[1];
            Ci->L[jks+offsets[1]] += L[2];
            Ci->L[jks+offsets[2]] += L[3];
        }
    }
}

// Only for tree code and hybrid tree method
void Kernel<Stokes>::M2P(C_iter Ci, C_iter Cj) const
{
    const complex I(0., 1.);                                      // Imaginary unit
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    int offsets[3] = {NTERM, 2*NTERM, 3*NTERM};
    for (B_iter B = Ci->LEAF; B != Ci->LEAF + Ci->NDLEAF; ++B)
    {
        vect dist = B->X - Cj->X - Xperiodic;
        real r, theta, phi, factor;
        vect gradient[4] = {0, 0, 0, 0};
        vect cartesian = 0;
        cart2sph(r, theta, phi, dist);
        evalLocalTheta(r, theta, phi, Ynm, YnmTheta);
        int offsets[3] = {NTERM, 2*NTERM, 3*NTERM};
        for (int n = 0; n != P; ++n)
        {
            int nm  = n * n + n;
            int nms = n * (n + 1) / 2;
            B->TRG[0] += std::real(Cj->M[nms] * Ynm[nm]);
            B->TRG[1] += std::real(Cj->M[nms+offsets[0]] * Ynm[nm]);
            B->TRG[2] += std::real(Cj->M[nms+offsets[1]] * Ynm[nm]);

            factor = 1.0 / r * (n + 1);
            gradient[0][0] -= std::real(Cj->M[nms] * Ynm[nm]) * factor;
            gradient[0][1] += std::real(Cj->M[nms] * YnmTheta[nm]);

            gradient[1][0] -= std::real(Cj->M[nms+offsets[0]] * Ynm[nm]) * factor;
            gradient[1][1] += std::real(Cj->M[nms+offsets[0]] * YnmTheta[nm]);

            gradient[2][0] -= std::real(Cj->M[nms+offsets[1]] * Ynm[nm]) * factor;
            gradient[2][1] += std::real(Cj->M[nms+offsets[1]] * YnmTheta[nm]);

            gradient[3][0] -= std::real(Cj->M[nms+offsets[2]] * Ynm[nm]) * factor;
            gradient[3][1] += std::real(Cj->M[nms+offsets[2]] * YnmTheta[nm]);
            for (int m = 1; m <= n; ++m)
            {
                nm  = n * n + n + m;
                nms = n * (n + 1) / 2 + m;
                B->TRG[0] += 2 * std::real(Cj->M[nms] * Ynm[nm]);
                B->TRG[1] += 2 * std::real(Cj->M[nms+offsets[0]] * Ynm[nm]);
                B->TRG[2] += 2 * std::real(Cj->M[nms+offsets[1]] * Ynm[nm]);

                gradient[0][0] -= 2 * std::real(Cj->M[nms] * Ynm[nm]) * factor;
                gradient[0][1] += 2 * std::real(Cj->M[nms] * YnmTheta[nm]);
                gradient[0][2] += 2 * std::real(Cj->M[nms] * Ynm[nm] * I) * m;

                gradient[1][0] -= 2 * std::real(Cj->M[nms+offsets[0]] * Ynm[nm]) * factor;
                gradient[1][1] += 2 * std::real(Cj->M[nms+offsets[0]] * YnmTheta[nm]);
                gradient[1][2] += 2 * std::real(Cj->M[nms+offsets[0]] * Ynm[nm] * I) * m;

                gradient[2][0] -= 2 * std::real(Cj->M[nms+offsets[1]] * Ynm[nm]) * factor;
                gradient[2][1] += 2 * std::real(Cj->M[nms+offsets[1]] * YnmTheta[nm]);
                gradient[2][2] += 2 * std::real(Cj->M[nms+offsets[1]] * Ynm[nm] * I) * m;

                gradient[3][0] -= 2 * std::real(Cj->M[nms+offsets[2]] * Ynm[nm]) * factor;
                gradient[3][1] += 2 * std::real(Cj->M[nms+offsets[2]] * YnmTheta[nm]);
                gradient[3][2] += 2 * std::real(Cj->M[nms+offsets[2]] * Ynm[nm] * I) * m;
            }
        }
        sph2cart(r, theta, phi, gradient[0], cartesian);
        cartesian *= -B->X[0];
        gradient[0] = cartesian;
        sph2cart(r, theta, phi, gradient[1], cartesian);
        cartesian *= -B->X[1];
        gradient[1] = cartesian;
        sph2cart(r, theta, phi, gradient[2], cartesian);
        cartesian *= -B->X[2];
        gradient[2] = cartesian;
        sph2cart(r, theta, phi, gradient[3], cartesian);
        gradient[3] = cartesian;

        B->TRG[0] += (gradient[0][0] + gradient[1][0] + gradient[2][0] + gradient[3][0]);
        B->TRG[1] += (gradient[0][1] + gradient[1][1] + gradient[2][1] + gradient[3][1]);
        B->TRG[2] += (gradient[0][2] + gradient[1][2] + gradient[2][2] + gradient[3][2]);
    }
}

void Kernel<Stokes>::L2L(C_iter Ci) const
{
    const complex I(0., 1.);                                      // Imaginary unit
    complex Ynm[4*P*P], factor;
    C_iter Cj = Ci0 + Ci->PARENT;
    vect dist = Ci->X - Cj->X;
    real rho, alpha, beta;
    cart2sph(rho, alpha, beta, dist);
    evalMultipole(rho, alpha, beta, Ynm);
    int offsets[3] = {NTERM, 2*NTERM, 3*NTERM};
    for (int j = 0; j != P; ++j)
    {
        for (int k = 0; k <= j; ++k)
        {
            const int jk = j * j + j + k;
            const int jks = j * (j + 1) / 2 + k;
            complex L[4] = {0};
            for (int n = j; n != P; ++n)
            {
                for (int m = j + k - n; m < 0; ++m)
                {
                    const int jnkm = (n - j) * (n - j) + n - j + m - k;
                    const int nm   = n * n + n - m;
                    const int nms  = n * (n + 1) / 2 - m;
                    factor = Ynm[jnkm] * real(ODDEVEN(k) * Anm[jnkm] * Anm[jk] / Anm[nm]);
                    L[0] += std::conj(Cj->L[nms]) * factor;
                    L[1] += std::conj(Cj->L[nms+offsets[0]]) * factor;
                    L[2] += std::conj(Cj->L[nms+offsets[1]]) * factor;
                    L[3] += std::conj(Cj->L[nms+offsets[2]]) * factor;
                }
                for (int m = 0; m <= n; ++m)
                {
                    if (n - j >= abs(m - k))
                    {
                        const int jnkm = (n - j) * (n - j) + n - j + m - k;
                        const int nm   = n * n + n + m;
                        const int nms  = n * (n + 1) / 2 + m;
                        factor = std::pow(I, real(m - k - abs(m - k))) * Ynm[jnkm] * Anm[jnkm] * Anm[jk] / Anm[nm];
                        L[0] += Cj->L[nms] * factor;
                        L[1] += Cj->L[nms+offsets[0]] * factor;
                        L[2] += Cj->L[nms+offsets[1]] * factor;
                        L[3] += Cj->L[nms+offsets[2]] * factor;
                    }
                }
            }
            Ci->L[jks] += L[0] * EPS;
            Ci->L[jks+offsets[0]] += L[1] * EPS;
            Ci->L[jks+offsets[1]] += L[2] * EPS;
            Ci->L[jks+offsets[2]] += L[3] * EPS;
        }
    }
}

void Kernel<Stokes>::L2P(C_iter Ci) const
{
    const complex I(0., 1.);                                      // Imaginary unit
    complex Ynm[4*P*P], YnmTheta[4*P*P];
    int offsets[3] = {NTERM, 2*NTERM, 3*NTERM};
    for (B_iter B = Ci->LEAF; B != Ci->LEAF + Ci->NDLEAF; ++B)
    {
        vect dist = B->X - Ci->X;
        vect gradient[4] = {0, 0, 0, 0};
        vect cartesian = 0;
        real r, theta, phi, factor;
        cart2sph(r, theta, phi, dist);
        evalMultipoleTheta(r, theta, phi, Ynm, YnmTheta);
        for (int n = 0; n != P; ++n)
        {
            int nm  = n * n + n;
            int nms = n * (n + 1) / 2;

            B->TRG[0] += std::real(Ci->L[nms] * Ynm[nm]);
            B->TRG[1] += std::real(Ci->L[nms+offsets[0]] * Ynm[nm]);
            B->TRG[2] += std::real(Ci->L[nms+offsets[1]] * Ynm[nm]);

            factor = 1.0 / r * n;
            gradient[0][0] += std::real(Ci->L[nms] * Ynm[nm]) * factor;
            gradient[0][1] += std::real(Ci->L[nms] * YnmTheta[nm]);

            gradient[1][0] += std::real(Ci->L[nms+offsets[0]] * Ynm[nm]) * factor;
            gradient[1][1] += std::real(Ci->L[nms+offsets[0]] * YnmTheta[nm]);

            gradient[2][0] += std::real(Ci->L[nms+offsets[1]] * Ynm[nm]) * factor;
            gradient[2][1] += std::real(Ci->L[nms+offsets[1]] * YnmTheta[nm]);

            gradient[3][0] += std::real(Ci->L[nms+offsets[2]] * Ynm[nm]) * factor;
            gradient[3][1] += std::real(Ci->L[nms+offsets[2]] * YnmTheta[nm]);
            for (int m = 1; m <= n; ++m)
            {
                nm  = n * n + n + m;
                nms = n * (n + 1) / 2 + m;
                B->TRG[0] += 2 * std::real(Ci->L[nms] * Ynm[nm]);
                B->TRG[1] += 2 * std::real(Ci->L[nms+offsets[0]] * Ynm[nm]);
                B->TRG[2] += 2 * std::real(Ci->L[nms+offsets[1]] * Ynm[nm]);

                gradient[0][0] += 2 * std::real(Ci->L[nms] * Ynm[nm]) * factor;
                gradient[0][1] += 2 * std::real(Ci->L[nms] * YnmTheta[nm]);
                gradient[0][2] += 2 * std::real(Ci->L[nms] * Ynm[nm] * I) * m;

                gradient[1][0] += 2 * std::real(Ci->L[nms+offsets[0]] * Ynm[nm]) * factor;
                gradient[1][1] += 2 * std::real(Ci->L[nms+offsets[0]] * YnmTheta[nm]);
                gradient[1][2] += 2 * std::real(Ci->L[nms+offsets[0]] * Ynm[nm] * I) * m;

                gradient[2][0] += 2 * std::real(Ci->L[nms+offsets[1]] * Ynm[nm]) * factor;
                gradient[2][1] += 2 * std::real(Ci->L[nms+offsets[1]] * YnmTheta[nm]);
                gradient[2][2] += 2 * std::real(Ci->L[nms+offsets[1]] * Ynm[nm] * I) * m;

                gradient[3][0] += 2 * std::real(Ci->L[nms+offsets[2]] * Ynm[nm]) * factor;
                gradient[3][1] += 2 * std::real(Ci->L[nms+offsets[2]] * YnmTheta[nm]);
                gradient[3][2] += 2 * std::real(Ci->L[nms+offsets[2]] * Ynm[nm] * I) * m;
            }
        }
        sph2cart(r, theta, phi, gradient[0], cartesian);
        cartesian *= -B->X[0];
        gradient[0] = cartesian;
        sph2cart(r, theta, phi, gradient[1], cartesian);
        cartesian *= -B->X[1];
        gradient[1] = cartesian;
        sph2cart(r, theta, phi, gradient[2], cartesian);
        cartesian *= -B->X[2];
        gradient[2] = cartesian;
        sph2cart(r, theta, phi, gradient[3], cartesian);
        gradient[3] = cartesian;

        B->TRG[0] += (gradient[0][0] + gradient[1][0] + gradient[2][0] + gradient[3][0]);
        B->TRG[1] += (gradient[0][1] + gradient[1][1] + gradient[2][1] + gradient[3][1]);
        B->TRG[2] += (gradient[0][2] + gradient[1][2] + gradient[2][2] + gradient[3][2]);
    }
}

