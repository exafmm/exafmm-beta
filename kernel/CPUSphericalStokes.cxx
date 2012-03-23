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
#ifdef STOKES
void Kernel<Stokes>::P2M(C_iter Ci) const
{
    complex Ynm[4*P*P];
    int offsets[3] = {NTERM,2*NTERM,3*NTERM};
    for (B_iter B = Ci->LEAF; B != Ci->LEAF + Ci->NDLEAF; ++B)
    {
        vect dist = B->X - Ci->X;
        real rho, alpha, beta;
        cart2sph(rho, alpha, beta, dist);
        evalMultipole(rho, alpha, -beta, Ynm);
        real f0 = B->FORCE[0];
        real f1 = B->FORCE[1];
        real f2 = B->FORCE[2];
        real fdotx = (B->X[0] * f0 + B->X[1] * f1 + B->X[2] * f2);
        for (int n = 0; n != P; ++n)
        {
            for (int m = 0; m <= n; ++m)
            {
                const int nm  = n * (n + 1) + m;
                const int nms = n * (n + 1) / 2 + m;
                Ci->M[nms] += f0 * Ynm[nm];
                Ci->M[nms+offsets[0]] += f1 * Ynm[nm];
                Ci->M[nms+offsets[1]] += f2 * Ynm[nm];
                Ci->M[nms+offsets[2]] += fdotx * Ynm[nm];
            }
        }
    }
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
        d[0] = B->FORCE[0]*B->X[2];
        d[1] = B->FORCE[1]*B->X[2];
        d[2] = -B->FORCE[2]*B->X[2];
        real u = -B->FORCE[2];
        for (int n = 0; n != P; ++n)
        {
            int nm  = n * n + n;
            int nms = n * (n + 1) / 2;
            
            factor = 1.0 / rho * n;
            gradient[0] += std::real( Ynm[nm]) * factor;
            gradient[1] += std::real( YnmTheta[nm]);
            Ci->M[nms] += u * (gradient[0]*d[0]+gradient[1]*d[1]);
            
            for (int m = 1; m <= n; ++m)
            {
                const int nm  = n * (n + 1) + m;
                const int nms = n * (n + 1) / 2 + m;
                gradient[0] += 2 * std::real( Ynm[nm]) * factor;
                gradient[1] += 2 * std::real( YnmTheta[nm]);
                gradient[2] += 2 * std::real( Ynm[nm] * I) * m;
                Ci->M[nms] += u * (gradient[0]*d[0]+gradient[1]*d[1]+gradient[2]*d[2]);
            }
        }
    }
}
#endif
void Kernel<Stokes>::M2M(C_iter Ci, C_iter Cj) const
{
    const complex I(0., 1.);                                      // Imaginary unit
    complex Ynm[4*P*P];
    vect dist = Ci->X - Cj->X;
    real rho, alpha, beta;
    cart2sph(rho, alpha, beta, dist);
    evalMultipole(rho, alpha, -beta, Ynm);
    int offsets[3] = {NTERM,2*NTERM,3*NTERM};
    for (int j = 0; j != P; ++j)
    {
        for (int k = 0; k <= j; ++k)
        {
            const int jk = j * j + j + k;
            const int jks = j * (j + 1) / 2 + k;
            complex M[4] = {0};
            for (int n = 0; n <= j; ++n)
            {
                for (int m = -n; m <= std::min(k - 1, n); ++m)
                {
                    if (j - n >= k - m)
                    {
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

void Kernel<Stokes>::M2L(C_iter Ci, C_iter Cj) const
{
    complex Ynm[4*P*P];
    vect dist = Ci->X - Cj->X - Xperiodic;
    real rho, alpha, beta;
    cart2sph(rho, alpha, beta, dist);
    evalLocal(rho, alpha, beta, Ynm);
    int offsets[3] = {NTERM,2*NTERM,3*NTERM};
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
    int offsets[3] = {NTERM,2*NTERM,3*NTERM};
    for (B_iter B = Ci->LEAF; B != Ci->LEAF + Ci->NDLEAF; ++B)
    {
        vect dist = B->X - Cj->X - Xperiodic;
        real r, theta, phi, factor;
        vect gradient[4] = {0,0,0,0};
        vect cartesian = 0;
        cart2sph(r, theta, phi, dist);
        evalLocalTheta(r, theta, phi, Ynm, YnmTheta);
        int offsets[3] = {NTERM,2*NTERM,3*NTERM};
        for (int n = 0; n != P; ++n)
        {
            int nm  = n * n + n;
            int nms = n * (n + 1) / 2;
            B->TRG[0] += std::real(Cj->M[nms] * Ynm[nm]);
            B->TRG[1] += std::real(Cj->M[nms+offsets[0]] * Ynm[nm]);
            B->TRG[2] += std::real(Cj->M[nms+offsets[1]] * Ynm[nm]);
            
            factor = 1.0/r*(n+1);
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
        
        B->TRG[0] += (gradient[0][0]+gradient[1][0]+gradient[2][0]+gradient[3][0]);
        B->TRG[1] += (gradient[0][1]+gradient[1][1]+gradient[2][1]+gradient[3][1]);
        B->TRG[2] += (gradient[0][2]+gradient[1][2]+gradient[2][2]+gradient[3][2]);
    }
}

void Kernel<Stokes>::L2L(C_iter Ci, C_iter Cj) const
{
    const complex I(0., 1.);                                      // Imaginary unit
    complex Ynm[4*P*P], factor;
    vect dist = Ci->X - Cj->X;
    real rho, alpha, beta;
    cart2sph(rho, alpha, beta, dist);
    evalMultipole(rho, alpha, beta, Ynm);
    int offsets[3] = {NTERM,2*NTERM,3*NTERM};
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
    int offsets[3] = {NTERM,2*NTERM,3*NTERM};
    for (B_iter B = Ci->LEAF; B != Ci->LEAF + Ci->NDLEAF; ++B)
    {
        vect dist = B->X - Ci->X;
        vect gradient[4] = {0,0,0,0};
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
        
        B->TRG[0] += (gradient[0][0]+gradient[1][0]+gradient[2][0]+gradient[3][0]);
        B->TRG[1] += (gradient[0][1]+gradient[1][1]+gradient[2][1]+gradient[3][1]);
        B->TRG[2] += (gradient[0][2]+gradient[1][2]+gradient[2][2]+gradient[3][2]);
    }
}

