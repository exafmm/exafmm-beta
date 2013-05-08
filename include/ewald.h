#ifndef ewald_h
#define ewald_h
#include <stack>
#include "types.h"

class Ewald {
//! Wave structure for Ewald summation
  struct Wave {
    vec3   K;                                                   //!< 3-D wave number vector
    real_t REAL;                                                //!< real part of wave
    real_t IMAG;                                                //!< imaginary part of wave
  };
  typedef std::vector<Wave>           Waves;                    //!< Vector of Wave types
  typedef std::vector<Wave>::iterator W_iter;                   //!< Iterator for Wave types

 private:
  real_t KSIZE;                                                 //!< Number of waves in Ewald summation
  real_t ALPHA;                                                 //!< Scaling parameter for Ewald summation
  real_t SIGMA;                                                 //!< Scaling parameter for Ewald summation
  real_t THETA;                                                 //!< Neighbor acceptance criteria
  real_t CYCLE;                                                 //!< Periodic cycle
  vec3   Xperiodic;                                             //!< Coordinate offset for periodic B.C.

 private:
//! Forward DFT
  void dft(Waves &waves, Bodies &bodies) {
    real_t scale = 2 * M_PI / CYCLE;                            // Scale conversion
    for (W_iter W=waves.begin(); W!=waves.end(); W++) {         // Loop over waves
      W->REAL = W->IMAG = 0;                                    //  Initialize waves
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {     //  Loop over bodies
	real_t th = 0;                                          //   Initialize phase
	for (int d=0; d<3; d++) th += W->K[d] * B->X[d] * scale;//   Determine phase
	W->REAL += B->SRC * cos(th);                            //   Accumulate real component
	W->IMAG += B->SRC * sin(th);                            //   Accumulate imaginary component
      }                                                         //  End loop over bodies
    }                                                           // End loop over waves
  }

//! Inverse DFT
  void idft(Waves &waves, Bodies &bodies) {
    real_t scale = 2 * M_PI / CYCLE;                            // Scale conversion
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      kvec4 TRG = 0;                                            //  Initialize target values
      for (W_iter W=waves.begin(); W!=waves.end(); W++) {       //   Loop over waves
	real_t th = 0;                                          //    Initialzie phase
	for (int d=0; d<3; d++) th += W->K[d] * B->X[d] * scale;//    Determine phase
	real_t dtmp = W->REAL * sin(th) - W->IMAG * cos(th);    //    Temporary value
	TRG[0]     += W->REAL * cos(th) + W->IMAG * sin(th);    //    Accumulate potential
	for (int d=0; d<3; d++) TRG[d+1] -= dtmp * W->K[d];     //    Accumulate force
      }                                                         //   End loop over waves
      B->TRG = TRG;                                             //  Copy results to bodies
    }                                                           // End loop over bodies
  }

//! Initialize wave vector
  Waves initWaves() {
    Waves waves;                                                // Initialzie wave vector
    real_t kmaxsq = KSIZE * KSIZE;                              // kmax squared
    int kmax = int(KSIZE);                                      // kmax as integer
    for (int l=0; l<=kmax; l++) {                               // Loop over x component
      int mmin = -kmax;                                         //  Determine minimum y component
      if (l==0) mmin = 0;                                       //  Exception for minimum y component
      for (int m=mmin; m<=kmax; m++) {                          //  Loop over y component
	int nmin = -kmax;                                       //   Determine minimum z component
	if (l==0 && m==0) nmin=1;                               //   Exception for minimum z component
	for (int n=nmin; n<=kmax; n++) {                        //   Loop over z component
	  real_t ksq = l * l + m * m + n * n;                   //    Wave number squared
	  if (ksq <= kmaxsq) {                                  //    If wave number is below kmax
	    Wave wave;                                          //     Initialzie wave structure
	    wave.K[0] = l;                                      //     x component of k
	    wave.K[1] = m;                                      //     y component of k
	    wave.K[2] = n;                                      //     z component of k
	    wave.REAL = wave.IMAG = 0;                          //     Initialize amplitude
	    waves.push_back(wave);                              //     Push wave to vector
	  }                                                     //    End if for wave number
	}                                                       //   End loop over z component
      }                                                         //  End loop over y component
    }                                                           // End loop over x component
    return waves;                                               // Return wave vector
  }
//! Ewald real part P2P kernel
  void P2P(C_iter Ci, C_iter Cj) {
    for (B_iter Bi=Ci->BODY; Bi!=Ci->BODY+Ci->NDBODY; Bi++) {   // Loop over target bodies
      for (B_iter Bj=Cj->BODY; Bj!=Cj->BODY+Cj->NDBODY; Bj++) { //  Loop over source bodies
	vec3 dist = Bi->X - Bj->X - Xperiodic;                  //   Distance vector from source to target
	for (int d=0; d<3; d++) {                               //   Loop over dimensions
	  if (dist[d] < -CYCLE/2) dist[d] += CYCLE;             //    Wrap domain so that target is always at
	  if (dist[d] >= CYCLE/2) dist[d] -= CYCLE;             //     the center of a [-CYCLE/2,CYCLE/2]^3 source cube
	}                                                       //   End loop over dimensions
	real_t R2 = norm(dist);                                 //   R^2
	if (R2 != 0) {                                          //   Exclude self interaction
	  real_t R2s = R2 * ALPHA * ALPHA;                      //    (R * alpha)^2
	  real_t Rs = std::sqrt(R2s);                           //    R * alpha
	  real_t invRs = 1 / Rs;                                //    1 / (R * alpha)
	  real_t invR2s = invRs * invRs;                        //    1 / (R * alpha)^2
	  real_t invR3s = invR2s * invRs;                       //    1 / (R * alpha)^3
	  real_t dtmp = Bj->SRC * (M_2_SQRTPI * exp(-R2s) * invR2s + erfc(Rs) * invR3s);
	  dtmp *= ALPHA * ALPHA * ALPHA;                        //    Scale temporary value
	  Bi->TRG[0] += Bj->SRC * erfc(Rs) * invRs * ALPHA;     //    Ewald real potential
	  Bi->TRG[1] -= dist[0] * dtmp;                         //    x component of Ewald real force
	  Bi->TRG[2] -= dist[1] * dtmp;                         //    y component of Ewald real force
	  Bi->TRG[3] -= dist[2] * dtmp;                         //    z component of Ewald real force
	}                                                       //   End if for self interaction
      }                                                         //  End loop over source bodies
    }                                                           // End loop over target bodies
  }

//! Traverse tree to find neighbors
  void traverse(C_iter Ci, C_iter C, C_iter C0) {
    std::stack<C_iter> cellStack;                               // Stack of cell iterators
    cellStack.push(C);                                          // Push pair to queue
    while (!cellStack.empty()) {                                // While traversal stack is not empty
      C = cellStack.top();                                      //  Get cell from top of stack
      cellStack.pop();                                          //  Pop traversal stack
      for (C_iter Cj=C0+C->CHILD; Cj!=C0+C->CHILD+C->NCHILD; Cj++) {// Loop over cell's children
        vec3 dX = Ci->X - Cj->X - Xperiodic;                    //   Distance vector from source to target
        real_t Rq = std::sqrt(norm(dX));                        //   Scalar distance
        if (Rq * THETA < Ci->R + Cj->R && Cj->NCHILD == 0) {    //   If twigs are close
          P2P(Ci,Cj);                                           //    Ewald real part
        } else if (Cj->NCHILD != 0) {                           //   If cells are not twigs
          cellStack.push(Cj);                                   //    Push source cell to stack
        }                                                       //   End if for twig cells
      }                                                         //  End loop over cell's children
    }                                                           // End while loop for traversal stack
  }

 public:
//! Constructor
 Ewald(real_t ksize, real_t alpha, real_t sigma, real_t theta, real_t cycle) :
  KSIZE(ksize), ALPHA(alpha), SIGMA(sigma), THETA(theta), CYCLE(cycle), Xperiodic(0) {}

//! Rescale positions and charges
  Bounds rescale(Bodies &bodies) {
    const int numBodies = bodies.size();                        // Number of bodies
    real_t average = 0;                                         // Initialize average charge
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      B->X *= CYCLE / (2 * M_PI);                               //  Rescale positions
      B->X += CYCLE / 2;                                        //  Shift positions
      B->SRC = drand48() / numBodies;                           //  Randomize charges
      B->IBODY = B-bodies.begin();                              //  Initial body numbering for sorting back
      average += B->SRC;                                        //  Accumulate charges
    }                                                           // End loop over bodies
    average /= numBodies;                                       // Divide by total to get average
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      B->SRC -= average;                                        //  Make average charge 0
    }                                                           // End loop over bodies
    Bounds bounds;                                              // Initialize bounds
    bounds.Xmin = 0;                                            // Force Xmin to be -CYCLE/2
    bounds.Xmax = CYCLE;                                        // Force Xmax to be CYCLE/2
    return bounds;                                              // Return bounds
  }

//! Ewald real part
  void realPart(Cells &cells, Cells &jcells) {
    C_iter Cj = jcells.begin();                                 // Set begin iterator for source cells
    for (C_iter Ci=cells.begin(); Ci!=cells.end(); Ci++) {      // Loop over target cells
      if (Ci->NCHILD == 0) {                                    //  If cell is a twig
        for (int ix=-1; ix<=1; ix++) {                          //   Loop over x periodic direction
          for (int iy=-1; iy<=1; iy++) {                        //    Loop over y periodic direction
            for (int iz=-1; iz<=1; iz++) {                      //     Loop over z periodic direction
              Xperiodic[0] = ix * CYCLE;                        //      Coordinate offset for x periodic direction
              Xperiodic[1] = iy * CYCLE;                        //      Coordinate offset for y periodic direction
              Xperiodic[2] = iz * CYCLE;                        //      Coordinate offset for z periodic direction
              traverse(Ci,Cj,Cj);                               //      Traverse the source tree
            }                                                   //     End loop over z periodic direction
          }                                                     //    End loop over y periodic direction
        }                                                       //   End loop over x periodic direction
        for (B_iter B=Ci->BODY; B!=Ci->BODY+Ci->NCBODY; B++) {  //   Loop over all bodies in cell
          B->TRG[0] -= M_2_SQRTPI * B->SRC * ALPHA;             //    Self term of Ewald real part
        }                                                       //   End loop over all bodies in cell
      }                                                         //  End if for twig cells
    }                                                           // End loop over target cells
  }

//! Ewald wave part
  void wavePart(Bodies &bodies) {
    Waves waves = initWaves();                                  // Initialize wave vector
    dft(waves,bodies);                                          // Apply DFT to bodies to get waves
    real_t scale = 2 * M_PI / CYCLE;                            // Scale conversion
    real_t coef = .5 / M_PI / M_PI / SIGMA / CYCLE;             // First constant
    real_t coef2 = scale * scale / (4 * ALPHA * ALPHA);         // Second constant
    for (W_iter W=waves.begin(); W!=waves.end(); W++) {         // Loop over waves
      real_t K2 = norm(W->K);                                   //  Wave number squared
      real_t factor = coef * exp(-K2 * coef2) / K2;             //  Wave factor
      W->REAL *= factor;                                        //  Apply wave factor to real part
      W->IMAG *= factor;                                        //  Apply wave factor to imaginary part
    }                                                           // End loop over waves
    idft(waves,bodies);                                         // Inverse DFT
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      for (int d=0; d<3; d++) B->TRG[d+1] *= scale;             //  Scale forces
    }                                                           // End loop over bodies

    vec3 dipole = 0;                                            // Initialize dipole correction
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      dipole += (B->X - CYCLE/2) * B->SRC;                      //  Calcuate dipole of the whole system
    }                                                           // End loop over bodies
    coef = 4 * M_PI / (3 * CYCLE * CYCLE * CYCLE);              // Precalcualte constant
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      B->TRG[0] += coef * norm(dipole) / bodies.size() / B->SRC;//  Dipole correction for potential
      for (int d=0; d!=3; d++) {                                //  Loop over dimensions
	B->TRG[d+1] += coef * dipole[d];                        //   Dipole correction for forces
      }                                                         //  End loop over dimensions
    }                                                           // End loop over bodies
  }

//! Evaluate relaitve L2 norm error
  void evalError(Bodies &bodies, Bodies &bodies2,
                 double &diff1, double &norm1, double &diff2, double &norm2) {
    double p = 0, p2 = 0;                                       // Initialize total energy
    B_iter B2 = bodies2.begin();                                // Set iterator for bodies2
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++, B2++) { // Loop over bodies & bodies2
      p += B->TRG[0] * B->SRC;                                  //  Accumulate energy from bodies
      p2 += B2->TRG[0] * B2->SRC;                               //  Accumulate energy from bodies2
      double df = 0, f = 0;                                     //  Initialize difference and value
      df += (B->TRG[1] - B2->TRG[1]) * (B->TRG[1] - B2->TRG[1]);//  Difference of x acceleration
      df += (B->TRG[2] - B2->TRG[2]) * (B->TRG[2] - B2->TRG[2]);//  Difference of y acceleration
      df += (B->TRG[3] - B2->TRG[3]) * (B->TRG[3] - B2->TRG[3]);//  Difference of z acceleration
      f += B2->TRG[1] * B2->TRG[1];                             //  Value of x acceleration
      f += B2->TRG[2] * B2->TRG[2];                             //  Value of y acceleration
      f += B2->TRG[3] * B2->TRG[3];                             //  Value of z acceleration
      diff2 += df;                                              //  Accumulate difference of force
      norm2 += f;                                               //  Accumulate value of force
    }                                                           // End loop over bodies & bodies2
    diff1 += (p - p2) * (p - p2);                               // Difference of potential
    norm1 += p2 * p2;                                           // Value of potential
  }
};
#endif
