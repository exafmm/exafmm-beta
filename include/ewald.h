#ifndef ewald_h
#define ewald_h
#include "types.h"

class Ewald : public Logger {
//! Wave structure for Ewald summation
  struct Wave {
    vec3   K;                                                   //!< 3-D wave number vector
    real_t REAL;                                                //!< real part of wave
    real_t IMAG;                                                //!< imaginary part of wave
  };
  typedef std::vector<Wave>           Waves;                    //!< Vector of Wave types
  typedef std::vector<Wave>::iterator W_iter;                   //!< Iterator for Wave types

 private:
  int ksize;                                                    //!< Number of waves in Ewald summation
  real_t alpha;                                                 //!< Scaling parameter for Ewald summation
  real_t sigma;                                                 //!< Scaling parameter for Ewald summation
  real_t cutoff;                                                //!< Neighbor acceptance criteria
  real_t cycle;                                                 //!< Periodic cycle

 private:
//! Forward DFT
  void dft(Waves &waves, Bodies &bodies) const {
    real_t scale = 2 * M_PI / cycle;                            // Scale conversion
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
  void idft(Waves &waves, Bodies &bodies) const {
    real_t scale = 2 * M_PI / cycle;                            // Scale conversion
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      kvec4 TRG = 0;                                            //  Initialize target values
      for (W_iter W=waves.begin(); W!=waves.end(); W++) {       //   Loop over waves
	real_t th = 0;                                          //    Initialzie phase
	for (int d=0; d<3; d++) th += W->K[d] * B->X[d] * scale;//    Determine phase
	real_t dtmp = W->REAL * sin(th) - W->IMAG * cos(th);    //    Temporary value
	TRG[0]     += W->REAL * cos(th) + W->IMAG * sin(th);    //    Accumulate potential
	for (int d=0; d<3; d++) TRG[d+1] -= dtmp * W->K[d];     //    Accumulate force
      }                                                         //   End loop over waves
      for (int d=0; d<3; d++) TRG[d+1] *= scale;                //   Scale forces
      B->TRG += TRG;                                            //  Copy results to bodies
    }                                                           // End loop over bodies
  }

//! Initialize wave vector
  Waves initWaves() const {
    Waves waves;                                                // Initialzie wave vector
    int kmaxsq = ksize * ksize;                                 // kmax squared
    int kmax = ksize;                                           // kmax as integer
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
  void P2P(C_iter Ci, C_iter Cj, vec3 Xperiodic) const {
    for (B_iter Bi=Ci->BODY; Bi!=Ci->BODY+Ci->NDBODY; Bi++) {   // Loop over target bodies
      for (B_iter Bj=Cj->BODY; Bj!=Cj->BODY+Cj->NDBODY; Bj++) { //  Loop over source bodies
	vec3 dist = Bi->X - Bj->X - Xperiodic;                  //   Distance vector from source to target
	for (int d=0; d<3; d++) {                               //   Loop over dimensions
	  if (dist[d] < -cycle/2) dist[d] += cycle;             //    Wrap domain so that target is always at
	  if (dist[d] >= cycle/2) dist[d] -= cycle;             //     the center of a [-cycle/2,cycle/2]^3 source cube
	}                                                       //   End loop over dimensions
	real_t R2 = norm(dist);                                 //   R^2
	if (R2 != 0) {                                          //   Exclude self interaction
	  real_t R2s = R2 * alpha * alpha;                      //    (R * alpha)^2
	  real_t Rs = std::sqrt(R2s);                           //    R * alpha
	  real_t invRs = 1 / Rs;                                //    1 / (R * alpha)
	  real_t invR2s = invRs * invRs;                        //    1 / (R * alpha)^2
	  real_t invR3s = invR2s * invRs;                       //    1 / (R * alpha)^3
	  real_t dtmp = Bj->SRC * (M_2_SQRTPI * exp(-R2s) * invR2s + erfc(Rs) * invR3s);
	  dtmp *= alpha * alpha * alpha;                        //    Scale temporary value
	  Bi->TRG[0] += Bj->SRC * erfc(Rs) * invRs * alpha;     //    Ewald real potential
	  Bi->TRG[1] -= dist[0] * dtmp;                         //    x component of Ewald real force
	  Bi->TRG[2] -= dist[1] * dtmp;                         //    y component of Ewald real force
	  Bi->TRG[3] -= dist[2] * dtmp;                         //    z component of Ewald real force
	}                                                       //   End if for self interaction
      }                                                         //  End loop over source bodies
    }                                                           // End loop over target bodies
  }

//! Traverse tree to find neighbors
  void traverse(C_iter Ci, C_iter Cj, C_iter C0, vec3 Xperiodic) const {
    vec3 dX = Ci->X - Cj->X - Xperiodic;                        // Distance vector from source to target
    real_t R = std::sqrt(norm(dX));                             // Scalar distance
    if (R < cutoff && Cj->NCHILD == 0) {                        // If cells are close
      P2P(Ci,Cj,Xperiodic);                                     //  Ewald real part
    } else {                                                    // If cells are far
      for (C_iter CC=C0+Cj->CHILD; CC!=C0+Cj->CHILD+Cj->NCHILD; CC++) {// Loop over cell's children
        traverse(Ci,CC,C0,Xperiodic);                           //   Recursively call traverse
      }                                                         //  End loop over cell's children
    }                                                           // End if for far cells
  }

//! Find neighbor cells
  void neighbor(C_iter Ci, C_iter Cj) const {
    vec3 Xperiodic;                                             //  Coordinate offset for periodic B.C.
    for (int ix=-1; ix<=1; ix++) {                              //  Loop over x periodic direction
      for (int iy=-1; iy<=1; iy++) {                            //   Loop over y periodic direction
	for (int iz=-1; iz<=1; iz++) {                          //    Loop over z periodic direction
	  Xperiodic[0] = ix * cycle;                            //     Coordinate offset for x periodic direction
          Xperiodic[1] = iy * cycle;                            //     Coordinate offset for y periodic direction
	  Xperiodic[2] = iz * cycle;                            //     Coordinate offset for z periodic direction
	  traverse(Ci,Cj,Cj,Xperiodic);                         //     Traverse the source tree
	}                                                       //    End loop over z periodic direction
      }                                                         //   End loop over y periodic direction
    }                                                           //  End loop over x periodic direction
  }

 public:
//! Constructor
 Ewald(int _ksize, real_t _alpha, real_t _sigma, real_t _cutoff, real_t _cycle) :
  ksize(_ksize), alpha(_alpha), sigma(_sigma), cutoff(_cutoff), cycle(_cycle) {}

//! Ewald real part
  void realPart(Cells &cells, Cells &jcells) {
    startTimer("Ewald real part");                              // Start timer
    C_iter Cj = jcells.begin();                                 // Set begin iterator for source cells
    spawn_tasks {                                               // Intitialize tasks
      for (C_iter Ci=cells.begin(); Ci!=cells.end(); Ci++) {    //  Loop over target cells
        spawn_task0(if (Ci->NCHILD == 0) neighbor(Ci,Cj));      //   Find neighbors of leaf cells
      }                                                         //  End loop over target cells
      sync_tasks;                                               //  Synchronize tasks
    }                                                           // Finalize tasks
    stopTimer("Ewald real part");                               // Stop timer
  }

//! Subtract self term
  void selfTerm(Bodies &bodies) {
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       //  Loop over all bodies
      B->TRG[0] -= M_2_SQRTPI * B->SRC * alpha;                 //   Self term of Ewald real part
    }                                                           //  End loop over all bodies in cell
  }

//! Ewald wave part
  void wavePart(Bodies &bodies, Bodies &jbodies) {
    startTimer("Ewald wave part");                              // Start timer
    Waves waves = initWaves();                                  // Initialize wave vector
    dft(waves,jbodies);                                         // Apply DFT to bodies to get waves
    real_t scale = 2 * M_PI / cycle;                            // Scale conversion
    real_t coef = .5 / M_PI / M_PI / sigma / cycle;             // First constant
    real_t coef2 = scale * scale / (4 * alpha * alpha);         // Second constant
    for (W_iter W=waves.begin(); W!=waves.end(); W++) {         // Loop over waves
      real_t K2 = norm(W->K);                                   //  Wave number squared
      real_t factor = coef * exp(-K2 * coef2) / K2;             //  Wave factor
      W->REAL *= factor;                                        //  Apply wave factor to real part
      W->IMAG *= factor;                                        //  Apply wave factor to imaginary part
    }                                                           // End loop over waves
    idft(waves,bodies);                                         // Inverse DFT
    stopTimer("Ewald wave part");                               // Stop timer
  }

  void print(int stringLength) {
    if (verbose) {
      std::cout << std::setw(stringLength) << std::fixed << std::left// Set format
                << "ksize" << " : " << ksize << std::endl       // Print ksize
                << std::setw(stringLength)                      // Set format
                << "alpha" << " : " << alpha << std::endl       // Print alpha
                << std::setw(stringLength)                      // Set format
                << "sigma" << " : " << sigma << std::endl       // Print sigma
                << std::setw(stringLength)                      // Set format
                << "cutoff" << " : " << cutoff << std::endl     // Print cutoff
                << std::setw(stringLength)                      // Set format
                << "cycle" << " : " << cycle << std::endl;      // Print cycle
    }
  }
};
#endif
