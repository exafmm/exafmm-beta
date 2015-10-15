#ifndef ewald_h
#define ewald_h
#include "logger.h"
#include "types.h"

namespace exafmm {
  class Ewald {
    //! Wave structure for Ewald summation
    struct Wave {
      vec3   K;                                                 //!< 3-D wave number vector
      real_t REAL;                                              //!< real part of wave
      real_t IMAG;                                              //!< imaginary part of wave
    };
    typedef std::vector<Wave> Waves;                            //!< Vector of Wave types
    typedef Waves::iterator   W_iter;                           //!< Iterator of Wave types

  private:
    const int ksize;                                            //!< Number of waves in Ewald summation
    const real_t alpha;                                         //!< Scaling parameter for Ewald summation
    const real_t sigma;                                         //!< Scaling parameter for Ewald summation
    const real_t cutoff;                                        //!< Cutoff distance
    const real_t cycle;                                         //!< Periodic cycle

  private:
    //! Forward DFT
    void dft(Waves & waves, Bodies & bodies) const {
      real_t scale = 2 * M_PI / cycle;                          // Scale conversion
#pragma omp parallel for
      for (int w=0; w<int(waves.size()); w++) {                 // Loop over waves
	W_iter W=waves.begin()+w;                               //  Wave iterator
	W->REAL = W->IMAG = 0;                                  //  Initialize waves
	for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {   //  Loop over bodies
	  real_t th = 0;                                        //   Initialize phase
	  for (int d=0; d<3; d++) th += W->K[d] * B->X[d] * scale;//   Determine phase
	  W->REAL += B->SRC * cos(th);                          //   Accumulate real component
	  W->IMAG += B->SRC * sin(th);                          //   Accumulate imaginary component
	}                                                       //  End loop over bodies
      }                                                         // End loop over waves
    }

    //! Inverse DFT
    void idft(Waves & waves, Bodies & bodies) const {
      real_t scale = 2 * M_PI / cycle;                          // Scale conversion
#pragma omp parallel for
      for (int b=0; b<int(bodies.size()); b++) {                // Loop over bodies
	B_iter B=bodies.begin()+b;                              //  Body iterator
	kvec4 TRG = kreal_t(0);                                 //  Initialize target values
	for (W_iter W=waves.begin(); W!=waves.end(); W++) {     //   Loop over waves
	  real_t th = 0;                                        //    Initialzie phase
	  for (int d=0; d<3; d++) th += W->K[d] * B->X[d] * scale;//    Determine phase
	  real_t dtmp = W->REAL * sin(th) - W->IMAG * cos(th);  //    Temporary value
	  TRG[0]     += W->REAL * cos(th) + W->IMAG * sin(th);  //    Accumulate potential
	  for (int d=0; d<3; d++) TRG[d+1] -= dtmp * W->K[d];   //    Accumulate force
	}                                                       //   End loop over waves
	for (int d=0; d<3; d++) TRG[d+1] *= scale;              //   Scale forces
	B->TRG += TRG;                                          //  Copy results to bodies
      }                                                         // End loop over bodies
    }

    //! Initialize wave vector
    Waves initWaves() const {
      Waves waves;                                              // Initialzie wave vector
      int kmaxsq = ksize * ksize;                               // kmax squared
      int kmax = ksize;                                         // kmax as integer
      for (int l=0; l<=kmax; l++) {                             // Loop over x component
	int mmin = -kmax;                                       //  Determine minimum y component
	if (l==0) mmin = 0;                                     //  Exception for minimum y component
	for (int m=mmin; m<=kmax; m++) {                        //  Loop over y component
	  int nmin = -kmax;                                     //   Determine minimum z component
	  if (l==0 && m==0) nmin=1;                             //   Exception for minimum z component
	  for (int n=nmin; n<=kmax; n++) {                      //   Loop over z component
	    real_t ksq = l * l + m * m + n * n;                 //    Wave number squared
	    if (ksq <= kmaxsq) {                                //    If wave number is below kmax
	      Wave wave;                                        //     Initialzie wave structure
	      wave.K[0] = l;                                    //     x component of k
	      wave.K[1] = m;                                    //     y component of k
	      wave.K[2] = n;                                    //     z component of k
	      wave.REAL = wave.IMAG = 0;                        //     Initialize amplitude
	      waves.push_back(wave);                            //     Push wave to vector
	    }                                                   //    End if for wave number
	  }                                                     //   End loop over z component
	}                                                       //  End loop over y component
      }                                                         // End loop over x component
      return waves;                                             // Return wave vector
    }

    //! Ewald real part P2P kernel
    void P2P(C_iter Ci, C_iter Cj, vec3 Xperiodic) const {
      for (B_iter Bi=Ci->BODY; Bi!=Ci->BODY+Ci->NBODY; Bi++) {  // Loop over target bodies
	for (B_iter Bj=Cj->BODY; Bj!=Cj->BODY+Cj->NBODY; Bj++) {//  Loop over source bodies
	  vec3 dX = Bi->X - Bj->X - Xperiodic;                  //   Distance vector from source to target
	  real_t R2 = norm(dX);                                 //   R^2
	  if (0 < R2 && R2 < cutoff * cutoff) {                 //   Exclude self interaction and cutoff
	    real_t R2s = R2 * alpha * alpha;                    //    (R * alpha)^2
	    real_t Rs = std::sqrt(R2s);                         //    R * alpha
	    real_t invRs = 1 / Rs;                              //    1 / (R * alpha)
	    real_t invR2s = invRs * invRs;                      //    1 / (R * alpha)^2
	    real_t invR3s = invR2s * invRs;                     //    1 / (R * alpha)^3
	    real_t dtmp = Bj->SRC * (M_2_SQRTPI * exp(-R2s) * invR2s + erfc(Rs) * invR3s);
	    dtmp *= alpha * alpha * alpha;                      //    Scale temporary value
	    Bi->TRG[0] += Bj->SRC * erfc(Rs) * invRs * alpha;   //    Ewald real potential
	    Bi->TRG[1] -= dX[0] * dtmp;                         //    x component of Ewald real force
	    Bi->TRG[2] -= dX[1] * dtmp;                         //    y component of Ewald real force
	    Bi->TRG[3] -= dX[2] * dtmp;                         //    z component of Ewald real force
	  }                                                     //   End if for self interaction
	}                                                       //  End loop over source bodies
      }                                                         // End loop over target bodies
    }

    //! Recursive functor for traversing tree to find neighbors
    struct Neighbor {
      Ewald * ewald;                                            //!< Ewald object
      C_iter Ci;                                                //!< Iterator of current target cell
      C_iter Cj;                                                //!< Iterator of current source cell
      C_iter C0;                                                //!< Iterator of first source cell
      Neighbor(Ewald * _ewald, C_iter _Ci, C_iter _Cj, C_iter _C0) :// Constructor
	ewald(_ewald), Ci(_Ci), Cj(_Cj), C0(_C0) {}             // Initialize variables
      void operator() () {                                      // Overload operator()
	vec3 dX = Ci->X - Cj->X;                                //  Distance vector from source to target
	wrap(dX, ewald->cycle);                                 //  Wrap around periodic domain
	vec3 Xperiodic = Ci->X - Cj->X - dX;                    //  Coordinate offset for periodic B.C.
	real_t R = std::sqrt(norm(dX));                         //  Scalar distance
	if (R < 3 * ewald->cutoff) {                            //  If cells are close
	  if(Cj->NCHILD == 0) ewald->P2P(Ci,Cj,Xperiodic);      //   Ewald real part
	  for (C_iter CC=C0+Cj->ICHILD; CC!=C0+Cj->ICHILD+Cj->NCHILD; CC++) {// Loop over cell's children
	    Neighbor neighbor(ewald, Ci, CC, C0);               //    Instantiate recursive functor
	    neighbor();                                         //    Recursive call
	  }                                                     //   End loop over cell's children
	}                                                       //  End if for far cells
      }                                                         // End overload operator()
    };

  public:
    //! Constructor
    Ewald(int _ksize, real_t _alpha, real_t _sigma, real_t _cutoff, real_t _cycle) :
      ksize(_ksize), alpha(_alpha), sigma(_sigma), cutoff(_cutoff), cycle(_cycle) {} // Initialize variables

    //! Ewald real part
    void realPart(Cells & cells, Cells & jcells) {
      logger::startTimer("Ewald real part");                    // Start timer
      C_iter Cj = jcells.begin();                               // Set begin iterator of source cells
      mk_task_group;                                            // Intitialize tasks
      for (C_iter Ci=cells.begin(); Ci!=cells.end(); Ci++) {    // Loop over target cells
	if (Ci->NCHILD == 0) {                                  //  If target cell is leaf
	  Neighbor neighbor(this, Ci, Cj, Cj);                  //   Instantiate recursive functor
	  create_taskc(neighbor);                               //   Create task for recursive call
	}                                                       //  End if for leaf target cell
      }                                                         // End loop over target cells
      wait_tasks;                                               // Synchronize tasks
      logger::stopTimer("Ewald real part");                     // Stop timer
    }

    //! Subtract self term
    void selfTerm(Bodies & bodies) {
      for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {     // Loop over all bodies
	B->TRG[0] -= M_2_SQRTPI * B->SRC * alpha;               //  Self term of Ewald real part
      }                                                         // End loop over all bodies in cell
    }

    //! Ewald wave part
    void wavePart(Bodies & bodies, Bodies & jbodies) {
      logger::startTimer("Ewald wave part");                    // Start timer
      Waves waves = initWaves();                                // Initialize wave vector
      dft(waves,jbodies);                                       // Apply DFT to bodies to get waves
      real_t scale = 2 * M_PI / cycle;                          // Scale conversion
      real_t coef = .5 / M_PI / M_PI / sigma / cycle;           // First constant
      real_t coef2 = scale * scale / (4 * alpha * alpha);       // Second constant
      for (W_iter W=waves.begin(); W!=waves.end(); W++) {       // Loop over waves
	real_t K2 = norm(W->K);                                 //  Wave number squared
	real_t factor = coef * exp(-K2 * coef2) / K2;           //  Wave factor
	W->REAL *= factor;                                      //  Apply wave factor to real part
	W->IMAG *= factor;                                      //  Apply wave factor to imaginary part
      }                                                         // End loop over waves
      idft(waves,bodies);                                       // Inverse DFT
      logger::stopTimer("Ewald wave part");                     // Stop timer
    }

    void print(int stringLength) {
      if (logger::verbose) {                                    // If verbose flag is true
	std::cout << std::setw(stringLength) << std::fixed << std::left// Set format
		  << "ksize" << " : " << ksize << std::endl     //  Print ksize
		  << std::setw(stringLength)                    //  Set format
		  << "alpha" << " : " << alpha << std::endl     //  Print alpha
		  << std::setw(stringLength)                    //  Set format
		  << "sigma" << " : " << sigma << std::endl     //  Print sigma
		  << std::setw(stringLength)                    //  Set format
		  << "cutoff" << " : " << cutoff << std::endl   //  Print cutoff
		  << std::setw(stringLength)                    //  Set format
		  << "cycle" << " : " << cycle << std::endl;    //  Print cycle
      }                                                         // End if for verbose flag
    }
  };
}
#endif
