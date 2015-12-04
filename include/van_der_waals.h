#ifndef van_der_waals_h
#define van_der_waals_h
#include "logger.h"
#include "types.h"

namespace exafmm {
  class VanDerWaals {
  private:
    const real_t cuton;                                         //!< Cuton distance
    const real_t cutoff;                                        //!< Cutoff distance
    const real_t cycle;                                         //!< Periodic cycle
    const int numTypes;                                         //!< Number of atom types
    std::vector<real_t> rscale;                                 //!< Distance scaling parameter for VdW potential
    std::vector<real_t> gscale;                                 //!< Value scaling parameter for VdW potential
    std::vector<real_t> fgscale;                                //!< Value scaling parameter for VdW force

  private:
    //! Van der Waals P2P kernel
    void P2P(C_iter Ci, C_iter Cj, vec3 Xperiodic) const {
      for (B_iter Bi=Ci->BODY; Bi!=Ci->BODY+Ci->NBODY; Bi++) {
	int atypei = int(Bi->SRC);
	for (B_iter Bj=Cj->BODY; Bj!=Cj->BODY+Cj->NBODY; Bj++) {
	  vec3 dX = Bi->X - Bj->X - Xperiodic;
	  real_t R2 = norm(dX);
	  if (R2 != 0) {
	    int atypej = int(Bj->SRC);
	    real_t rs = rscale[atypei*numTypes+atypej];
	    real_t gs = gscale[atypei*numTypes+atypej];
	    real_t fgs = fgscale[atypei*numTypes+atypej];
	    real_t R2s = R2 * rs;
	    real_t invR2 = 1.0 / R2s;
	    real_t invR6 = invR2 * invR2 * invR2;
	    real_t cuton2 = cuton * cuton;
	    real_t cutoff2 = cutoff * cutoff;
	    if (R2 < cutoff2) {
	      real_t tmp = 0, dtmp = 0;
	      if (cuton2 < R2) {
		real_t tmp1 = (cutoff2 - R2) / ((cutoff2-cuton2)*(cutoff2-cuton2)*(cutoff2-cuton2));
		real_t tmp2 = tmp1 * (cutoff2 - R2) * (cutoff2 - 3 * cuton2 + 2 * R2);
		tmp = invR6 * (invR6 - 1) * tmp2;
		dtmp = invR6 * (invR6 - 1) * 12 * (cuton2 - R2) * tmp1
		  - 6 * invR6 * (invR6 + (invR6 - 1) * tmp2) * tmp2 / R2;
	      } else {
		tmp = invR6 * (invR6 - 1);
		dtmp = invR2 * invR6 * (2 * invR6 - 1);
	      }
	      dtmp *= fgs;
	      Bi->TRG[0] += gs * tmp;
	      Bi->TRG[1] -= dX[0] * dtmp;
	      Bi->TRG[2] -= dX[1] * dtmp;
	      Bi->TRG[3] -= dX[2] * dtmp;
	    }
	  }
	}
      }
    }

    //! Recursive functor for traversing tree to find neighbors
    struct Neighbor {
      VanDerWaals * VdW;                                        //!< VanDerWaals object
      C_iter Ci;                                                //!< Iterator of current target cell
      C_iter Cj;                                                //!< Iterator of current source cell
      C_iter C0;                                                //!< Iterator of first source cell
      vec3 Xperiodic;                                           //!< Coordinate offset for periodic B.C.
      Neighbor(VanDerWaals * _VdW, C_iter _Ci, C_iter _Cj, C_iter _C0, vec3 _Xperiodic) :// Constructor
	VdW(_VdW), Ci(_Ci), Cj(_Cj), C0(_C0), Xperiodic(_Xperiodic) {} // Initialize variables
      void operator() () {                                      // Overload operator()
	vec3 dX = Ci->X - Cj->X - Xperiodic;                    //  Coordinate offset for periodic B.C.
	real_t R = std::sqrt(norm(dX));                         //  Scalar distance
	if (R < 3 * VdW->cutoff) {                              //  If cells are close
	  if(Cj->NCHILD == 0) VdW->P2P(Ci, Cj, Xperiodic);      //   Van der Waals kernel
	  for (C_iter CC=C0+Cj->ICHILD; CC!=C0+Cj->ICHILD+Cj->NCHILD; CC++) {// Loop over cell's children
	    Neighbor neighbor(VdW, Ci, CC, C0, Xperiodic);      //    Instantiate recursive functor
	    neighbor();                                         //    Find neighbors recursively
	  }                                                     //   End loop over cell's children
	}                                                       //  End if for far cells
      }                                                         // End overload operator()
    };

  public:
    //! Constructor
    VanDerWaals(double _cuton, double _cutoff, double _cycle, int _numTypes,
		double * _rscale, double * _gscale, double * _fgscale) :
      cuton(_cuton), cutoff(_cutoff), cycle(_cycle), numTypes(_numTypes) {// Initialize variables
      rscale.resize(numTypes*numTypes);
      gscale.resize(numTypes*numTypes);
      fgscale.resize(numTypes*numTypes);
      for (int i=0; i<numTypes*numTypes; i++) {
	rscale[i] = _rscale[i];
	gscale[i] = _gscale[i];
	fgscale[i] = _fgscale[i];
      }
    }

    //! Evaluate Van Der Waals potential and force
    void evaluate(Cells & cells, Cells & jcells) {
      logger::startTimer("Van der Waals");                      // Start timer
      C_iter Cj = jcells.begin();                               // Set begin iterator of source cells
      vec3 Xperiodic;                                           // Coordinate offset for periodic B.C.
      for (int ix=-1; ix<=1; ix++) {                            // Loop over x periodic direction
	for (int iy=-1; iy<=1; iy++) {                          //  Loop over y periodic direction
	  for (int iz=-1; iz<=1; iz++) {                        //   Loop over z periodic direction
	    Xperiodic[0] = ix * cycle;                          //    Coordinate offset for x periodic direction
	    Xperiodic[1] = iy * cycle;                          //    Coordinate offset for y periodic direction
	    Xperiodic[2] = iz * cycle;                          //    Coordinate offset for z periodic direction
	    mk_task_group;                                      //    Intitialize tasks
	    for (C_iter Ci=cells.begin(); Ci!=cells.end(); Ci++) {//  Loop over target cells
	      if (Ci->NCHILD == 0) {                            //     If target cell is leaf
		Neighbor neighbor(this, Ci, Cj, Cj, Xperiodic); //      Instantiate recursive functor
		create_taskc(neighbor);                         //      Create task for recursive call
	      }                                                 //     End if for leaf target cell
	    }                                                   //    End loop over target cells
	    wait_tasks;                                         //    Synchronize tasks
	  }                                                     //   End loop over z periodic direction
	}                                                       //  End loop over y periodic direction
      }                                                         // End loop over x periodic direction
      logger::stopTimer("Van der Waals");                       // Stop timer
    }

    void print(int stringLength) {
      if (logger::verbose) {                                    // If verbose flag is true
	std::cout << std::setw(stringLength) << std::fixed << std::left// Set format
		  << "cuton" << " : " << cuton << std::endl     //  Print cuton
		  << std::setw(stringLength)                    //  Set format
		  << "cutoff" << " : " << cutoff << std::endl   //  Print cutoff
		  << std::setw(stringLength)                    //  Set format
		  << "cycle" << " : " << cycle << std::endl;    //  Print cycle
      }                                                         // End if for verbose flag
    }
  };
}
#endif
