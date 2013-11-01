#ifndef vanderwaals_h
#define vanderwaals_h
#include "types.h"

class VanDerWaals : public Logger {
 private:
  real_t cuton;                                                 //!< Cuton distance
  real_t cutoff;                                                //!< Cutoff distance
  real_t cycle;                                                 //!< Periodic cycle

 private:
//! Leonard-Jones P2P kernel
  void P2P(C_iter Ci, C_iter Cj, vec3 Xperiodic) const {
    for (B_iter Bi=Ci->BODY; Bi!=Ci->BODY+Ci->NBODY; Bi++) {    // Loop over target bodies
      for (B_iter Bj=Cj->BODY; Bj!=Cj->BODY+Cj->NBODY; Bj++) {  //  Loop over source bodies
	vec3 dist = Bi->X - Bj->X - Xperiodic;                  //   Distance vector from source to target
	real_t R2 = norm(dist);                                 //   R^2
	if (0 < R2 && R2 < cutoff * cutoff) {                   //   Exclude self interaction
	}                                                       //   End if for self interaction
      }                                                         //  End loop over source bodies
    }                                                           // End loop over target bodies
  }

//! Traverse tree to find neighbors
  void traverse(C_iter Ci, C_iter Cj, C_iter C0, vec3 Xperiodic) const {
    vec3 dX = Ci->X - Cj->X - Xperiodic;                        // Distance vector from source to target
    real_t R = std::sqrt(norm(dX));                             // Scalar distance
    if (R < 3 * cutoff) {                                       // If cells are close
      if(Cj->NCHILD == 0) P2P(Ci,Cj,Xperiodic);                 //  L-J kernel
      for (C_iter CC=C0+Cj->ICHILD; CC!=C0+Cj->ICHILD+Cj->NCHILD; CC++) {// Loop over cell's children
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
  VanDerWaals(int _cuton, real_t _cutoff, real_t _cycle) :
    cuton(_cuton), cutoff(_cutoff), cycle(_cycle) {}

//! Evaluate Van Der Waals potential and force
  void evaluate(Cells &cells, Cells &jcells) {
    startTimer("Van Der Waals");                                // Start timer
    C_iter Cj = jcells.begin();                                 // Set begin iterator for source cells
    spawn_tasks {                                               // Intitialize tasks
      for (C_iter Ci=cells.begin(); Ci!=cells.end(); Ci++) {    //  Loop over target cells
        spawn_task0(if (Ci->NCHILD == 0) neighbor(Ci,Cj));      //   Find neighbors of leaf cells
      }                                                         //  End loop over target cells
      sync_tasks;                                               //  Synchronize tasks
    }                                                           // Finalize tasks
    stopTimer("Van Der Waals");                                 // Stop timer
  }

  void print(int stringLength) {
    if (verbose) {
      std::cout << std::setw(stringLength) << std::fixed << std::left// Set format
                << "cuton" << " : " << cuton << std::endl       // Print cuton
                << std::setw(stringLength)                      // Set format
                << "cutoff" << " : " << cutoff << std::endl     // Print cutoff
                << std::setw(stringLength)                      // Set format
                << "cycle" << " : " << cycle << std::endl;      // Print cycle
    }
  }
};
#endif
