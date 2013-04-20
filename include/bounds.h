#ifndef bounds_h
#define bounds_h
#include "logger.h"
#include "thread.h"
#include "types.h"

class Bounds : public Logger {
 private:
  typedef std::pair<vec3,vec3> vec3Pair;                        //!< Pair of vec3
  int NSPAWN;                                                   //!< Threshold of NDBODY for spawning new threads

 protected:
  int    IMAGES;                                                //!< Number of periodic image sublevels
  fvec3  localXmin;                                             //!< Local Xmin for a given rank
  fvec3  localXmax;                                             //!< Local Xmax for a given rank

 public:
  Box    localBox;                                              //!< Center and radius of local box
  real_t CYCLE;                                                 //!< Periodic cycle

 private:
//! Get Xmin and Xmax of domain
  vec3Pair getLocal(B_iter BiBegin, B_iter BiEnd) {
    assert(BiEnd - BiBegin > 0);
    if (BiEnd - BiBegin < NSPAWN) {                             // If number of elements is small enough
      vec3 Xmin = BiBegin->X, Xmax = BiBegin->X;                //  Initialize Xmin and Xmax with first value
      for (B_iter B=BiBegin; B!=BiEnd; B++) {                   //  Loop over range of bodies
        Xmin = min(B->X, Xmin);                                 //   Update Xmin
        Xmax = max(B->X, Xmax);                                 //   Update Xmax
      }                                                         //  End loop over range of bodies
      return vec3Pair(Xmin, Xmax);                              //  Return Xmin and Xmax as pair
    } else {                                                    // Else if number of elements are large
      B_iter BiMid = BiBegin + (BiEnd - BiBegin) / 2;           //  Middle iterator
      vec3Pair bounds0, bounds1;                                //  Pair : first Xmin : second Xmax
      spawn_tasks {                                             //  Initialize tasks
	spawn_task1(bounds0, bounds0 = getLocal(BiBegin, BiMid));//  Recursive call with new task
	bounds1 = getLocal(BiMid, BiEnd);                       //  Recursive call with old task
	sync_tasks;                                             //  Synchronize tasks
	bounds0.first = min(bounds0.first, bounds1.first);      //  Minimum of the two Xmins
	bounds0.second = max(bounds0.second, bounds1.second);   //  Maximum of the two Xmaxs
      }
      return bounds0;                                           //  Return Xmin and Xmax
    }                                                           // End if for number fo elements
  }

 public:
  Bounds(int nspawn, int images) : NSPAWN(nspawn), IMAGES(images) {}
  ~Bounds() {}
//! Set center and size of root cell
  void setLocal(Bodies &bodies) {
    startTimer("Set bounds");                                   // Start timer
    vec3Pair bounds = getLocal(bodies.begin(), bodies.end());   // Get Xmin (first) and Xmax (second) of domain
    for (int d=0; d<3; d++) localXmin[d] = bounds.first[d];     // Set local Xmin
    for (int d=0; d<3; d++) localXmax[d] = bounds.second[d];    // Set local Xmax
    for (int d=0; d<3; d++) localBox.X[d] = (localXmax[d] + localXmin[d]) / 2;// Calculate center of domain
    localBox.R = 0;                                             // Initialize localRadius
    for (int d=0; d<3; d++) {                                   // Loop over dimensions
      localBox.R = std::max(localBox.X[d] - localXmin[d], localBox.R);// Calculate min distance from center
      localBox.R = std::max(localXmax[d] - localBox.X[d], localBox.R);// Calculate max distance from center
    }                                                           // End loop over dimensions
    localBox.R *= 1.00001;                                      // Add some leeway to radius
    if (IMAGES == 0) {                                          // If non-periodic boundary condition
      CYCLE = 2 * localBox.R;                                   //  Set global radius for serial run
    } else {                                                    // If periodic boundary condition
      CYCLE = 2 * M_PI;                                         //  Set global radius to 2 * pi
    }                                                           // End if for periodic boundary condition
    stopTimer("Set bounds",printNow);
  }
};
#endif
