#ifndef boundbox_h
#define boundbox_h
#include "logger.h"
#include "thread.h"
#include "types.h"

class BoundBox : public Logger {
 private:
  int NSPAWN;                                                   //!< Threshold of NDBODY for spawning new threads

//! Recursively get Xmin and Xmax of domain
  Bounds boundsRecursion(B_iter BiBegin, B_iter BiEnd) {
    assert(BiEnd - BiBegin > 0);
    if (BiEnd - BiBegin < NSPAWN) {                             // If number of elements is small enough
      Bounds bounds;
      bounds.Xmin = BiBegin->X;                                 //  Initialize Xmin with first value
      bounds.Xmax = BiBegin->X;                                 //  Initialize Xmax with first value
      for (B_iter B=BiBegin; B!=BiEnd; B++) {                   //  Loop over range of bodies
        bounds.Xmin = min(B->X, bounds.Xmin);                   //   Update Xmin
        bounds.Xmax = max(B->X, bounds.Xmax);                   //   Update Xmax
      }                                                         //  End loop over range of bodies
      return bounds;                                            //  Return Xmin and Xmax as pair
    } else {                                                    // Else if number of elements are large
      B_iter BiMid = BiBegin + (BiEnd - BiBegin) / 2;           //  Middle iterator
      Bounds bounds0, bounds1;                                  //  Pair : first Xmin : second Xmax
      spawn_tasks {                                             //  Initialize tasks
	spawn_task1(bounds0, bounds0 = boundsRecursion(BiBegin, BiMid));//  Recursive call with new task
	bounds1 = boundsRecursion(BiMid, BiEnd);                //  Recursive call with old task
	sync_tasks;                                             //  Synchronize tasks
	bounds0.Xmin = min(bounds0.Xmin, bounds1.Xmin);         //  Minimum of the two Xmins
	bounds0.Xmax = max(bounds0.Xmax, bounds1.Xmax);         //  Maximum of the two Xmaxs
      }
      return bounds0;                                           //  Return Xmin and Xmax
    }                                                           // End if for number fo elements
  }

 public:
  BoundBox(int nspawn) : NSPAWN(nspawn) {}
  ~BoundBox() {}

  // ! Get Xmin and Xmax of domain
  Bounds getBounds(Bodies bodies) {
    startTimer("Get bounds");                                   // Start timer
    Bounds bounds = boundsRecursion(bodies.begin(),bodies.end());// Recursive call for bounds calculation
    stopTimer("Get bounds",printNow);                           // Stop timer
    return bounds;                                              // Return Xmin and Xmax
  }

//! Transform Xmin & Xmax to X (center) & R (radius)
  Box bounds2box(Bounds bounds) {
    vec3 Xmin = bounds.Xmin;                                   // Set local Xmin
    vec3 Xmax = bounds.Xmax;                                  // Set local Xmax
    Box box;                                                    // Bounding box
    for (int d=0; d<3; d++) box.X[d] = (Xmax[d] + Xmin[d]) / 2; // Calculate center of domain
    box.R = 0;                                                  // Initialize localRadius
    for (int d=0; d<3; d++) {                                   // Loop over dimensions
      box.R = std::max(box.X[d] - Xmin[d], box.R);              //  Calculate min distance from center
      box.R = std::max(Xmax[d] - box.X[d], box.R);              //  Calculate max distance from center
    }                                                           // End loop over dimensions
    box.R *= 1.00001;                                           // Add some leeway to radius
    return box;                                                 // Return box.X and box.R
  }
};
#endif
