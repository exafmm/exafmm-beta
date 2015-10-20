#ifndef boundbox_h
#define boundbox_h
#include "logger.h"
#include "thread.h"
#include "types.h"

class BoundBox : public Logger {
 private:
  int nspawn;                                                   //!< Threshold of NDBODY for spawning new threads

//! Recursively get Xmin and Xmax of domain
  Bounds boundsRecursion(B_iter BiBegin, B_iter BiEnd, Bounds bounds) {
    assert(BiEnd - BiBegin > 0);
    if (BiEnd - BiBegin < nspawn) {                             // If number of elements is small enough
      for (B_iter B=BiBegin; B!=BiEnd; B++) {                   //  Loop over range of bodies
        bounds.Xmin = min(B->X, bounds.Xmin);                   //   Update Xmin
        bounds.Xmax = max(B->X, bounds.Xmax);                   //   Update Xmax
      }                                                         //  End loop over range of bodies
      return bounds;                                            //  Return Xmin and Xmax as pair
    } else {                                                    // Else if number of elements are large
      B_iter BiMid = BiBegin + (BiEnd - BiBegin) / 2;           //  Middle iterator
      spawn_tasks {                                             //  Initialize tasks
	spawn_task1(bounds, bounds = boundsRecursion(BiBegin, BiMid, bounds));// Recursive call with new task
	Bounds bounds2 = boundsRecursion(BiMid, BiEnd, bounds); //  Recursive call with old task
	sync_tasks;                                             //  Synchronize tasks
	bounds.Xmin = min(bounds.Xmin, bounds2.Xmin);           //  Minimum of the two Xmins
	bounds.Xmax = max(bounds.Xmax, bounds2.Xmax);           //  Maximum of the two Xmaxs
      }
      return bounds;                                            //  Return Xmin and Xmax
    }                                                           // End if for number fo elements
  }

 public:
  BoundBox(int _nspawn) : nspawn(_nspawn) {}

  // ! Get Xmin and Xmax of domain
  Bounds getBounds(Bodies bodies) {
    startTimer("Get bounds");                                   // Start timer
    Bounds bounds;                                              // Bounds : Contains Xmin, Xmax
    if (bodies.empty()) {                                       // If body vector is empty
      bounds.Xmin = bounds.Xmax = 0;                            //  Set bounds to 0
    } else {                                                    // If body vector is not empty
      bounds.Xmin = bounds.Xmax = bodies.front().X;             //  Initialize Xmin, Xmax
      bounds = boundsRecursion(bodies.begin(),bodies.end(),bounds);// Recursive call for bounds calculation
    }                                                           // End if for empty body vector
    stopTimer("Get bounds");                                    // Stop timer
    return bounds;                                              // Return Xmin and Xmax
  }

  // ! Update Xmin and Xmax of domain
  Bounds getBounds(Bodies bodies, Bounds bounds) {
    startTimer("Get bounds");                                   // Start timer
    bounds = boundsRecursion(bodies.begin(),bodies.end(),bounds);// Recursive call for bounds calculation
    stopTimer("Get bounds");                                    // Stop timer
    return bounds;                                              // Return Xmin and Xmax
  }
};
#endif
