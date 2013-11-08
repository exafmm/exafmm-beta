#ifndef boundbox_h
#define boundbox_h
#include "logger.h"
#include "thread.h"
#include "types.h"

class BoundBox : public Logger {
 private:
  int nspawn;                                                   //!< Threshold of NBODY for spawning new threads

  //! Recursive functor for calculating bounds
  struct boundsFunctor {
    B_iter BiBegin; B_iter BiEnd; Bounds & bounds; int nspawn;  // State variables
    boundsFunctor(B_iter _BiBegin, B_iter _BiEnd, Bounds & _bounds, int _nspawn) : // Constructor
      BiBegin(_BiBegin), BiEnd(_BiEnd), bounds(_bounds), nspawn(_nspawn) {}// Initialize state variables
    Bounds operator() () {                                      // Overload operator()
      assert(BiEnd - BiBegin > 0);                              // Validate range
      if (BiEnd - BiBegin < nspawn) {                           // If number of elements is small enough
	for (B_iter B=BiBegin; B!=BiEnd; B++) {                 //  Loop over range of bodies
	  bounds.Xmin = min(B->X, bounds.Xmin);                 //   Update Xmin
	  bounds.Xmax = max(B->X, bounds.Xmax);                 //   Update Xmax
	}                                                       //  End loop over range of bodies
	return bounds;                                          //  Return Xmin and Xmax as pair
      } else {                                                  // Else if number of elements are large
	B_iter BiMid = BiBegin + (BiEnd - BiBegin) / 2;         //  Middle iterator
	task_group;                                             //  Initialize tasks
        boundsFunctor leftBounds(BiBegin, BiMid, bounds, nspawn);// Recursion for left branch
	create_taskc(leftBounds);                               //  Create task for left branch
        boundsFunctor rightBounds(BiMid, BiEnd, bounds, nspawn);// Recursion for right branch
	Bounds bounds2 = rightBounds();                         //  Use same task for right branch
	wait_tasks;                                             //  Synchronize tasks
	bounds.Xmin = min(bounds.Xmin, bounds2.Xmin);           //  Minimum of the two Xmins
	bounds.Xmax = max(bounds.Xmax, bounds2.Xmax);           //  Maximum of the two Xmaxs
	return bounds;                                          //  Return Xmin and Xmax
      }                                                         // End if for number fo elements
    }
  };

 public:
  BoundBox(int _nspawn) : nspawn(_nspawn) {}

  //! Get Xmin and Xmax of domain
  Bounds getBounds(Bodies bodies) {
    startTimer("Get bounds");                                   // Start timer
    Bounds bounds;                                              // Bounds : Contains Xmin, Xmax
    if (bodies.empty()) {                                       // If body vector is empty
      bounds.Xmin = bounds.Xmax = 0;                            //  Set bounds to 0
    } else {                                                    // If body vector is not empty
      bounds.Xmin = bounds.Xmax = bodies.front().X;             //  Initialize Xmin, Xmax
      boundsFunctor boundsRecursion(bodies.begin(),bodies.end(),bounds,nspawn);// Initialize recursion
      bounds = boundsRecursion();                               // Recursive call for bounds calculation
    }                                                           // End if for empty body vector
    stopTimer("Get bounds");                                    // Stop timer
    return bounds;                                              // Return Xmin and Xmax
  }

  //! Update Xmin and Xmax of domain
  Bounds getBounds(Bodies bodies, Bounds bounds) {
    startTimer("Get bounds");                                   // Start timer
    boundsFunctor boundsRecursion(bodies.begin(),bodies.end(),bounds,nspawn);// Initialize recursion
    bounds = boundsRecursion();                                 // Recursive call for bounds calculation
    stopTimer("Get bounds");                                    // Stop timer
    return bounds;                                              // Return Xmin and Xmax
  }
};
#endif
