#ifndef bound_box_h
#define bound_box_h
#include "logger.h"
#include "thread.h"
#include "types.h"

namespace exafmm {
  class BoundBox {
  private:
    const int nspawn;                                           //!< Threshold of NBODY for spawning new threads

    //! Recursive functor for calculating bounds of bodies
    struct BodiesRecursion {
      B_iter BiBegin;                                           //!< Begin iterator of bodies
      B_iter BiEnd;                                             //!< End iterator of bodies
      Bounds & bounds;                                          //!< Bounds : Contains Xmin, Xmax
      int nspawn;                                               //!< Threshold of NBODY for spawning new threads
      BodiesRecursion(B_iter _BiBegin, B_iter _BiEnd, Bounds & _bounds, int _nspawn) : // Constructor
	BiBegin(_BiBegin), BiEnd(_BiEnd), bounds(_bounds), nspawn(_nspawn) {}// Initialize variables
      void operator() () {                                      // Overload operator()
	assert(BiEnd - BiBegin > 0);                            //  Validate range
	if (BiEnd - BiBegin < nspawn) {                         //  If number of elements is small enough
	  for (B_iter B=BiBegin; B!=BiEnd; B++) {               //   Loop over range of bodies
	    bounds.Xmin = min(B->X, bounds.Xmin - 1e-5);        //    Update Xmin
	    bounds.Xmax = max(B->X, bounds.Xmax + 1e-5);        //    Update Xmax
	  }                                                     //   End loop over range of bodies
	} else {                                                //  Else if number of elements are large
	  B_iter BiMid = BiBegin + (BiEnd - BiBegin) / 2;       //   Middle iterator
	  Bounds bounds2 = bounds;                              //   Copy bounds
	  mk_task_group;                                        //   Initialize tasks
	  BodiesRecursion leftBranch(BiBegin, BiMid, bounds, nspawn);// Instantiate recursive functor
	  create_taskc(leftBranch);                             //   Create new task for left branch
	  BodiesRecursion rightBranch(BiMid, BiEnd, bounds2, nspawn);// Instantiate recursive functor
	  rightBranch();                                        //   Use old task for right branch
	  wait_tasks;                                           //   Synchronize tasks
	  bounds.Xmin = min(bounds.Xmin, bounds2.Xmin);         //   Minimum of the two Xmins
	  bounds.Xmax = max(bounds.Xmax, bounds2.Xmax);         //   Maximum of the two Xmaxs
	}                                                       //  End if for number fo elements
      }                                                         // End overload operator()
    };

    //! Recursive functor for calculating bounds of cells
    struct CellsRecursion {
      C_iter CiBegin;                                           //!< Begin iterator of cells
      C_iter CiEnd;                                             //!< End iterator of cells
      Bounds & bounds;                                          //!< Bounds : Contains Xmin, Xmax
      int nspawn;                                               //!< Threshold of NBODY for spawning new threads
      CellsRecursion(C_iter _CiBegin, C_iter _CiEnd, Bounds & _bounds, int _nspawn) : // Constructor
	CiBegin(_CiBegin), CiEnd(_CiEnd), bounds(_bounds), nspawn(_nspawn) {}// Initialize variables
      void operator() () {                                      // Overload operator()
	assert(CiEnd - CiBegin > 0);                            //  Validate range
	if (CiEnd - CiBegin < nspawn) {                         //  If number of elements is small enough
	  for (C_iter C=CiBegin; C!=CiEnd; C++) {               //   Loop over range of cells
	    bounds.Xmin = min(C->X, bounds.Xmin - 1e-5);        //    Update Xmin
	    bounds.Xmax = max(C->X, bounds.Xmax + 1e-5);        //    Update Xmax
	  }                                                     //   End loop over range of bodies
	} else {                                                //  Else if number of elements are large
	  C_iter CiMid = CiBegin + (CiEnd - CiBegin) / 2;       //   Middle iterator
	  Bounds bounds2 = bounds;                              //   Copy bounds
	  mk_task_group;                                        //   Initialize tasks
	  CellsRecursion leftBranch(CiBegin, CiMid, bounds, nspawn);// Instantiate recursive functor
	  create_taskc(leftBranch);                             //   Create new task for left branch
	  CellsRecursion rightBranch(CiMid, CiEnd, bounds2, nspawn);// Instantiate recursive functor
	  rightBranch();                                        //   Use old task for right branch
	  wait_tasks;                                           //   Synchronize tasks
	  bounds.Xmin = min(bounds.Xmin, bounds2.Xmin);         //   Minimum of the two Xmins
	  bounds.Xmax = max(bounds.Xmax, bounds2.Xmax);         //   Maximum of the two Xmaxs
	}                                                       //  End if for number fo elements
      }                                                         // End overload operator()
    };

  public:
    //! Constructor
    BoundBox(int _nspawn) : nspawn(_nspawn) {}                  // Initialize variables

    //! Get Xmin and Xmax of bodies
    Bounds getBounds(Bodies & bodies) {
      logger::startTimer("Get bounds");                         // Start timer
      Bounds bounds;                                            // Bounds : Contains Xmin, Xmax
      if (bodies.empty()) {                                     // If body vector is empty
	bounds.Xmin = bounds.Xmax = 0;                          //  Set bounds to 0
      } else {                                                  // If body vector is not empty
	bounds.Xmin = bounds.Xmax = bodies.front().X;           //  Initialize Xmin, Xmax
	BodiesRecursion bodiesRecursion(bodies.begin(),bodies.end(),bounds,nspawn);// Instantiate recursive functor
	bodiesRecursion();                                      //  Recursive call for bounds calculation
      }                                                         // End if for empty body vector
      logger::stopTimer("Get bounds");                          // Stop timer
      return bounds;                                            // Return Xmin and Xmax
    }

    //! Update Xmin and Xmax of bodies
    Bounds getBounds(Bodies bodies, Bounds bounds) {
      logger::startTimer("Get bounds");                         // Start timer
      BodiesRecursion bodiesRecursion(bodies.begin(),bodies.end(),bounds,nspawn);// Instantiate recursive functor
      bodiesRecursion();                                        // Recursive call for bounds calculation
      logger::stopTimer("Get bounds");                          // Stop timer
      return bounds;                                            // Return Xmin and Xmax
    }

    //! Get Xmin and Xmax of cells
    Bounds getBounds(Cells cells) {
      logger::startTimer("Get bounds");                         // Start timer
      Bounds bounds;                                            // Bounds : Contains Xmin, Xmax
      if (cells.empty()) {                                      // If cell vector is empty
	bounds.Xmin = bounds.Xmax = 0;                          //  Set bounds to 0
      } else {                                                  // If cell vector is not empty
	bounds.Xmin = bounds.Xmax = cells.front().X;            //  Initialize Xmin, Xmax
	CellsRecursion cellsRecursion(cells.begin(),cells.end(),bounds,nspawn);// Instantiate recursive functor
	cellsRecursion();                                       //  Recursive call for bounds calculation
      }                                                         // End if for empty body vector
      logger::stopTimer("Get bounds");                          // Stop timer
      return bounds;                                            // Return Xmin and Xmax
    }

    //! Update Xmin and Xmax of cells
    Bounds getBounds(Cells cells, Bounds bounds) {
      logger::startTimer("Get bounds");                         // Start timer
      CellsRecursion cellsRecursion(cells.begin(),cells.end(),bounds,nspawn);// Instantiate recursive functor
      cellsRecursion();                                         // Recursive call for bounds calculation
      logger::stopTimer("Get bounds");                          // Stop timer
      return bounds;                                            // Return Xmin and Xmax
    }
  };
}
#endif
