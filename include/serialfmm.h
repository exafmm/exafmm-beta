#ifndef serialfmm_h
#define serialfmm_h
#include "treebuilder.h"

class SerialFMM : public TreeBuilder {
protected:
  real_t globalRadius;                                          //!< Radius of global root cell
  vec3   globalCenter;                                          //!< Center of global root cell
  fvec3  globalXmin;                                            //!< Global Xmin for a given rank
  fvec3  globalXmax;                                            //!< Global Xmax for a given rank

private:
//! Error optimization of Rcrit
  void setRcrit(C_iter C, C_iter C0, real_t c) {
    __init_tasks__;                                             // Initialize tasks
    for (C_iter CC=C0+C->CHILD; CC!=C0+C->CHILD+C->NCHILD; CC++) {// Loop over child cells
      spawn_task0(setRcrit(CC, C0, c));                         //  Recursive call with new task
    }                                                           // End loop over child cells
    __sync_tasks__;                                             // Synchronize tasks
#if Cartesian
    for( int i=1; i<MTERM; ++i ) C->M[i] /= C->M[0];            // Normalize multipole expansion coefficients
#endif
    real_t x = 1.0 / THETA;                                     // Inverse of theta
#if ERROR_OPT
    real_t a = c * pow(std::abs(C->M[0]),1.0/3);                // Cell coefficient
    for (int i=0; i<5; i++) {                                   // Newton-Rhapson iteration
      real_t f = x * x - 2 * x + 1 - a * pow(x,-P);             //  Function value
      real_t df = (P + 2) * x - 2 * (P + 1) + P / x;            //  Function derivative value
      x -= f / df;                                              //  Increment x
    }                                                           // End Newton-Rhapson iteration
#endif
    C->RCRIT *= x;                                              // Multiply Rcrit by error optimized parameter x
  }

  //! Get Xmin and Xmax of domain
  vec3Pair getBounds(B_iter BiBegin, B_iter BiEnd) {
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
      __init_tasks__;                                           //  Initialize tasks
      vec3Pair bounds0, bounds1;                                //  Pair : first Xmin : second Xmax
      spawn_task1(bounds0, bounds0 = getBounds(BiBegin, BiMid));//  Recursive call with new task
      bounds1 = getBounds(BiMid, BiEnd);                        //  Recursive call with old task
      __sync_tasks__;                                           //  Synchronize tasks
      bounds0.first = min(bounds0.first, bounds1.first);        //  Minimum of the two Xmins
      bounds0.second = max(bounds0.second, bounds1.second);     //  Maximum of the two Xmaxs
      return bounds0;                                           //  Return Xmin and Xmax
    }                                                           // End if for number fo elements
  }

//! Recursive call for upward pass
  void upwardRecursion(C_iter C, C_iter C0) {
    __init_tasks__;                                             // Initialize tasks
    for (C_iter CC=C0+C->CHILD; CC!=C0+C->CHILD+C->NCHILD; CC++) {// Loop over child cells
      spawn_task0(upwardRecursion(CC, C0));                      //  Recursive call with new task
    }                                                           // End loop over child cells
    __sync_tasks__;                                             // Synchronize tasks
    real_t Rmax = 0;                                            //  Initialize Rmax
    setCenter(C);                                               //  Set center of cell to center of mass
    C->M = 0;                                                   //  Initialize multipole expansion coefficients
    C->L = 0;                                                   //  Initialize local expansion coefficients
    P2M(C,Rmax);                                                //  P2M kernel
    M2M(C,Rmax);                                                //  M2M kernel
  }

//! Recursive call for downward pass
  void downwardRecursion(C_iter C, C_iter C0) const {
    L2L(C);                                                     // L2L kernel
    L2P(C);                                                     // L2P kernel
    __init_tasks__;                                             // Initialize tasks
    for (C_iter CC=C0+C->CHILD; CC!=C0+C->CHILD+C->NCHILD; CC++) {// Loop over child cells
      spawn_task0(downwardRecursion(CC, C0));                    //  Recursive call with new task
    }                                                           // End loop over chlid cells
    __sync_tasks__;                                             // Synchronize tasks
  }

public:
//! Set center and size of root cell
  void setBounds(Bodies &bodies) {
    startTimer("Set bounds");                                   // Start timer
    vec3Pair bounds = getBounds(bodies.begin(), bodies.end());  // Get Xmin (first) and Xmax (second) of domain
    localXmin = bounds.first;                                   // Set local Xmin
    localXmax = bounds.second;                                  // Set local Xmax
    localCenter = (localXmax + localXmin) / 2;                  // Calculate center of domain
    localRadius = 0;                                            // Initialize localRadius
    for (int d=0; d<3; d++) {                                   // Loop over dimensions
      localRadius = std::max(localCenter[d] - localXmin[d], localRadius);// Calculate min distance from center
      localRadius = std::max(localXmax[d] - localCenter[d], localRadius);// Calculate max distance from center 
    }                                                           // End loop over dimensions
    localRadius *= 1.00001;                                     // Add some leeway to radius
    if (IMAGES == 0) {                                          // If non-periodic boundary condition
      globalRadius = localRadius;                               //  Set global radius for serial run
      globalCenter = localCenter;                               //  Set global center for serial run
      globalXmin = localXmin;                                   //  Set global Xmin for serial run
      globalXmax = localXmax;                                   //  Set global Xmax for serial run
    } else {                                                    // If periodic boundary condition
      globalRadius = M_PI;                                      //  Set global radius
      globalCenter = 0;                                         //  Set global radius
      globalXmin = -M_PI;                                       //  Set global Xmin
      globalXmax = M_PI;                                        //  Set global Xmax
    }                                                           // End if for periodic boundary condition
    stopTimer("Set bounds",printNow);
  }

//! Build tree structure top down
  void buildTree(Bodies &bodies, Cells &cells) {
    growTree(bodies);                                           // Grow tree from root
    linkTree(cells);                                            // Form parent-child links in tree
  }

//! Upward pass (P2M, M2M)
  void upwardPass(Cells &cells) {
    startTimer("Upward pass");                                  // Start timer
    Ci0 = cells.begin();                                        // Set iterator of target root cell
    Cj0 = cells.begin();                                        // Set iterator of source root cell
    upwardRecursion(Ci0, Ci0);                                  // Recursive call for upward pass
    real_t c = (1 - THETA) * (1 - THETA) / pow(THETA,P+2) / pow(std::abs(Ci0->M[0]),1.0/3); // Root coefficient
    setRcrit(Ci0, Ci0, c);                                      // Error optimization of Rcrit
    for (C_iter C=cells.begin(); C!=cells.begin()+9; C++) {     // Loop over top 2 levels of cells
      C->RCRIT *= 10;                                           //  Prevent approximation
    }                                                           // End loop over top 2 levels of cells
    stopTimer("Upward pass",printNow);                          // Stop timer
  }

//! Evaluate P2P and M2L using tree traversal
  void evaluate(Cells &icells, Cells &jcells, bool mutual=false) {
    Ci0 = icells.begin();                                       // Set iterator of target root cell
    Cj0 = jcells.begin();                                       // Set iterator of source root cell
    startTimer("Traverse");                                     // Start timer
    if (IMAGES == 0) {                                          // If non-periodic boundary condition
      Xperiodic = 0;                                            //  No periodic shift
      traverse(Ci0,Cj0,mutual);                                 //  Traverse the tree
    } else {                                                    // If periodic boundary condition
      for (int ix=-1; ix<=1; ix++) {                            //  Loop over x periodic direction
        for (int iy=-1; iy<=1; iy++) {                          //   Loop over y periodic direction
          for (int iz=-1; iz<=1; iz++) {                        //    Loop over z periodic direction
            Xperiodic[0] = ix * 2 * globalRadius;               //     Coordinate shift for x periodic direction
            Xperiodic[1] = iy * 2 * globalRadius;               //     Coordinate shift for y periodic direction
            Xperiodic[2] = iz * 2 * globalRadius;               //     Coordinate shift for z periodic direction
            traverse(Ci0,Cj0,false);                            //     Traverse the tree for this periodic image
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      traversePeriodic(globalRadius);                           //  Traverse tree for periodic images
    }                                                           // End if for periodic boundary condition
    stopTimer("Traverse",printNow);                             // Stop timer
  }

//! Downward pass (L2L, L2P)
  void downwardPass(Cells &cells) { 
    startTimer("Downward pass");                                // Start timer
    C_iter C0 = cells.begin();                                  // Root cell
    L2P(C0);                                                    // If root is the only cell do L2P
    __init_tasks__;                                             // Initialize tasks
    for (C_iter CC=C0+C0->CHILD; CC!=C0+C0->CHILD+C0->NCHILD; CC++) {// Loop over child cells
      spawn_task0(downwardRecursion(CC, C0));                   //  Recursive call for downward pass
    }                                                           // End loop over child cells
    __sync_tasks__;                                             // Synchronize tasks
    stopTimer("Downward pass",printNow);                        // Stop timer
    if(printNow) printTreeData(cells);                          // Print tree data
  }

//! Direct summation
  void direct(Bodies &ibodies, Bodies &jbodies) {
    Cells cells(2);                                             // Define a pair of cells to pass to P2P kernel
    C_iter Ci = cells.begin(), Cj = cells.begin()+1;            // First cell is target, second cell is source
    Ci->BODY = ibodies.begin();                                 // Iterator of first target body
    Ci->NDBODY = ibodies.size();                                // Number of target bodies
    Cj->BODY = jbodies.begin();                                 // Iterator of first source body
    Cj->NDBODY = jbodies.size();                                // Number of source bodies
    int prange = 0;                                             // Range of periodic images
    for (int i=0; i<IMAGES; i++) {                              // Loop over periodic image sublevels
      prange += int(powf(3,i));                                 //  Accumulate range of periodic images
    }                                                           // End loop over perioidc image sublevels
    for (int ix=-prange; ix<=prange; ix++) {                    // Loop over x periodic direction
      for (int iy=-prange; iy<=prange; iy++) {                  //  Loop over y periodic direction
        for (int iz=-prange; iz<=prange; iz++) {                //   Loop over z periodic direction
          Xperiodic[0] = ix * 2 * globalRadius;                 //    Coordinate shift for x periodic direction
          Xperiodic[1] = iy * 2 * globalRadius;                 //    Coordinate shift for y periodic direction
          Xperiodic[2] = iz * 2 * globalRadius;                 //    Coordinate shift for z periodic direction
          P2P(Ci,Cj,false);                                     //    Evaluate P2P kernel
        }                                                       //   End loop over z periodic direction
      }                                                         //  End loop over y periodic direction
    }                                                           // End loop over x periodic direction
  }

//! Normalize bodies after direct summation
  void normalize(Bodies &bodies) {
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {       // Loop over bodies
      B->TRG /= B->SRC;                                         //  Normalize by target charge
    }                                                           // End loop over bodies
  }
};

#endif
