#ifndef serialfmm_h
#define serialfmm_h
#include "treebuilder.h"

class SerialFMM : public TreeBuilder {
protected:
  real_t globalRadius;                                          //!< Radius of global root cell
  vec3 globalCenter;                                            //!< Center of global root cell
  vec3 globalXmin;                                              //!< Global Xmin for a given rank
  vec3 globalXmax;                                              //!< Global Xmax for a given rank

private:
//! Error optimization of Rcrit
  void setRcrit(C_iter C, real_t c) {
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

  //! Y := component-wise minimum of X and Y
  inline void vec3_min(vec3& X, vec3& Y) {
    if (X[0] < Y[0]) Y[0] = X[0];
    if (X[1] < Y[1]) Y[1] = X[1];
    if (X[2] < Y[2]) Y[2] = X[2];
  }

  //! Y := component-wise maximum of X and Y
  inline void vec3_max(vec3& X, vec3& Y) {
    if (X[0] > Y[0]) Y[0] = X[0];
    if (X[1] > Y[1]) Y[1] = X[1];
    if (X[2] > Y[2]) Y[2] = X[2];
  }

  //! get bounds of bodies in [B0,B1), in parallel
  //! return the pair of minimum and maximum
  vec3Pair getBoundsRec(B_iter B0, B_iter B1) {
    assert(B1 - B0 > 0);
    if (B1 - B0 < 1000) {
      vec3 xmin = B0->X, xmax = B0->X;
      for (B_iter B=B0; B!=B1; B++) {
        vec3_min(B->X, xmin);
        vec3_max(B->X, xmax);
      }
      return vec3Pair(xmin, xmax);
    } else {
      int nh = (B1 - B0) / 2;
      __init_tasks__;
      vec3Pair vt0, vt1;
      spawn_task1(vt0, vt0 = getBoundsRec(B0, B0 + nh));
      vt1 = getBoundsRec(B0 + nh, B1);
      __sync_tasks__;
      vec3_min(vt1.first, vt0.first);
      vec3_max(vt1.second, vt0.second);
      return vt0;
    }
  }

  void upwardPassRec1(C_iter C, C_iter C0) {
    __init_tasks__;                                             // Initialize tasks
    for (C_iter CC=C0+C->CHILD; CC!=C0+C->CHILD+C->NCHILD; CC++) {// Loop over child cells
      spawn_task0(upwardPassRec1(CC, C0));                      //  Recursive tree traversal
    }                                                           // End loop over child cells
    __sync_tasks__;                                             // Synchronize tasks
    real_t Rmax = 0;                                            //  Initialize Rmax
    setCenter(C);                                               //  Set center of cell to center of mass
    C->M = 0;                                                   //  Initialize multipole expansion coefficients
    C->L = 0;                                                   //  Initialize local expansion coefficients
    P2M(C,Rmax);                                                //  P2M kernel
    M2M(C,Rmax);                                                //  M2M kernel
  }

  void upwardPassRec2(C_iter C, C_iter C0, real_t c) {
    __init_tasks__;                                             // Initialize tasks
    for (C_iter CC=C0+C->CHILD; CC!=C0+C->CHILD+C->NCHILD; CC++) {// Loop over child cells
      spawn_task0(upwardPassRec2(CC, C0, c));                   //  Recursive tree traversal
    }                                                           // End loop over child cells
    __sync_tasks__;                                             // Synchronize tasks
    setRcrit(C, c);                                             // Error optimization of Rcrit
  }

  void downwardPassRec1(C_iter C, C_iter C0) const {
    L2L(C);                                                     // L2L kernel
    L2P(C);                                                     // L2P kernel
    __init_tasks__;                                             // Initialize tasks
    for (C_iter CC=C0+C->CHILD; CC!=C0+C->CHILD+C->NCHILD; CC++) {// Loop over child cells
      spawn_task0(downwardPassRec1(CC, C0));                    //  Recursive tree traversal
    }                                                           // End loop over chlid cells
    __sync_tasks__;                                             // Synchronize tasks
  }

public:
//! Set center and size of root cell
  void setBounds(Bodies &bodies) {
    startTimer("Set bounds");                                   // Start timer
    vec3Pair vt = getBoundsRec(bodies.begin(), bodies.end());   // Get Xmin and Xmax of domain
    localXmin = vt.first;                                       // Set local Xmin
    localXmax = vt.second;                                      // Set local Xmax
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
    growTreeRec(bodies);                                        // Grow tree from root
    linkTreeRec(bodies,cells);                                  // Form parent-child links in tree
  }

//! Upward pass (P2M, M2M)
  void upwardPass(Cells &cells) {
    startTimer("Upward pass");                                  // Start timer
    Ci0 = cells.begin();                                        // Set iterator of target root cell
    Cj0 = cells.begin();                                        // Set iterator of source root cell
    upwardPassRec1(Ci0, Ci0);
    real_t c = (1 - THETA) * (1 - THETA) / pow(THETA,P+2) / pow(std::abs(Ci0->M[0]),1.0/3);
    upwardPassRec2(Ci0, Ci0, c);
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
      spawn_task0(downwardPassRec1(CC, C0));                    //  Recursive tree traversal
    }                                                           // End loop over child cells
    __sync_tasks__;                                             // Synchronize tasks
    stopTimer("Downward pass",printNow);                        // Stop timer
    if(printNow) printTreeData(cells);                          // Print tree data
  }

//! Direct summation
  void direct(Bodies &ibodies, Bodies &jbodies) {
    Cells cells(2);                                             // Define a pair of cells to pass to P2P kernel
    C_iter Ci = cells.begin(), Cj = cells.begin()+1;            // First cell is target, second cell is source
    Ci->LEAF = ibodies.begin();                                 // Iterator of first target leaf
    Ci->NDLEAF = ibodies.size();                                // Number of target leafs
    Cj->LEAF = jbodies.begin();                                 // Iterator of first source leaf
    Cj->NDLEAF = jbodies.size();                                // Number of source leafs
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
