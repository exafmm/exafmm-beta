#ifndef evaluator_h
#define evaluator_h
#include "kernel.h"

class Evaluator : public Kernel {                               // Evaluator is the interface between tree and kernel
protected:
  C_iter      CI0;                                              // icells.begin()
  C_iter      CIB;                                              // icells begin per call
  C_iter      CIE;                                              // icells end per call
  C_iter      CJ0;                                              // jcells.begin()
  C_iter      CJB;                                              // jcells begin per call
  C_iter      CJE;                                              // jcells end per call
  Pairs       pairs;                                            // Stack of interacting cell pairs
  Lists       listM2L;                                          // M2L interaction list
  Lists       listM2P;                                          // M2P interaction list
  Lists       listP2P;                                          // P2P interaction list

  int         Iperiodic;                                        // Periodic image flag (using each bit for 27 images)
  const int   Icenter;                                          // Periodic image flag at center
  Maps        flagM2L;                                          // Flag indicating existance of periodic image for M2L
  Maps        flagM2P;                                          // Flag indicating existance of periodic image for M2P
  Maps        flagP2P;                                          // Flag indicating existance of periodic image for P2P

private:
  void treecode(C_iter Ci, C_iter Cj) {                         // Tree walk for treecode
    if( Ci->NCHILD == 0 && Cj->NCHILD == 0) {                   // If both cells are twigs
      if( Cj->NLEAF != 0 ) {                                    // If the twig has leafs
        tryP2P(Ci,Cj);                                          //  Try to evaluate P2P kernel
      } else {                                                  // If the twig has no leafs
#ifdef DEBUG
        std::cout << "Cj->ICELL=" << Cj->ICELL << " has no leaf. Doing M2P instead of P2P." << std::endl;
#endif
        listM2P[Ci-CI0].push_back(Cj);                          // Push source cell into M2P interaction list
      }                                                         // Endif for twigs with leafs
    } else if ( Ci->NCHILD != 0 ) {                             // If target is not twig
      for( int i=0; i<Ci->NCHILD; i++ ) {                       //  Loop over child cells of target
        tryM2P(CI0+Ci->CHILD[i],Cj);                            //   Try to evaluate M2P kernel
      }                                                         //  End loop over child cells of target
    } else {                                                    // If target is twig
      for( int i=0; i<Cj->NCHILD; i++ ) {                       //  Loop over child cells of source
        tryM2P(Ci,CJ0+Cj->CHILD[i]);                            //   Try to evaluate M2P kernel
      }                                                         //  End loop over child cells of source
    }                                                           // Endif for type of interaction
  }

  void FMM(C_iter Ci, C_iter Cj) {                              // Tree walk for FMM
    if( Ci->NCHILD == 0 && Cj->NCHILD == 0 ) {                  // If both cells are twigs
      if( Cj->NLEAF != 0 ) {                                    // If the twig has leafs
        tryP2P(Ci,Cj);                                          //  Try to evaluate P2P kernel
      } else {                                                  // If the twig has no leafs
//#ifdef DEBUG
        std::cout << "Cj->ICELL=" << Cj->ICELL << " has no leaf. Doing M2P instead of P2P." << std::endl;
//#endif
        listM2P[Ci-CI0].push_back(Cj);                          // Push source cell into M2P interaction list
      }                                                         // Endif for twigs with leafs
    } else if ( Cj->NCHILD == 0 || (Ci->NCHILD != 0 && Ci->R > Cj->R) ) {// If source is twig or target is larger
      for( int i=0; i<Ci->NCHILD; i++ ) {                       //  Loop over child cells of target
        tryM2L(CI0+Ci->CHILD[i],Cj);                            //   Try to evaluate M2L kernel
      }                                                         //  End loop over child cells of target
    } else {                                                    // If target is twig or source is larger
      for( int i=0; i<Cj->NCHILD; i++ ) {                       //  Loop over child cells of source
        tryM2L(Ci,CJ0+Cj->CHILD[i]);                            //   Try to evaluate M2L kernel
      }                                                         //  End loop over child cells of source
    }                                                           // Endif for type of interaction
  }

public:
  Evaluator() : Icenter(1 << 13) {}                             // Constructor
  ~Evaluator() {}                                               // Destructor

  void setKernel(std::string name) {                            // Set kernel name
    kernelName = name;                                          // Set class variable kernelName
  }

  void addM2L(C_iter Cj) {                                      // Add single list for kernel unit test
    listM2L.resize(1);                                          // Resize vector of M2L interation lists
    flagM2L.resize(1);                                          // Resize vector of M2L periodic image flags
    listM2L[0].push_back(Cj);                                   // Push single cell into list
    flagM2L[0][Cj] |= Icenter;                                  // Flip bit of periodic image flag
  }

  void addM2P(C_iter Cj) {                                      // Add single list for kernel unit test
    listM2P.resize(1);                                          // Resize vector of M2P interation lists
    flagM2P.resize(1);                                          // Resize vector of M2L periodic image flags
    listM2P[0].push_back(Cj);                                   // Push single cell into list
    flagM2P[0][Cj] |= Icenter;                                  // Flip bit of periodic image flag
  }

  int getPeriodicRange() {                                      // Get range of periodic images
    int prange = 0;                                             //  Range of periodic images
    for( int i=0; i!=IMAGES; ++i ) {                            //  Loop over periodic image sublevels
      prange += int(pow(3,i));                                  //   Accumulate range of periodic images
    }                                                           //  End loop over perioidc image sublevels
    return prange;                                              // Return range of periodic images
  }

  Bodies periodicBodies(Bodies &bodies) {                       // Create periodic images of bodies
    Bodies jbodies;                                             // Vector for periodic images of bodies
    int prange = getPeriodicRange();                            // Get range of periodic images
    for( int ix=-prange; ix<=prange; ++ix ) {                   // Loop over x periodic direction
      for( int iy=-prange; iy<=prange; ++iy ) {                 //  Loop over y periodic direction
        for( int iz=-prange; iz<=prange; ++iz ) {               //   Loop over z periodic direction
          for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {//    Loop over bodies
            Body body = *B;                                     //     Copy current body
            body.X[0] += ix * 2 * R0;                           //     Shift x position
            body.X[1] += iy * 2 * R0;                           //     Shift y position
            body.X[2] += iz * 2 * R0;                           //     Shift z position
            jbodies.push_back(body);                            //     Push shifted body into jbodies
          }                                                     //    End loop over bodies
        }                                                       //   End loop over z periodic direction
      }                                                         //  End loop over y periodic direction
    }                                                           // End loop over x periodic direction
    return jbodies;                                             // Return vector for periodic images of bodies
  }

  void traverse(Cells &cells, Cells &jcells, int method) {      // Traverse tree to get interaction list
    C_iter root = cells.end() - 1;                              // Iterator for root target cell
    C_iter jroot = jcells.end() - 1;                            // Iterator for root source cell
    if( IMAGES != 0 ) {                                         // If periodic boundary condition
      jroot = jcells.end() - 1 - 26 * (IMAGES - 1);             //  The root is not at the end
    }                                                           // Endif for periodic boundary condition
    CI0 = cells.begin();                                        // Set begin iterator for target cells
    CJ0 = jcells.begin();                                       // Set begin iterator for source cells
    listM2L.resize(cells.size());                               // Resize M2L interaction list
    listM2P.resize(cells.size());                               // Resize M2P interaction list
    listP2P.resize(cells.size());                               // Resize P2P interaction list
    flagM2L.resize(cells.size());                               // Resize M2L periodic image flag
    flagM2P.resize(cells.size());                               // Resize M2P periodic image flag
    flagP2P.resize(cells.size());                               // Resize P2P periodic image flag
    if( IMAGES == 0 ) {                                         // If free boundary condition
      Iperiodic = Icenter;                                      //  Set periodic image flag to center
      Xperiodic = 0;                                            //  Set periodic coordinate offset
      Pair pair(root,jroot);                                    //  Form pair of root cells
      pairs.push(pair);                                         //  Push pair to stack
      while( !pairs.empty() ) {                                 //  While interaction stack is not empty
        pair = pairs.top();                                     //   Get interaction pair from top of stack
        pairs.pop();                                            //   Pop interaction stack
        switch (method) {                                       //   Swtich between methods
        case 0 : treecode(pair.first,pair.second); break;       //    0 : treecode
        case 1 : FMM(pair.first,pair.second);      break;       //    1 : FMM
        }                                                       //   End switch between methods
      }                                                         //  End while loop for interaction stack
    } else {                                                    // If periodic boundary condition
      int I = 0;                                                //  Initialize index of periodic image
      for( int ix=-1; ix<=1; ++ix ) {                           //  Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //   Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz, ++I ) {                  //    Loop over z periodic direction
            Iperiodic = 1 << I;                                 //     Set periodic image flag
            Xperiodic[0] = ix * 2 * R0;                         //     Coordinate offset for x periodic direction
            Xperiodic[1] = iy * 2 * R0;                         //     Coordinate offset for y periodic direction
            Xperiodic[2] = iz * 2 * R0;                         //     Coordinate offset for z periodic direction
            Pair pair(root,jroot);                              //     Form pair of root cells
            pairs.push(pair);                                   //     Push pair to stack
            while( !pairs.empty() ) {                           //     While interaction stack is not empty
              pair = pairs.top();                               //      Get interaction pair from top of stack
              pairs.pop();                                      //      Pop interaction stack
              switch (method) {                                 //      Swtich between methods
              case 0 : treecode(pair.first,pair.second); break; //      0 : treecode
              case 1 : FMM(pair.first,pair.second);      break; //      1 : FMM
              }                                                 //      End switch between methods
            }                                                   //     End while loop for interaction stack
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      for( C_iter Ci=cells.begin(); Ci!=cells.end(); ++Ci ) {   //  Loop over target cells
        listM2L[Ci-CI0].sort();                                 //  Sort interaction list
        listM2L[Ci-CI0].unique();                               //  Eliminate duplicate periodic entries
        listM2P[Ci-CI0].sort();                                 //  Sort interaction list
        listM2P[Ci-CI0].unique();                               //  Eliminate duplicate periodic entries
        listP2P[Ci-CI0].sort();                                 //  Sort interaction list
        listP2P[Ci-CI0].unique();                               //  Eliminate duplicate periodic entries
      }                                                         //  End loop over target cells
    }                                                           // Endif for periodic boundary condition
  }

  void upwardPeriodic(Cells &jcells) {                          // Upward phase for periodic cells
    Cells pccells, pjcells;                                     // Periodic jcells for M2L/M2P & M2M
    pccells.push_back(jcells.back());                           // Root cell is first periodic cell
    for( int level=0; level<IMAGES-1; ++level ) {               // Loop over sublevels of tree
      C_iter C = pccells.end() - 1;                             //  Set previous periodic cell as source
      for( int ix=-1; ix<=1; ++ix ) {                           //  Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //   Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz ) {                       //    Loop over z periodic direction
            Cell cell;                                          //     New periodic jcell for M2M
            cell.X[0] = C->X[0] + ix * 2 * C->R;                //     Set new x coordinate for periodic image
            cell.X[1] = C->X[1] + iy * 2 * C->R;                //     Set new y cooridnate for periodic image
            cell.X[2] = C->X[2] + iz * 2 * C->R;                //     Set new z coordinate for periodic image
            cell.M = C->M;                                      //     Copy multipoles to new periodic image
            pjcells.push_back(cell);                            //     Push cell into periodic jcell vector
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      Cell cell;                                                //  New periodic cell at next sublevel
      cell.X = C->X;                                            //  This is the center cell
      cell.R = 3 * C->R;                                        //  The cell size increase three times
      pccells.push_back(cell);                                  //  Push cell into periodic cell vector
      CI = pccells.end() - 1;                                   //  Set current cell as target for M2M
      while( !pjcells.empty() ) {                               //  While there are periodic jcells remaining
        CJ = pjcells.end() - 1;                                 //   Set current jcell as source for M2M
        selectM2M_CPU();                                        //   Select M2M_CPU kernel
        pjcells.pop_back();                                     //   Pop last element from periodic jcell vector
      }                                                         //  End while for remaining periodic jcells
      for( int ix=-1; ix<=1; ++ix ) {                           //  Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //   Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz ) {                       //    Loop over z periodic direction
            if( ix != 0 || iy != 0 || iz != 0 ) {               //     If periodic cell is not at center
              cell.X[0]  = CI->X[0] + ix * 2 * CI->R;           //      Set new x coordinate for periodic image
              cell.X[1]  = CI->X[1] + iy * 2 * CI->R;           //      Set new y cooridnate for periodic image
              cell.X[2]  = CI->X[2] + iz * 2 * CI->R;           //      Set new z coordinate for periodic image
              cell.M     = CI->M;                               //      Copy multipoles to new periodic image
              cell.NLEAF = cell.NCHILD = 0;                     //      Initialize NLEAF & NCHILD
              jcells.push_back(cell);                           //      Push cell into periodic jcell vector
            }                                                   //     Endif for periodic center cell
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
    }                                                           // End loop over sublevels of tree
  }

  void initialize() {                                           // Initialize GPU
    if( kernelName == "Laplace" ) {                             // If Laplace kernel
      LaplaceInit();                                            //  Initialize GPU
    } else if ( kernelName == "BiotSavart" ) {                  // If Biot Savart kernel
      BiotSavartInit();                                         //  Initialize GPU
    } else if ( kernelName == "Stretching" ) {                  // If Stretching kernel
      StretchingInit();                                         //  Initialize GPU
    } else if ( kernelName == "Gaussian" ) {                    // If Gaussian kernel
      GaussianInit();                                           //  Initialize GPU
    } else if ( kernelName == "CoulombVdW" ) {                  // If CoulombVdW kernel
      CoulombVdWInit();                                         //  Initialize GPU
    } else {                                                    // If kernel is none of the above
      if(MPIRANK == 0) std::cout << "Invalid kernel type in initialize" << std::endl;// Invalid kernel type
      abort();                                                  //  Abort execution
    }                                                           // Endif for kernel type
  }

  void finalize() {                                             // Finalize GPU
    if( kernelName == "Laplace" ) {                             // If Laplace kernel
      LaplaceFinal();                                           //  Finalize GPU
    } else if ( kernelName == "BiotSavart" ) {                  // If Biot Savart kernel
      BiotSavartFinal();                                        //  Finalize GPU
    } else if ( kernelName == "Stretching" ) {                  // If Stretching kernel
      StretchingFinal();                                        //  Finalize GPU
    } else if ( kernelName == "Gaussian" ) {                    // If Gaussian kernel
      GaussianFinal();                                          //  Finalize GPU
    } else if ( kernelName == "CoulombVdW" ) {                  // If CoulombVdW kernel
      CoulombVdWFinal();                                        //  Finalize GPU
    } else {                                                    // If kernel is none of the above
      if(MPIRANK == 0) std::cout << "Invalid kernel type in finalize" << std::endl;// Invalid kernel type
      abort();                                                  //  Abort execution
    }                                                           // Endif for kernel type
  }

  void setSourceBody();                                         // Set source buffer for bodies
  void setSourceCell(bool isM);                                 // Set source buffer for cells
  void setTargetBody(Lists lists, Maps flags);                  // Set target buffer for bodies
  void setTargetCell(Lists lists, Maps flags);                  // Set target buffer for cells
  void getTargetBody(Lists &lists);                             // Get body values from target buffer
  void getTargetCell(Lists &lists, bool isM);                   // Get cell values from target buffer
  void clearBuffers();                                          // Clear GPU buffers

  void tryP2P(C_iter Ci, C_iter Cj);                            // Interface for P2P kernel
  void tryM2L(C_iter Ci, C_iter Cj);                            // Interface for M2L kernel
  void tryM2P(C_iter Ci, C_iter Cj);                            // Interface for M2P kernel
  void traversePeriodic(Cells &cells, Cells &jcells, int method);// Traverse tree for periodic cells
  void evalP2P(Bodies &ibodies, Bodies &jbodies, bool onCPU=false);// Evaluate P2P kernel (all pairs)
  void evalP2M(Cells &twigs);                                   // Evaluate P2M kernel
  void evalM2M(Cells &cells);                                   // Evaluate M2M kernel
  void evalM2L(Cells &cells, bool kernel=false);                // Evaluate M2L kernel
  void evalM2P(Cells &cells, bool kernel=false);                // Evaluate M2P kernel
  void evalP2P(Cells &cells, bool kernel=false);                // Evaluate P2P kernel (near field)
  void evalL2L(Cells &cells);                                   // Evaluate L2L kernel
  void evalL2P(Cells &cells);                                   // Evaluate L2P kernel
};
#if cpu
#include "../kernel/cpuEvaluator.cxx"
#elif gpu
#include "../kernel/gpuEvaluator.cxx"
#endif

#endif
