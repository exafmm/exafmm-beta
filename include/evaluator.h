#ifndef evaluator_h
#define evaluator_h
#include "periodic.h"

class Evaluator : public Periodic {                             // Evaluator is the interface between tree and kernel
private:
  C_iter    CI0;                                                // icells.begin()
  C_iter    CJ0;                                                // jcells.begin()
  Pairs     pairs;                                              // Stack of interacting cell pairs
  Lists     listM2L;                                            // M2L interaction list
  Lists     listM2P;                                            // M2P interaction list
  Lists     listP2P;                                            // P2P interaction list

private:
  void tryM2L(C_iter Ci, C_iter Cj) {                           // Interface for M2L kernel
    vect dist = Ci->X - Cj->X - Xperiodic;                      // Distance vector between cells
    real R = std::sqrt(norm(dist));                             // Distance between cells
    if( Ci->R + Cj->R > THETA*R ) {                             // If cell is too large
      Pair pair(Ci,Cj);                                         //  Form pair of interacting cells
      pairs.push(pair);                                         //  Push interacting pair into stack
    } else {                                                    // If cell is small enough
      listM2L[Ci-CI0].push_back(Cj);                            // Push source cell into M2L interaction list
      flagM2L[Ci-CI0][Cj] |= Iperiodic;                         // Flip bit of periodic image flag
    }                                                           // Endif for interaction
  }

  void tryM2P(C_iter Ci, C_iter Cj) {                           // Interface for M2P kernel
    vect dist = Ci->X - Cj->X - Xperiodic;                      // Distance vector between cells
    real R = std::sqrt(norm(dist));                             // Distance between cells
    if( Ci->NCHILD != 0 || Ci->R + Cj->R > THETA*R ) {          // If target is not twig or cell is too large
      Pair pair(Ci,Cj);                                         //  Form pair of interacting cells
      pairs.push(pair);                                         //  Push interacting pair into stack
    } else {                                                    // If target is twig and cell is small enough
      listM2P[Ci-CI0].push_back(Cj);                            // Push source cell into M2P interaction list
      flagM2P[Ci-CI0][Cj] |= Iperiodic;                         // Flip bit of periodic image flag
    }                                                           // Endif for interaction
  }

  void treecode(C_iter Ci, C_iter Cj) {                         // Tree walk for treecode
    if( Ci->NCHILD == 0 && Cj->NCHILD == 0) {                   // If both cells are twigs
      if( Cj->NLEAF != 0 ) {                                    // If the twig has leafs
        listP2P[Ci-CI0].push_back(Cj);                          // Push source cell into P2P interaction list
        flagP2P[Ci-CI0][Cj] |= Iperiodic;                       // Flip bit of periodic image flag
      } else {                                                  // If the twig has no leafs
#ifdef DEBUG
        std::cout << "Cj->I=" << Cj->I << " has no leaf. Doing M2P instead of P2P." << std::endl;
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
        listP2P[Ci-CI0].push_back(Cj);                          // Push source cell into P2P interaction list
        flagP2P[Ci-CI0][Cj] |= Iperiodic;                       // Flip bit of periodic image flag
      } else {                                                  // If the twig has no leafs
#ifdef DEBUG
        std::cout << "Cj->I=" << Cj->I << " has no leaf. Doing M2P instead of P2P." << std::endl;
#endif
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

  void setSourceBody();                                         // Set source buffer for bodies
  void setSourceCell(bool isM);                                 // Set source buffer for cells
  void setTargetBody(Cells &cells, Lists lists, Maps flags);    // Set target buffer for bodies
  void setTargetCell(Cells &cells, Lists lists, Maps flags);    // Set target buffer for cells
  void getTargetBody(Cells &cells, Lists &lists);               // Get body values from target buffer
  void getTargetCell(Cells &cells, Lists &lists, bool isM);     // Get cell values from target buffer
  void clearBuffers();                                          // Clear GPU buffers

  void evalP2P(Bodies &ibodies, Bodies &jbodies);               // Evaluate P2P kernel (all pairs)
  void evalP2M(Cells &twigs);                                   // Evaluate P2M kernel
  void evalM2M(Cells &cells);                                   // Evaluate M2M kernel
  void evalM2L(Cells &cells);                                   // Evaluate M2L kernel
  void evalM2P(Cells &cells);                                   // Evaluate M2P kernel
  void evalP2P(Cells &cells);                                   // Evaluate P2P kernel (near field)
  void evalL2L(Cells &cells);                                   // Evaluate L2L kernel
  void evalL2P(Cells &cells);                                   // Evaluate L2P kernel

  void traverse(Cells &cells, Cells &jcells, int method) {      // Traverse tree to get interaction list
    C_iter root = cells.end()-1;                                // Iterator for root cell
    C_iter jroot = jcells.end()-1;                              // Iterator for root cell
    CI0 = cells.begin();                                        // Set begin iterator for icells
    CJ0 = jcells.begin();                                       // Set begin iterator for jcells
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
          case 0 : treecode(pair.first,pair.second); break;     //    0 : treecode
          case 1 : FMM(pair.first,pair.second);      break;     //    1 : FMM
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
                case 0 : treecode(pair.first,pair.second); break;//      0 : treecode
                case 1 : FMM(pair.first,pair.second);      break;//      1 : FMM
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

};

#endif
