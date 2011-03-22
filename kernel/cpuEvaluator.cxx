#define EVALUATOR
#include "evaluator.h"
#undef EVALUATOR

void Evaluator::evalP2P(Bodies &ibodies, Bodies &jbodies, bool onCPU) {// Evaluate P2P
  BI0 = ibodies.begin();                                        // Set target bodies begin iterator
  BIN = ibodies.end();                                          // Set target bodies end iterator
  BJ0 = jbodies.begin();                                        // Set source bodies begin iterator
  BJN = jbodies.end();                                          // Set source bodies end iterator
  Xperiodic = 0 * onCPU;                                        // Set periodic coordinate offset (onCPU is dummy)
  P2P_CPU();                                                    // Evaluate P2P kernel
}

void Evaluator::evalP2M(Cells &cells) {                         // Evaluate P2M
  startTimer("evalP2M      ");                                  // Start timer
  for( CJ=cells.begin(); CJ!=cells.end(); ++CJ ) {              // Loop over cells
    CJ->M = CJ->L = 0;                                          //  Initialize multipole & local coefficients
    P2M();                                                      //  Evaluate P2M kernel
  }                                                             // End loop over cells
  stopTimer("evalP2M      ");                                   // Stop timer
}

void Evaluator::evalM2M(Cells &cells) {                         // Evaluate M2M
  startTimer("evalM2M      ");                                  // Start timer
  CJ0 = cells.begin();                                          // Set begin iterator
  for( CJ=cells.begin(); CJ!=cells.end()-1; ++CJ ) {            // Loop over cells bottomup (except root cell)
    CI = CJ0 + CJ->PARENT;                                      //  Set target cell iterator
    M2M_CPU();                                                  //  Evaluate M2M kernel
  }                                                             // End loop over cells
  stopTimer("evalM2M      ");                                   // Stop timer
}

void Evaluator::evalM2L(Cells &cells) {                         // Evaluate M2L
  startTimer("evalM2L      ");                                  // Start timer
  CI0 = cells.begin();                                          // Set begin iterator
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {              // Loop over cells
    while( !listM2L[CI-CI0].empty() ) {                         //  While M2L interaction list is not empty
      CJ = listM2L[CI-CI0].back();                              //   Set source cell iterator
      Iperiodic = flagM2L[CI-CI0][CJ];                          //   Set periodic image flag
      int I = 0;                                                //   Initialize index of periodic image
      for( int ix=-1; ix<=1; ++ix ) {                           //   Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //    Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz, ++I ) {                  //     Loop over z periodic direction
            if( Iperiodic & (1 << I) ) {                        //      If periodic flag is on
              Xperiodic[0] = ix * 2 * R0;                       //       Coordinate offset for x periodic direction
              Xperiodic[1] = iy * 2 * R0;                       //       Coordinate offset for y periodic direction
              Xperiodic[2] = iz * 2 * R0;                       //       Coordinate offset for z periodic direction
              M2L();                                            //       Evaluate M2L kernel
            }                                                   //      Endif for periodic flag
          }                                                     //     End loop over x periodic direction
        }                                                       //    End loop over y periodic direction
      }                                                         //   End loop over z periodic direction
      listM2L[CI-CI0].pop_back();                               //   Pop last element from M2L interaction list
    }                                                           //  End while for M2L interaction list
  }                                                             // End loop over cells topdown
  listM2L.clear();                                              // Clear interaction lists
  flagM2L.clear();                                              // Clear periodic image flags
  stopTimer("evalM2L      ");                                   // Stop timer
}

void Evaluator::evalM2P(Cells &cells) {                         // Evaluate M2P
  startTimer("evalM2P      ");                                  // Start timer
  CI0 = cells.begin();                                          // Set begin iterator
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {              // Loop over cells
    while( !listM2P[CI-CI0].empty() ) {                         //  While M2P interaction list is not empty
      CJ = listM2P[CI-CI0].back();                              //   Set source cell iterator
      Iperiodic = flagM2P[CI-CI0][CJ];                          //   Set periodic image flag
      int I = 0;                                                //   Initialize index of periodic image
      for( int ix=-1; ix<=1; ++ix ) {                           //   Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //    Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz, ++I ) {                  //     Loop over z periodic direction
            if( Iperiodic & (1 << I) ) {                        //      If periodic flag is on
              Xperiodic[0] = ix * 2 * R0;                       //       Coordinate offset for x periodic direction
              Xperiodic[1] = iy * 2 * R0;                       //       Coordinate offset for y periodic direction
              Xperiodic[2] = iz * 2 * R0;                       //       Coordinate offset for z periodic direction
              M2P();                                            //       Evaluate M2P kernel
            }                                                   //      Endif for periodic flag
          }                                                     //     End loop over x periodic direction
        }                                                       //    End loop over y periodic direction
      }                                                         //   End loop over z periodic direction
      listM2P[CI-CI0].pop_back();                               //   Pop last element from M2P interaction list
    }                                                           //  End while for M2P interaction list
  }                                                             // End loop over cells topdown
  listM2P.clear();                                              // Clear interaction lists
  flagM2P.clear();                                              // Clear periodic image flags
  stopTimer("evalM2P      ");                                   // Stop timer
}

void Evaluator::evalP2P(Cells &cells) {                         // Evaluate P2P
  startTimer("evalP2P      ");                                  // Start timer
  CI0 = cells.begin();                                          // Set begin iterator
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {              // Loop over cells
    BI0 = CI->LEAF;                                             //  Set target bodies begin iterator
    BIN = CI->LEAF + CI->NLEAF;                                 //  Set target bodies end iterator
    while( !listP2P[CI-CI0].empty() ) {                         //  While M2P interaction list is not empty
      CJ = listP2P[CI-CI0].back();                              //   Set source cell iterator
      BJ0 = CJ->LEAF;                                           //   Set source bodies begin iterator
      BJN = CJ->LEAF + CJ->NLEAF;                               //   Set source bodies end iterator
      Iperiodic = flagP2P[CI-CI0][CJ];                          //   Set periodic image flag
      int I = 0;                                                //   Initialize index of periodic image
      for( int ix=-1; ix<=1; ++ix ) {                           //   Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //    Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz, ++I ) {                  //     Loop over z periodic direction
            if( Iperiodic & (1 << I) ) {                        //      If periodic flag is on
              Xperiodic[0] = ix * 2 * R0;                       //       Coordinate offset for x periodic direction
              Xperiodic[1] = iy * 2 * R0;                       //       Coordinate offset for y periodic direction
              Xperiodic[2] = iz * 2 * R0;                       //       Coordinate offset for z periodic direction
              P2P_CPU();                                        //       Evaluate P2P kernel
            }                                                   //      Endif for periodic flag
          }                                                     //     End loop over x periodic direction
        }                                                       //    End loop over y periodic direction
      }                                                         //   End loop over z periodic direction
      listP2P[CI-CI0].pop_back();                               //   Pop last element from M2P interaction list
    }                                                           //  End while for M2P interaction list
  }                                                             // End loop over cells topdown
  listP2P.clear();                                              // Clear interaction lists
  flagP2P.clear();                                              // Clear periodic image flags
  stopTimer("evalP2P      ");                                   // Stop timer
}

void Evaluator::evalL2L(Cells &cells) {                         // Evaluate L2L
  startTimer("evalL2L      ");                                  // Start timer
  CI0 = cells.begin();                                          // Set begin iterator
  for( CI=cells.end()-2; CI!=cells.begin()-1; --CI ) {          // Loop over cells topdown (except root cell)
    CJ = CI0 + CI->PARENT;                                      //  Set source cell iterator
    L2L();                                                      //  Evaluate L2L kernel
  }                                                             // End loop over cells topdown
  stopTimer("evalL2L      ");                                   // Stop timer
}

void Evaluator::evalL2P(Cells &cells) {                         // Evaluate L2P
  startTimer("evalL2P      ");                                  // Start timer
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {              // Loop over cells
    if( CI->NCHILD == 0 ) L2P();                                //  If cell is a twig evaluate L2P kernel
  }                                                             // End loop over cells topdown
  stopTimer("evalL2P      ");                                   // Stop timer
}
