/*
Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
template<Equation equation>
void Evaluator<equation>::testMACP2P(C_iter CI, C_iter CJ) {    // Test multipole acceptance criteria for P2P kernel
  P2P_CPU(CI,CJ);                                               // Perform P2P_CPU kernel
  NP2P++;                                                       // Count P2P kernel execution
}

template<Equation equation>
void Evaluator<equation>::testMACM2L(C_iter CI, C_iter CJ) {    // Test multipole acceptance criteria for M2L kernel
  vect dist = CI->X - CJ->X - Xperiodic;                        // Distance vector between cells
  real R = std::sqrt(norm(dist));                               // Distance between cells
  if( CI->R + CJ->R > THETA*R ) {                               // If cell is too large
    Pair pair(CI,CJ);                                           //  Form pair of interacting cells
    pairs.push(pair);                                           //  Push interacting pair into stack
  } else {                                                      // If cell is small enough
    M2L(CI,CJ);                                                 //  Perform M2L kernel
    NM2L++;                                                     //  Count M2L kernel execution
  }                                                             // Endif for interaction
}

template<Equation equation>
void Evaluator<equation>::testMACM2P(C_iter CI, C_iter CJ) {    // Test multipole acceptance criteria for M2P kernel
  vect dist = CI->X - CJ->X - Xperiodic;                        // Distance vector between cells
  real R = std::sqrt(norm(dist));                               // Distance between cells
  if( CI->NCHILD != 0 || CI->R + CJ->R > THETA*R ) {            // If target is not twig or cell is too large
    Pair pair(CI,CJ);                                           //  Form pair of interacting cells
    pairs.push(pair);                                           //  Push interacting pair into stack
  } else {                                                      // If target is twig and cell is small enough
    M2P(CI,CJ);                                                 //  Perform M2P kernel
    NM2P++;                                                     //  Count M2P kernel execution
  }                                                             // Endif for interaction
}

template<Equation equation>
void Evaluator<equation>::timeKernels() {                       // Time all kernels for auto-tuning
  Bodies ibodies(1000), jbodies(1000);                          // Artificial bodies
  for( B_iter Bi=ibodies.begin(),Bj=jbodies.begin(); Bi!=ibodies.end(); ++Bi, ++Bj ) {// Loop over artificial bodies
    Bi->X = 0;                                                  //  Set coordinates of target body
    Bj->X = 1;                                                  //  Set coordinates of source body
  }                                                             // End loop over artificial bodies
  Cells cells;                                                  // Artificial cells
  cells.resize(2);                                              // Two artificial cells
  C_iter CI = cells.begin(), CJ = cells.begin()+1;              // Artificial target & source cell
  CI->X = 0;                                                    // Set coordinates of target cell
  CI->NDLEAF = 10;                                              // Number of leafs in target cell
  CI->LEAF = ibodies.begin();                                   // Leaf iterator in target cell
  CJ->X = 1;                                                    // Set coordinates of source cell
  CJ->NDLEAF = 1000;                                            // Number of leafs in source cell
  CJ->LEAF = jbodies.begin();                                   // Leaf iterator in source cell
  startTimer("P2P kernel   ");                                  // Start timer
  for( int i=0; i!=1; ++i ) P2P_CPU(CI,CJ);                     // Perform P2P_CPU kernel
  timeP2P = stopTimer("P2P kernel   ") / 10000;                 // Stop timer
  startTimer("M2L kernel   ");                                  // Start timer
  for( int i=0; i!=1000; ++i ) M2L(CI,CJ);                      // Perform M2L kernel
  timeM2L = stopTimer("M2L kernel   ") / 1000;                  // Stop timer
  startTimer("M2P kernel   ");                                  // Start timer
  for( int i=0; i!=100; ++i ) M2P(CI,CJ);                       // Perform M2P kernel
  timeM2P = stopTimer("M2P kernel   ") / 1000;                  // Stop timer
}

template<Equation equation>
void Evaluator<equation>::traversePeriodic(Cells &cells, Cells &jcells, int method) {// Traverse tree for periodic cells
  Xperiodic = 0;                                                // Set periodic coordinate offset
  C_iter CJ = jcells.end()-1;                                   // Initialize iterator for periodic source cell
  for( int level=0; level<IMAGES-1; ++level ) {                 // Loop over sublevels of tree
    for( int I=0; I!=26; ++I, --CJ ) {                          //  Loop over periodic images (exclude center)
      switch (method) {                                         //   Switch between method
      case 0 :                                                  //   0 : treecode
        for( C_iter CI=cells.begin(); CI!=cells.end(); ++CI ) { //   Loop over cells
          if( CI->NCHILD == 0 ) {                               //     If cell is twig
            M2P(CI,CJ);                                         //      Perform M2P kernel
          }                                                     //     Endif for twig
        }                                                       //    End loop over cells
        break;                                                  //    Terminate this case
      case 1 :                                                  //   1 : FMM
      case 2 :                                                  //   2 : hybrid
        C_iter CI = cells.end() - 1;                            //    Set root cell as target iterator
        M2L(CI,CJ);                                             //    Perform M2P kernel
        break;                                                  //    Terminate this case
      }                                                         //   End switch between methods
    }                                                           //  End loop over x periodic direction
  }                                                             // End loop over sublevels of tree
}

template<Equation equation>
void Evaluator<equation>::evalP2P(Bodies &ibodies, Bodies &jbodies, bool onCPU) {// Evaluate P2P
  Xperiodic = 0 * onCPU;                                        // Set periodic coordinate offset (onCPU is dummy)
  Cells cells;                                                  // Cells to put target and source bodies
  cells.resize(2);                                              // Resize cells to put target and source bodies
  cells[0].LEAF = ibodies.begin();                              // Iterator of first target leaf
  cells[0].NDLEAF = ibodies.size();                             // Number of target leafs
  cells[1].LEAF = jbodies.begin();                              // Iterator of first source leaf
  cells[1].NDLEAF = jbodies.size();                             // Number of source leafs
  C_iter CI = cells.begin(), CJ = cells.begin()+1;              // Iterator of target and source cells
  P2P_CPU(CI,CJ);                                               // Perform P2P_CPU kernel
}

template<Equation equation>
void Evaluator<equation>::evalP2M(Cells &cells) {               // Evaluate P2M
  startTimer("evalP2M      ");                                  // Start timer
  for( C_iter CI=cells.begin(); CI!=cells.end(); ++CI ) {       // Loop over cells
    CI->M = 0;                                                  //  Initialize multipole coefficients
    CI->L = 0;                                                  //  Initialize local coefficients
    if( CI->NCHILD == 0 ) {                                     //  If cell is a twig
      P2M(CI);                                                  //   Perform P2M kernel
    }                                                           //  Endif for twig
  }                                                             // End loop over cells
  stopTimer("evalP2M      ");                                   // Stop timer
}

template<Equation equation>
void Evaluator<equation>::evalM2M(Cells &cells) {               // Evaluate M2M
  startTimer("evalM2M      ");                                  // Start timer
  CJ0 = cells.begin();                                          // Set begin iterator
  for( C_iter CJ=cells.begin(); CJ!=cells.end()-1; ++CJ ) {     // Loop over cells bottomup (except root cell)
    C_iter CI = CJ0 + CJ->PARENT;                               //  Set target cell iterator
    M2M_CPU(CI,CJ);                                             //  Perform M2M_CPU kernel
  }                                                             // End loop over cells
  stopTimer("evalM2M      ");                                   // Stop timer
}

template<Equation equation>
void Evaluator<equation>::evalM2L(Cells &cells, bool kernel) {  // Evaluate M2L
  if( kernel ) {                                                // If this is a kernel unit test
    startTimer("evalM2L      ");                                // Start timer
    CI0 = cells.begin();                                        // Set begin iterator
    for( C_iter CI=cells.begin(); CI!=cells.end(); ++CI ) {     // Loop over cells
      while( !listM2L[CI-CI0].empty() ) {                       //  While M2L interaction list is not empty
        C_iter CJ = listM2L[CI-CI0].back();                     //   Set source cell iterator
        Iperiodic = flagM2L[CI-CI0][CJ];                        //   Set periodic image flag
        int I = 0;                                              //   Initialize index of periodic image
        for( int ix=-1; ix<=1; ++ix ) {                         //   Loop over x periodic direction
          for( int iy=-1; iy<=1; ++iy ) {                       //    Loop over y periodic direction
            for( int iz=-1; iz<=1; ++iz, ++I ) {                //     Loop over z periodic direction
              if( Iperiodic & (1 << I) ) {                      //      If periodic flag is on
                Xperiodic[0] = ix * 2 * R0;                     //       Coordinate offset for x periodic direction
                Xperiodic[1] = iy * 2 * R0;                     //       Coordinate offset for y periodic direction
                Xperiodic[2] = iz * 2 * R0;                     //       Coordinate offset for z periodic direction
                M2L(CI,CJ);                                     //       Perform M2L kernel
              }                                                 //      Endif for periodic flag
            }                                                   //     End loop over x periodic direction
          }                                                     //    End loop over y periodic direction
        }                                                       //   End loop over z periodic direction
        listM2L[CI-CI0].pop_back();                             //   Pop last element from M2L interaction list
      }                                                         //  End while for M2L interaction list
    }                                                           // End loop over cells topdown
    listM2L.clear();                                            // Clear interaction lists
    flagM2L.clear();                                            // Clear periodic image flags
    stopTimer("evalM2L      ");                                 // Stop timer
  }                                                             // End if for kernel unit test
}

template<Equation equation>
void Evaluator<equation>::evalM2P(Cells &cells, bool kernel) {  // Evaluate M2P
  if( kernel ) {                                                // If this is a kernel unit test
    startTimer("evalM2P      ");                                // Start timer
    CI0 = cells.begin();                                        // Set begin iterator
    for( C_iter CI=cells.begin(); CI!=cells.end(); ++CI ) {     // Loop over cells
      while( !listM2P[CI-CI0].empty() ) {                       //  While M2P interaction list is not empty
        C_iter CJ = listM2P[CI-CI0].back();                     //   Set source cell iterator
        Iperiodic = flagM2P[CI-CI0][CJ];                        //   Set periodic image flag
        int I = 0;                                              //   Initialize index of periodic image
        for( int ix=-1; ix<=1; ++ix ) {                         //   Loop over x periodic direction
          for( int iy=-1; iy<=1; ++iy ) {                       //    Loop over y periodic direction
            for( int iz=-1; iz<=1; ++iz, ++I ) {                //     Loop over z periodic direction
              if( Iperiodic & (1 << I) ) {                      //      If periodic flag is on
                Xperiodic[0] = ix * 2 * R0;                     //       Coordinate offset for x periodic direction
                Xperiodic[1] = iy * 2 * R0;                     //       Coordinate offset for y periodic direction
                Xperiodic[2] = iz * 2 * R0;                     //       Coordinate offset for z periodic direction
                M2P(CI,CJ);                                     //       Perform M2P kernel
              }                                                 //      Endif for periodic flag
            }                                                   //     End loop over x periodic direction
          }                                                     //    End loop over y periodic direction
        }                                                       //   End loop over z periodic direction
        listM2P[CI-CI0].pop_back();                             //   Pop last element from M2P interaction list
      }                                                         //  End while for M2P interaction list
    }                                                           // End loop over cells topdown
    listM2P.clear();                                            // Clear interaction lists
    flagM2P.clear();                                            // Clear periodic image flags
    stopTimer("evalM2P      ");                                 // Stop timer
  }                                                             // End if for kernel unit test
}

template<Equation equation>
void Evaluator<equation>::evalP2P(Cells &cells, bool kernel) {  // Evaluate P2P
  if( kernel ) {                                                // If this is a kernel unit test
    startTimer("evalP2P      ");                                // Start timer
    CI0 = cells.begin();                                        // Set begin iterator
    for( C_iter CI=cells.begin(); CI!=cells.end(); ++CI ) {     // Loop over cells
      while( !listP2P[CI-CI0].empty() ) {                       //  While M2P interaction list is not empty
        C_iter CJ = listP2P[CI-CI0].back();                     //   Set source cell iterator
        Iperiodic = flagP2P[CI-CI0][CJ];                        //   Set periodic image flag
        int I = 0;                                              //   Initialize index of periodic image
        for( int ix=-1; ix<=1; ++ix ) {                         //   Loop over x periodic direction
          for( int iy=-1; iy<=1; ++iy ) {                       //    Loop over y periodic direction
            for( int iz=-1; iz<=1; ++iz, ++I ) {                //     Loop over z periodic direction
              if( Iperiodic & (1 << I) ) {                      //      If periodic flag is on
                Xperiodic[0] = ix * 2 * R0;                     //       Coordinate offset for x periodic direction
                Xperiodic[1] = iy * 2 * R0;                     //       Coordinate offset for y periodic direction
                Xperiodic[2] = iz * 2 * R0;                     //       Coordinate offset for z periodic direction
                P2P_CPU(CI,CJ);                                 //       Perform P2P_CPU kernel
              }                                                 //      Endif for periodic flag
            }                                                   //     End loop over x periodic direction
          }                                                     //    End loop over y periodic direction
        }                                                       //   End loop over z periodic direction
        listP2P[CI-CI0].pop_back();                             //   Pop last element from M2P interaction list
      }                                                         //  End while for M2P interaction list
    }                                                           // End loop over cells topdown
    listP2P.clear();                                            // Clear interaction lists
    flagP2P.clear();                                            // Clear periodic image flags
    stopTimer("evalP2P      ");                                 // Stop timer
  }                                                             // End if for kernel unit test
}

template<Equation equation>
void Evaluator<equation>::evalL2L(Cells &cells) {               // Evaluate L2L
  startTimer("evalL2L      ");                                  // Start timer
  CI0 = cells.begin();                                          // Set begin iterator
  for( C_iter CI=cells.end()-2; CI!=cells.begin()-1; --CI ) {   // Loop over cells topdown (except root cell)
    C_iter CJ = CI0 + CI->PARENT;                               //  Set source cell iterator
    L2L(CI,CJ);                                                 //  Perform L2L kernel
  }                                                             // End loop over cells topdown
  stopTimer("evalL2L      ");                                   // Stop timer
}

template<Equation equation>
void Evaluator<equation>::evalL2P(Cells &cells) {               // Evaluate L2P
  startTimer("evalL2P      ");                                  // Start timer
  for( C_iter CI=cells.begin(); CI!=cells.end(); ++CI ) {       // Loop over cells
    if( CI->NCHILD == 0 ) {                                     //  If cell is a twig
      L2P(CI);                                                  //   Perform L2P kernel
    }                                                           //  Endif for twig
  }                                                             // End loop over cells topdown
  stopTimer("evalL2P      ");                                   // Stop timer
}
