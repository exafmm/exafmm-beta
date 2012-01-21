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
void Evaluator<equation>::evalP2P(Bodies &ibodies, Bodies &jbodies, bool) {// Evaluate all P2P kernels
  Xperiodic = 0;                                                // Set periodic coordinate offset
  Cells cells;                                                  // Cells to put target and source bodies
  cells.resize(2);                                              // Resize cells to put target and source bodies
  cells[0].LEAF = ibodies.begin();                              // Iterator of first target leaf
  cells[0].NDLEAF = ibodies.size();                             // Number of target leafs
  cells[1].LEAF = jbodies.begin();                              // Iterator of first source leaf
  cells[1].NDLEAF = jbodies.size();                             // Number of source leafs
  C_iter Ci = cells.begin(), Cj = cells.begin()+1;              // Iterator of target and source cells
  P2P(Ci,Cj);                                                   // Perform P2P kernel
}

template<Equation equation>
void Evaluator<equation>::evalP2M(Cells &cells) {               // Evaluate all P2M kernels
  startTimer("evalP2M      ");                                  // Start timer
  for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {          // Loop over cells
    C->M = 0;                                                   //  Initialize multipole coefficients
    C->L = 0;                                                   //  Initialize local coefficients
    if( C->NCHILD == 0 ) {                                      //  If cell is a twig
      P2M(C);                                                   //   Perform P2M kernel
    }                                                           //  Endif for twig
  }                                                             // End loop over cells
  stopTimer("evalP2M      ");                                   // Stop timer
}

template<Equation equation>
void Evaluator<equation>::evalM2M(Cells &cells, Cells &jcells) {// Evaluate all M2M kernels
  startTimer("evalM2M      ");                                  // Start timer
  Cj0 = jcells.begin();                                         // Set begin iterator
  for( C_iter Ci=cells.begin(); Ci!=cells.end(); ++Ci ) {       // Loop over target cells bottomup
    for( C_iter Cj=Cj0+Ci->CHILD; Cj!=Cj0+Ci->CHILD+Ci->NCHILD; ++Cj ) {// Loop over child cells
      M2M(Ci,Cj);                                               //   Perform M2M kernel
    }                                                           //  End loop over child cells
  }                                                             // End loop target over cells
  stopTimer("evalM2M      ");                                   // Stop timer
}

template<Equation equation>
void Evaluator<equation>::evalM2L(C_iter Ci, C_iter Cj) {       // Evaluate single M2L kernel
  M2L(Ci,Cj);                                                   // Perform M2L kernel
  NM2L++;                                                       // Count M2L kernel execution
}

template<Equation equation>
void Evaluator<equation>::evalM2L(Cells &cells) {               // Evaluate queued M2L kernels
  startTimer("evalM2L      ");                                  // Start timer
  Ci0 = cells.begin();                                          // Set begin iterator
  for( C_iter Ci=cells.begin(); Ci!=cells.end(); ++Ci ) {       // Loop over cells
    while( !listM2L[Ci-Ci0].empty() ) {                         //  While M2L interaction list is not empty
      C_iter Cj = listM2L[Ci-Ci0].back();                       //   Set source cell iterator
      Iperiodic = flagM2L[Ci-Ci0][Cj];                          //   Set periodic image flag
      int I = 0;                                                //   Initialize index of periodic image
      for( int ix=-1; ix<=1; ++ix ) {                           //   Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //    Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz, ++I ) {                  //     Loop over z periodic direction
            if( Iperiodic & (1 << I) ) {                        //      If periodic flag is on
              Xperiodic[0] = ix * 2 * R0;                       //       Coordinate offset for x periodic direction
              Xperiodic[1] = iy * 2 * R0;                       //       Coordinate offset for y periodic direction
              Xperiodic[2] = iz * 2 * R0;                       //       Coordinate offset for z periodic direction
              M2L(Ci,Cj);                                       //       Perform M2L kernel
            }                                                   //      Endif for periodic flag
          }                                                     //     End loop over x periodic direction
        }                                                       //    End loop over y periodic direction
      }                                                         //   End loop over z periodic direction
      listM2L[Ci-Ci0].pop_back();                               //   Pop last element from M2L interaction list
    }                                                           //  End while for M2L interaction list
  }                                                             // End loop over cells topdown
  listM2L.clear();                                              // Clear interaction lists
  flagM2L.clear();                                              // Clear periodic image flags
  stopTimer("evalM2L      ");                                   // Stop timer
}

template<Equation equation>
void Evaluator<equation>::evalM2P(C_iter Ci, C_iter Cj) {       // Evaluate single M2P kernel
  M2P(Ci,Cj);                                                   // Perform M2P kernel
  NM2P++;                                                       // Count M2P kernel execution
}

template<Equation equation>
void Evaluator<equation>::evalM2P(Cells &cells) {               // Evaluate queued M2P kernels
  startTimer("evalM2P      ");                                  // Start timer
  Ci0 = cells.begin();                                          // Set begin iterator
  for( C_iter Ci=cells.begin(); Ci!=cells.end(); ++Ci ) {       // Loop over cells
    while( !listM2P[Ci-Ci0].empty() ) {                         //  While M2P interaction list is not empty
      C_iter Cj = listM2P[Ci-Ci0].back();                       //   Set source cell iterator
      Iperiodic = flagM2P[Ci-Ci0][Cj];                          //   Set periodic image flag
      int I = 0;                                                //   Initialize index of periodic image
      for( int ix=-1; ix<=1; ++ix ) {                           //   Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //    Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz, ++I ) {                  //     Loop over z periodic direction
            if( Iperiodic & (1 << I) ) {                        //      If periodic flag is on
              Xperiodic[0] = ix * 2 * R0;                       //       Coordinate offset for x periodic direction
              Xperiodic[1] = iy * 2 * R0;                       //       Coordinate offset for y periodic direction
              Xperiodic[2] = iz * 2 * R0;                       //       Coordinate offset for z periodic direction
              M2P(Ci,Cj);                                       //       Perform M2P kernel
            }                                                   //      Endif for periodic flag
          }                                                     //     End loop over x periodic direction
        }                                                       //    End loop over y periodic direction
      }                                                         //   End loop over z periodic direction
      listM2P[Ci-Ci0].pop_back();                               //   Pop last element from M2P interaction list
    }                                                           //  End while for M2P interaction list
  }                                                             // End loop over cells topdown
  listM2P.clear();                                              // Clear interaction lists
  flagM2P.clear();                                              // Clear periodic image flags
  stopTimer("evalM2P      ");                                   // Stop timer
}

template<Equation equation>
void Evaluator<equation>::evalP2P(C_iter Ci, C_iter Cj) {       // Evaluate single P2P kernel
  P2P(Ci,Cj);                                                   // Perform P2P kernel
  NP2P++;                                                       // Count P2P kernel execution
}

template<Equation equation>
void Evaluator<equation>::evalP2P(Cells &cells) {               // Evaluate queued P2P kernels
  startTimer("evalP2P      ");                                  // Start timer
  Ci0 = cells.begin();                                          // Set begin iterator
  for( C_iter Ci=cells.begin(); Ci!=cells.end(); ++Ci ) {       // Loop over cells
    while( !listP2P[Ci-Ci0].empty() ) {                         //  While M2P interaction list is not empty
      C_iter Cj = listP2P[Ci-Ci0].back();                       //   Set source cell iterator
      Iperiodic = flagP2P[Ci-Ci0][Cj];                          //   Set periodic image flag
      int I = 0;                                                //   Initialize index of periodic image
      for( int ix=-1; ix<=1; ++ix ) {                           //   Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //    Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz, ++I ) {                  //     Loop over z periodic direction
            if( Iperiodic & (1 << I) ) {                        //      If periodic flag is on
              Xperiodic[0] = ix * 2 * R0;                       //       Coordinate offset for x periodic direction
              Xperiodic[1] = iy * 2 * R0;                       //       Coordinate offset for y periodic direction
              Xperiodic[2] = iz * 2 * R0;                       //       Coordinate offset for z periodic direction
              P2P(Ci,Cj);                                       //       Perform P2P kernel
            }                                                   //      Endif for periodic flag
          }                                                     //     End loop over x periodic direction
        }                                                       //    End loop over y periodic direction
      }                                                         //   End loop over z periodic direction
      listP2P[Ci-Ci0].pop_back();                               //   Pop last element from M2P interaction list
    }                                                           //  End while for M2P interaction list
  }                                                             // End loop over cells topdown
  listP2P.clear();                                              // Clear interaction lists
  flagP2P.clear();                                              // Clear periodic image flags
  stopTimer("evalP2P      ");                                   // Stop timer
}

template<Equation equation>
void Evaluator<equation>::evalL2L(Cells &cells) {               // Evaluate all L2L kernels
  startTimer("evalL2L      ");                                  // Start timer
  Ci0 = cells.begin();                                          // Set begin iterator
  for( C_iter Ci=cells.end()-2; Ci!=cells.begin()-1; --Ci ) {   // Loop over cells topdown (except root cell)
    C_iter Cj = Ci0 + Ci->PARENT;                               //  Set source cell iterator
    L2L(Ci,Cj);                                                 //  Perform L2L kernel
  }                                                             // End loop over cells topdown
  stopTimer("evalL2L      ");                                   // Stop timer
}

template<Equation equation>
void Evaluator<equation>::evalL2P(Cells &cells) {               // Evaluate all L2P kernels
  startTimer("evalL2P      ");                                  // Start timer
  for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {          // Loop over cells
    if( C->NCHILD == 0 ) {                                      //  If cell is a twig
      L2P(C);                                                   //   Perform L2P kernel
    }                                                           //  Endif for twig
  }                                                             // End loop over cells topdown
  stopTimer("evalL2P      ");                                   // Stop timer
}

template<Equation equation>
void Evaluator<equation>::evalEwaldReal(C_iter Ci, C_iter Cj) { // Evaluate single Ewald real kernel
  EwaldReal(Ci,Cj);                                             // Perform Ewald real kernel
}

template<Equation equation>
void Evaluator<equation>::evalEwaldReal(Cells &cells) {         // Evaluate queued Ewald real kernels
  startTimer("evalEwaldReal");                                  // Start timer
  Ci0 = cells.begin();                                          // Set begin iterator
  for( C_iter Ci=cells.begin(); Ci!=cells.end(); ++Ci ) {       // Loop over cells
    while( !listP2P[Ci-Ci0].empty() ) {                         //  While M2P interaction list is not empty
      C_iter Cj = listP2P[Ci-Ci0].back();                       //   Set source cell iterator
      EwaldReal(Ci,Cj);                                         //   Perform Ewald real kernel
      listP2P[Ci-Ci0].pop_back();                               //   Pop last element from M2P interaction list
    }                                                           //  End while for M2P interaction list
  }                                                             // End loop over cells topdown
  listP2P.clear();                                              // Clear interaction lists
  flagP2P.clear();                                              // Clear periodic image flags
  stopTimer("evalEwaldReal");                                   // Stop timer
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
  C_iter Ci = cells.begin(), Cj = cells.begin()+1;              // Artificial target & source cell
  Ci->X = 0;                                                    // Set coordinates of target cell
  Ci->NDLEAF = 10;                                              // Number of leafs in target cell
  Ci->LEAF = ibodies.begin();                                   // Leaf iterator in target cell
  Cj->X = 1;                                                    // Set coordinates of source cell
  Cj->NDLEAF = 1000;                                            // Number of leafs in source cell
  Cj->LEAF = jbodies.begin();                                   // Leaf iterator in source cell
  startTimer("P2P kernel   ");                                  // Start timer
  for( int i=0; i!=1; ++i ) P2P(Ci,Cj);                         // Perform P2P kernel
  timeP2P = stopTimer("P2P kernel   ") / 10000;                 // Stop timer
  startTimer("M2L kernel   ");                                  // Start timer
  for( int i=0; i!=1000; ++i ) M2L(Ci,Cj);                      // Perform M2L kernel
  timeM2L = stopTimer("M2L kernel   ") / 1000;                  // Stop timer
  startTimer("M2P kernel   ");                                  // Start timer
  for( int i=0; i!=100; ++i ) M2P(Ci,Cj);                       // Perform M2P kernel
  timeM2P = stopTimer("M2P kernel   ") / 1000;                  // Stop timer
}

#if QUARK
template<Equation equation>
inline void interactQuark(Quark *quark) {
  Evaluator<equation> *E;
  C_iter CI, CJ, Ci0, Cj0;
  quark_unpack_args_5(quark,E,CI,CJ,Ci0,Cj0);
  ThreadTrace beginTrace;
  E->startTracer(beginTrace);
  PairQueue privateQueue;
  Pair pair(CI,CJ);
  privateQueue.push(pair);
  while( !privateQueue.empty() ) {
    Pair Cij = privateQueue.front();
    privateQueue.pop();
    if(splitFirst(Cij.first,Cij.second)) {
      C_iter C = Cij.first;
      for( C_iter Ci=Ci0+C->CHILD; Ci!=Ci0+C->CHILD+C->NCHILD; ++Ci ) {
        E->interact(Ci,Cij.second,privateQueue);
      }
    } else {
      C_iter C = Cij.second;
      for( C_iter Cj=Cj0+C->CHILD; Cj!=Cj0+C->CHILD+C->NCHILD; ++Cj ) {
        E->interact(Cij.first,Cj,privateQueue);
      }
    }
  }
  E->stopTracer(beginTrace,0x0000ff);
}

template<Equation equation>
void Evaluator<equation>::interact(C_iter Ci, C_iter Cj, Quark *quark) {
  char string[256];
  sprintf(string,"%d %d",int(Ci-Ci0),int(Cj-Cj0));
  Quark_Task_Flags tflags = Quark_Task_Flags_Initializer;
  QUARK_Task_Flag_Set(&tflags,TASK_LABEL,intptr_t(string) );
  QUARK_Insert_Task(quark,interactQuark<equation>,&tflags,
                    sizeof(Evaluator),this,NODEP,
                    sizeof(Cell),&*Ci,OUTPUT,
                    sizeof(Cell),&*Cj,NODEP,
                    sizeof(Cell),&*Ci0,NODEP,
                    sizeof(Cell),&*Cj0,NODEP,
                    0);
}
#endif
