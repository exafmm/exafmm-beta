#include "evaluator.h"

void Evaluator::setSourceBody() {                               // Set source buffer for bodies
  startTimer("Set sourceB  ");                                  // Start timer
  for( M_iter M=sourceSize.begin(); M!=sourceSize.end(); ++M ) {// Loop over source map
    CJ = M->first;                                              //  Set source cell
    sourceBegin[CJ] = sourceHost.size() / 4;                    //  Key : iterator, Value : offset of source leafs
    for( B_iter B=CJ->LEAF; B!=CJ->LEAF+CJ->NLEAF; ++B ) {      //  Loop over leafs in source cell
      sourceHost.push_back(B->pos[0]);                          //   Copy x position to GPU buffer
      sourceHost.push_back(B->pos[1]);                          //   Copy y position to GPU buffer
      sourceHost.push_back(B->pos[2]);                          //   Copy z position to GPU buffer
      sourceHost.push_back(B->scal);                            //   Copy mass/charge to GPU buffer
    }                                                           //  End loop over leafs
  }                                                             // End loop over source map
  stopTimer("Set sourceB  ");                                   // Stop timer
}

void Evaluator::setSourceCell(bool isM=true) {                  // Set source buffer for cells
  startTimer("Set sourceC  ");                                  // Start timer
  for( M_iter M=sourceSize.begin(); M!=sourceSize.end(); ++M ) {// Loop over source map
    CJ = M->first;                                              //  Set source cell
    sourceBegin[CJ] = sourceHost.size();                        //  Key : iterator, Value : offset of sources
    sourceHost.push_back(CJ->X[0]);                             //  Copy x position to GPU buffer
    sourceHost.push_back(CJ->X[1]);                             //  Copy y position to GPU buffer
    sourceHost.push_back(CJ->X[2]);                             //  Copy z position to GPU buffer
    if( isM ) {                                                 //  If source is M
      for( int i=0; i!=NCOEF; ++i ) {                           //   Loop over coefs in source cell
        sourceHost.push_back((CJ->M[i]).real());                //    Copy real multipole to GPU buffer
        sourceHost.push_back((CJ->M[i]).imag());                //    Copy imaginary multipole to GPU buffer
      }                                                         //    End loop over coefs
    } else {                                                    //  If source is L
      for( int i=0; i!=NCOEF; ++i ) {                           //   Loop over coefs in source cell
        sourceHost.push_back((CJ->L[i]).real());                //    Copy real multipole to GPU buffer
        sourceHost.push_back((CJ->L[i]).imag());                //    Copy imaginary multipole to GPU buffer
      }                                                         //    End loop over coefs
    }                                                           //  Endif for source type
  }                                                             // End loop over source map
  stopTimer("Set sourceC  ");                                   // Stop timer
}

void Evaluator::setTargetBody(Cells &cells, Lists lists, Maps flags) {// Set target buffer for bodies
  startTimer("Set targetB  ");                                  // Start timer
  int key = 0;                                                  // Initialize key to range of coefs in source cells
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {              // Loop over target cells
    if( !lists[CI-CI0].empty() ) {                              //  If the interation list is not empty
      BI0 = CI->LEAF;                                           //   Set target bodies begin iterator
      BIN = CI->LEAF + CI->NLEAF;                               //   Set target bodies end iterator
      int blocks = (BIN - BI0 - 1) / THREADS + 1;               //   Number of thread blocks needed for this target cell
      for( int i=0; i!=blocks; ++i ) {                          //   Loop over thread blocks
        keysHost.push_back(key);                                //    Save key to range of leafs in source cells
      }                                                         //   End loop over thread blocks
      key += 3*lists[CI-CI0].size()+1;                          //   Increment key counter
      rangeHost.push_back(lists[CI-CI0].size());                //   Save size of interaction list
      for( L_iter L=lists[CI-CI0].begin(); L!=lists[CI-CI0].end(); ++L ) {//  Loop over interaction list
        CJ = *L;                                                //   Set source cell
        rangeHost.push_back(sourceBegin[CJ]);                   //    Set begin index of coefs in source cell
        rangeHost.push_back(sourceSize[CJ]);                    //    Set number of coefs in source cell
        rangeHost.push_back(flags[CI-CI0][CJ]);                 //    Set periodic image flag of source cell
      }                                                         //   End loop over interaction list
      targetBegin[CI] = targetHost.size() / 4;                  //   Key : iterator, Value : offset of target leafs
      for( B_iter B=BI0; B!=BIN; ++B ) {                        //   Loop over leafs in target cell
        targetHost.push_back(B->pos[0]);                        //    Copy x position to GPU buffer
        targetHost.push_back(B->pos[1]);                        //    Copy y position to GPU buffer
        targetHost.push_back(B->pos[2]);                        //    Copy z position to GPU buffer
        targetHost.push_back(0);                                //    Initialize target value of GPU buffer
      }                                                         //   End loop over leafs
      int numPad = blocks * THREADS - (BIN - BI0);              //   Number of elements to pad in target GPU buffer
      for( int i=0; i!=numPad; ++i ) {                          //   Loop over elements to pad
        targetHost.push_back(0);                                //    Pad x position in GPU buffer
        targetHost.push_back(0);                                //    Pad y position in GPU buffer
        targetHost.push_back(0);                                //    Pad z position in GPU buffer
        targetHost.push_back(0);                                //    Pad target value in GPU buffer
      }                                                         //   End loop over elements to pad
    }                                                           //  End if for empty interation list
  }                                                             // End loop over target cells
  stopTimer("Set targetB  ");                                   // Stop timer
}

void Evaluator::setTargetCell(Cells &cells, Lists lists, Maps flags) {// Set target buffer for cells
  startTimer("Set targetC  ");                                  // Start timer
  int key = 0;                                                  // Initialize key to range of coefs in target cells
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {              // Loop over target cells
    if( !lists[CI-CI0].empty() ) {                              //  If the interation list is not empty
      keysHost.push_back(key);                                  //   Save key to range of coefs in target cells
      key += 3*lists[CI-CI0].size()+1;                          //   Increment key counter
      rangeHost.push_back(lists[CI-CI0].size());                //   Save size of interaction list
      for( L_iter L=lists[CI-CI0].begin(); L!=lists[CI-CI0].end(); ++L ) {//  Loop over interaction list
        CJ = *L;                                                //    Set target cell
        int begin = sourceBegin[CJ];                            //    Get begin index of coefs in source cell
        int size = sourceSize[CJ];                              //    Get number of coefs in source cell
        rangeHost.push_back(begin);                             //    Set begin index of coefs in source cell
        rangeHost.push_back(size);                              //    Set number of coefs in source cell
        rangeHost.push_back(flags[CI-CI0][CJ]);                 //    Set periodic image flag of source cell
      }                                                         //   End loop over interaction list
      targetBegin[CI] = targetHost.size();                      //   Key : iterator, Value : offset of target coefs
      targetHost.push_back(CI->X[0]);                           //   Copy x position to GPU buffer
      targetHost.push_back(CI->X[1]);                           //   Copy y position to GPU buffer
      targetHost.push_back(CI->X[2]);                           //   Copy z position to GPU buffer
      for( int i=0; i!=2*NCOEF; ++i ) {                         //   Loop over coefs in target cell
        targetHost.push_back(0);                                //    Pad GPU buffer
      }                                                         //   End loop over coefs
      int numPad = 2 * THREADS - 2 * NCOEF - 3;                 //   Number of elements to pad in target GPU buffer
      assert(numPad >= 0);                                      //   THREADS must be large enough
      for( int i=0; i!=numPad; ++i ) {                          //   Loop over elements to pad
        targetHost.push_back(0);                                //    Pad GPU buffer
      }                                                         //   End loop over elements to pad
    }                                                           //  End if for empty interation list
  }                                                             // End loop over target cells
  stopTimer("Set targetC  ");                                   // Stop timer
}

void Evaluator::getTargetBody(Cells &cells, Lists &lists) {     // Get body values from target buffer
  startTimer("Get targetB  ");                                  // Start timer
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {              // Loop over target cells
    if( !lists[CI-CI0].empty() ) {                              //  If the interation list is not empty
      BI0 = CI->LEAF;                                           //   Set target bodies begin iterator
      BIN = CI->LEAF + CI->NLEAF;                               //   Set target bodies end iterator
      int begin = targetBegin[CI];                              //   Offset of target leafs
      for( B_iter B=BI0; B!=BIN; ++B ) {                        //   Loop over target bodies
        B->pot    += targetHost[4*(begin+B-BI0)+0];             //    Copy potential from GPU buffer
        B->acc[0] += targetHost[4*(begin+B-BI0)+1];             //    Copy acceleration from GPU buffer
        B->acc[1] += targetHost[4*(begin+B-BI0)+2];             //    Copy acceleration from GPU buffer
        B->acc[2] += targetHost[4*(begin+B-BI0)+3];             //    Copy acceleration from GPU buffer
      }                                                         //   End loop over target bodies
      lists[CI-CI0].clear();                                    //   Clear interaction list
    }                                                           //  End if for empty interation list
  }                                                             // End loop over target cells
  stopTimer("Get targetB  ");                                   // Stop timer
}

void Evaluator::getTargetCell(Cells &cells, Lists &lists, bool isM=true) {// Get body values from target buffer
  startTimer("Get targetC  ");                                  // Start timer
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {              // Loop over target cells
    if( !lists[CI-CI0].empty() ) {                              //  If the interation list is not empty
      int begin = targetBegin[CI];                              //   Offset of target coefs
      if( isM ) {                                               //   If target is M
        for( int i=0; i!=NCOEF; ++i ) {                         //    Loop over coefs in target cell
          CI->M[i].real() += targetHost[begin+2*i+0];           //     Copy real target values from GPU buffer
          CI->M[i].imag() += targetHost[begin+2*i+1];           //     Copy imaginary target values from GPU buffer
        }                                                       //    End loop over coefs
      } else {                                                  //   If target is L
        for( int i=0; i!=NCOEF; ++i ) {                         //    Loop over coefs in target cell
          CI->L[i].real() += targetHost[begin+2*i+0];           //     Copy real target values from GPU buffer
          CI->L[i].imag() += targetHost[begin+2*i+1];           //     Copy imaginary target values from GPU buffer
        }                                                       //    End loop over coefs
      }                                                         //   Endif for target type
      lists[CI-CI0].clear();                                    //   Clear interaction list
    }                                                           //  End if for empty interation list
  }                                                             // End loop over target cells
  stopTimer("Get targetC  ");                                   // Stop timer
}

void Evaluator::clearBuffers() {                                // Clear GPU buffers
  startTimer("Clear buffer ");                                  // Start timer
  keysHost.clear();                                             // Clear keys vector
  rangeHost.clear();                                            // Clear range vector
  constHost.clear();                                            // Clear const vector
  sourceHost.clear();                                           // Clear source vector
  targetHost.clear();                                           // Clear target vector
  sourceBegin.clear();                                          // Clear map for offset of source cells
  sourceSize.clear();                                           // Clear map for size of source cells
  targetBegin.clear();                                          // Clear map for offset of target cells
  stopTimer("Clear buffer ");                                   // Stop timer
}

void Evaluator::evalP2P(Bodies &ibodies, Bodies &jbodies, bool isPeriodic) {
  BI0 = ibodies.begin();                                        // Set target bodies begin iterator
  BIN = ibodies.end();                                          // Set target bodies end iterator
  BJ0 = jbodies.begin();                                        // Set source bodies begin iterator
  BJN = jbodies.end();                                          // Set source bodies end iterator
  constHost.push_back(2*R0);                                    // Copy domain size to GPU buffer
  for( B_iter B=BJ0; B!=BJN; ++B ) {                            // Loop over source bodies
    sourceHost.push_back(B->pos[0]);                            // Copy x position to GPU buffer
    sourceHost.push_back(B->pos[1]);                            // Copy y position to GPU buffer
    sourceHost.push_back(B->pos[2]);                            // Copy z position to GPU buffer
    sourceHost.push_back(B->scal);                              // Copy mass/charge to GPU buffer
  }                                                             // End loop over source bodies
  int key = 0;                                                  // Initialize key to range of leafs in source cells
  int blocks = (BIN - BI0 - 1) / THREADS + 1;                   // Number of thread blocks needed for this target cell
  for( int i=0; i!=blocks; ++i ) {                              // Loop over thread blocks
    keysHost.push_back(key);                                    //  Save key to range of leafs in source cells
  }                                                             // End loop over thread blocks
  rangeHost.push_back(1);                                       // Save size of interaction list
  rangeHost.push_back(0);                                       // Set begin index of leafs
  rangeHost.push_back(BJN-BJ0);                                 // Set number of leafs
  if( isPeriodic ) {
    rangeHost.push_back((1 << 27) - 1);                         // Set periodic image flag
  } else {
    rangeHost.push_back(1 << Icenter);                          // Set periodic image flag
  }
  for( B_iter B=BI0; B!=BIN; ++B ) {                            // Loop over target bodies
    targetHost.push_back(B->pos[0]);                            //  Copy x position to GPU buffer
    targetHost.push_back(B->pos[1]);                            //  Copy y position to GPU buffer
    targetHost.push_back(B->pos[2]);                            //  Copy z position to GPU buffer
    targetHost.push_back(0);                                    //  Initialize target value of GPU buffer
  }                                                             // End loop over target bodies
  int numPad = blocks * THREADS - (BIN - BI0);                  // Number of elements to pad in target GPU buffer
  for( int i=0; i!=numPad; ++i ) {                              // Loop over elements to pad
    targetHost.push_back(0);                                    //  Pad x position in GPU buffer
    targetHost.push_back(0);                                    //  Pad y position in GPU buffer
    targetHost.push_back(0);                                    //  Pad z position in GPU buffer
    targetHost.push_back(0);                                    //  Pad target value in GPU buffer
  }                                                             // End loop over elements to pad
  P2P();                                                        // Evaluate P2P kernel
  for( B_iter B=BI0; B!=BIN; ++B ) {                            // Loop over target bodies
    B->pot += targetHost[4*(B-BI0)+0];                          //  Copy potential from GPU buffer
    B->acc[0] += targetHost[4*(B-BI0)+1];                       //  Copy acceleration from GPU buffer
    B->acc[1] += targetHost[4*(B-BI0)+2];                       //  Copy acceleration from GPU buffer
    B->acc[2] += targetHost[4*(B-BI0)+3];                       //  Copy acceleration from GPU buffer
  }                                                             // End loop over target bodies
  keysHost.clear();                                             // Clear keys vector
  rangeHost.clear();                                            // Clear range vector
  constHost.clear();                                            // Clear const vector
  targetHost.clear();                                           // Clear target vector
  sourceHost.clear();                                           // Clear source vector
}

void Evaluator::evalP2M(Cells &cells) {                         // Evaluate P2M
  startTimer("Get list     ");                                  // Start timer
  CI0 = cells.begin();                                          // Set begin iterator for target
  CJ0 = cells.begin();                                          // Set begin iterator for source
  constHost.push_back(2*R0);                                    // Copy domain size to GPU buffer
  Lists listP2M(cells.size());                                  // Define P2M interation list vector
  Maps  flagP2M(cells.size());                                  // Define P2M periodic image flag
  for( CJ=cells.begin(); CJ!=cells.end(); ++CJ ) {              // Loop over target cells
    CJ->M = CJ->L = 0;                                          //  Initialize multipole & local coefficients
    listP2M[CJ-CJ0].push_back(CJ);                              //  Push source cell into P2M interaction list
    flagP2M[CJ-CJ0][CJ] |= 1 << Icenter;                        //  Flip bit of periodic image flag
    sourceSize[CJ] = CJ->NLEAF;                                 //  Key : iterator, Value : number of leafs
  }                                                             // End loop over source map
  stopTimer("Get list     ");                                   // Stop timer
  setSourceBody();                                              // Set source buffer for bodies
  setTargetCell(cells,listP2M,flagP2M);                         // Set target buffer for cells
  P2M();                                                        // Evaluate P2M kernel
  getTargetCell(cells,listP2M);                                 // Get body values from target buffer
  clearBuffers();                                               // Clear GPU buffers
}

void Evaluator::evalM2M(Cells &cells) {                         // Evaluate M2M
  CI0 = cells.begin();                                          // Set begin iterator for target
  CJ0 = cells.begin();                                          // Set begin iterator for source
  constHost.push_back(2*R0);                                    // Copy domain size to GPU buffer
  Lists listM2M(cells.size());                                  // Define M2M interation list vector
  Maps  flagM2M(cells.size());                                  // Define M2M periodic image flag
  int level = getLevel(CJ0->I);                                 // Level of twig
  while( level != 0 ) {                                         // While level of source is not root level
    startTimer("Get list     ");                                //  Start timer
    for( CJ=cells.begin(); CJ!=cells.end(); ++CJ ) {            //  Loop over cells bottomup (except root cell)
      if( getLevel(CJ->I) == level ) {                          //   If source cell is at current level
        CI = CJ0 + CJ->PARENT;                                  //    Set target cell iterator
        listM2M[CI-CI0].push_back(CJ);                          //    Push source cell into M2M interaction list
        flagM2M[CI-CI0][CJ] |= 1 << Icenter;                    //    Flip bit of periodic image flag
        sourceSize[CJ] = 2 * NCOEF;                             //    Key : iterator, Value : number of coefs
      }                                                         //   Endif for current level
    }                                                           //  End loop over cells
    stopTimer("Get list     ");                                 //  Stop timer
    setSourceCell();                                            //  Set source buffer for cells
    setTargetCell(cells,listM2M,flagM2M);                       //  Set target buffer for cells
    M2M();                                                      //  Evaluate M2M kernel
    getTargetCell(cells,listM2M);                               //  Get body values from target buffer
    clearBuffers();                                             //  Clear GPU buffers
    level--;                                                    //  Decrement level
  }                                                             // End while loop over levels
}

void Evaluator::evalM2L(Cells &cells) {                         // Evaluate M2L
  startTimer("Get list     ");                                  // Start timer
  CI0 = cells.begin();                                          // Set begin iterator
  constHost.push_back(2*R0);                                    // Copy domain size to GPU buffer
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {              // Loop over target cells
    for( L_iter L=listM2L[CI-CI0].begin(); L!=listM2L[CI-CI0].end(); ++L ) {//  Loop over interaction list
      CJ = *L;                                                  //   Set source cell
      sourceSize[CJ] = 2 * NCOEF;                               //   Key : iterator, Value : number of coefs
    }                                                           //  End loop over interaction list
  }                                                             // End loop over target cells
  stopTimer("Get list     ");                                   // Stop timer
  setSourceCell();                                              // Set source buffer for cells
  setTargetCell(cells,listM2L,flagM2L);                         // Set target buffer for cells
  M2L();                                                        // Evaluate M2L kernel
  getTargetCell(cells,listM2L,false);                           // Get body values from target buffer
  clearBuffers();                                               // Clear GPU buffers
}

void Evaluator::evalM2P(Cells &cells) {                         // Evaluate M2P
  startTimer("Get list     ");                                  // Start timer
  CI0 = cells.begin();                                          // Set begin iterator for target
  CJ0 = cells.begin();                                          // Set begin iterator for source
  constHost.push_back(2*R0);                                    // Copy domain size to GPU buffer
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {              // Loop over target cells
    for( L_iter L=listM2P[CI-CI0].begin(); L!=listM2P[CI-CI0].end(); ++L ) {//  Loop over interaction list
      CJ = *L;                                                  //   Set source cell
      sourceSize[CJ] = 2 * NCOEF;                               //   Key : iterator, Value : number of coefs
    }                                                           //  End loop over interaction list
  }                                                             // End loop over target cells
  stopTimer("Get list     ");                                   // Stop timer
  setSourceCell();                                              // Set source buffer for cells
  setTargetBody(cells,listM2P,flagM2P);                         // Set target buffer for bodies
  M2P();                                                        // Evaluate M2P kernel
  getTargetBody(cells,listM2P);                                 // Get body values from target buffer
  clearBuffers();                                               // Clear GPU buffers
}

void Evaluator::evalP2P(Cells &cells) {                         // Evaluate P2P
  startTimer("Get list     ");                                  // Start timer
  CI0 = cells.begin();                                          // Set begin iterator
  constHost.push_back(2*R0);                                    // Copy domain size to GPU buffer
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {              // Loop over target cells
    for( L_iter L=listP2P[CI-CI0].begin(); L!=listP2P[CI-CI0].end(); ++L ) {//  Loop over interaction list
      CJ = *L;                                                  //   Set source cell
      sourceSize[CJ] = CJ->NLEAF;                               //   Key : iterator, Value : number of leafs
    }                                                           //  End loop over interaction list
  }                                                             // End loop over target cells
  stopTimer("Get list     ");                                   // Stop timer
  setSourceBody();                                              // Set source buffer for bodies
  setTargetBody(cells,listP2P,flagP2P);                         // Set target buffer for bodies
  P2P();                                                        // Evaluate P2P kernel
  getTargetBody(cells,listP2P);                                 // Get body values from target buffer
  clearBuffers();                                               // Clear GPU buffers
}

void Evaluator::evalL2L(Cells &cells) {                         // Evaluate L2L
  CI0 = cells.begin();                                          // Set begin iterator
  constHost.push_back(2*R0);                                    // Copy domain size to GPU buffer
  Lists listL2L(cells.size());                                  // Define L2L interation list vector
  Maps  flagL2L(cells.size());                                  // Define L2L periodic image flag
  int maxLevel = getLevel(CI0->I);                              // Level of twig
  int level = 1;                                                // Start level from 1
  while( level != maxLevel+1 ) {                                // While level of source is not root level
    startTimer("Get list     ");                                //  Start timer
    for( CI=cells.end()-2; CI!=cells.begin()-1; --CI ) {        //  Loop over cells topdown (except root cell)
      if( getLevel(CI->I) == level ) {                          //   If target cell is at current level
        CJ = CI0 + CI->PARENT;                                  //    Set source cell iterator
        listL2L[CI-CI0].push_back(CJ);                          //    Push source cell into L2L interaction list
        flagL2L[CI-CI0][CJ] |= 1 << Icenter;                    //    Flip bit of periodic image flag
        if( sourceSize[CJ] == 0 ) {                             //    If the source cell has not been stored yet
          sourceSize[CJ] = 2 * NCOEF;                           //     Key : iterator, Value : number of coefs
        }                                                       //    Endif for current level
      }                                                         //   Endif for stored source cell
    }                                                           //  End loop over cells topdown
    stopTimer("Get list     ");                                 //  Stop timer
    setSourceCell(false);                                       //  Set source buffer for cells
    setTargetCell(cells,listL2L,flagL2L);                       //  Set target buffer for cells
    L2L();                                                      //  Evaluate L2L kernel
    getTargetCell(cells,listL2L,false);                         //  Get body values from target buffer
    clearBuffers();                                             //  Clear GPU buffers
    level++;                                                    //  Increment level
  }                                                             // End while loop over levels
}

void Evaluator::evalL2P(Cells &cells) {                         // Evaluate L2P
  startTimer("Get list     ");                                  // Start timer
  CI0 = cells.begin();                                          // Set begin iterator
  constHost.push_back(2*R0);                                    // Copy domain size to GPU buffer
  Lists listL2P(cells.size());                                  // Define L2P interation list vector
  Maps  flagL2P(cells.size());                                  // Define L2P periodic image flag
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {              // Loop over cells
    if( CI->NCHILD == 0 ) {                                     //  If cell is a twig evaluate L2P kernel
      listL2P[CI-CI0].push_back(CI);                            //   Push source cell into L2P interaction list
      flagL2P[CI-CI0][CI] |= 1 << Icenter;                      //   Flip bit of periodic image flag
      sourceSize[CI] = 2 * NCOEF;                               //   Key : iterator, Value : number of coefs
    }                                                           //  Endif for twig cells
  }                                                             // End loop over cells topdown
  stopTimer("Get list     ");                                   // Stop timer
  setSourceCell(false);                                         // Set source buffer for cells
  setTargetBody(cells,listL2P,flagL2P);                         // Set target buffer for bodies
  L2P();                                                        // Evaluate L2P kernel
  getTargetBody(cells,listL2P);                                 // Get body values from target buffer
  clearBuffers();                                               // Clear GPU buffers
}
