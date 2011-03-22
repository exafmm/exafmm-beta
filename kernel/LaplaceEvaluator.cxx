#include "evaluator.h"

void Evaluator::setSourceBody() {                               // Set source buffer for bodies
  startTimer("Set sourceB  ");                                  // Start timer
  for( M_iter M=sourceSize.begin(); M!=sourceSize.end(); ++M ) {// Loop over source map
    CJ = M->first;                                              //  Set source cell
    sourceBegin[CJ] = sourceHost.size() / 4;                    //  Key : iterator, Value : offset of source leafs
    for( B_iter B=CJ->LEAF; B!=CJ->LEAF+CJ->NLEAF; ++B ) {      //  Loop over leafs in source cell
      sourceHost.push_back(B->X[0]);                            //   Copy x position to GPU buffer
      sourceHost.push_back(B->X[1]);                            //   Copy y position to GPU buffer
      sourceHost.push_back(B->X[2]);                            //   Copy z position to GPU buffer
      sourceHost.push_back(B->Q);                               //   Copy mass/charge to GPU buffer
    }                                                           //  End loop over leafs
  }                                                             // End loop over source map
  stopTimer("Set sourceB  ");                                   // Stop timer
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
        targetHost.push_back(B->X[0]);                          //    Copy x position to GPU buffer
        targetHost.push_back(B->X[1]);                          //    Copy y position to GPU buffer
        targetHost.push_back(B->X[2]);                          //    Copy z position to GPU buffer
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

void Evaluator::evalP2P(Bodies &ibodies, Bodies &jbodies, bool onCPU) {// Evaluate P2P
  BI0 = ibodies.begin();                                        // Set target bodies begin iterator
  BIN = ibodies.end();                                          // Set target bodies end iterator
  BJ0 = jbodies.begin();                                        // Set source bodies begin iterator
  BJN = jbodies.end();                                          // Set source bodies end iterator
  if( onCPU ) {                                                 // If calculation is to be done on CPU
    Xperiodic = 0;                                              //  Set periodic coordinate offset
    P2P_CPU();                                                  //  Evaluate P2P kernel
  } else {                                                      // If calculation is to be done on GPU
    constHost.push_back(2*R0);                                  //  Copy domain size to GPU buffer
    for( B_iter B=BJ0; B!=BJN; ++B ) {                          //  Loop over source bodies
      sourceHost.push_back(B->X[0]);                            //  Copy x position to GPU buffer
      sourceHost.push_back(B->X[1]);                            //  Copy y position to GPU buffer
      sourceHost.push_back(B->X[2]);                            //  Copy z position to GPU buffer
      sourceHost.push_back(B->Q);                               //  Copy mass/charge to GPU buffer
    }                                                           //  End loop over source bodies
    int key = 0;                                                //  Initialize key to range of leafs in source cells
    int blocks = (BIN - BI0 - 1) / THREADS + 1;                 //  Number of thread blocks needed for this target cell
    for( int i=0; i!=blocks; ++i ) {                            //  Loop over thread blocks
      keysHost.push_back(key);                                  //   Save key to range of leafs in source cells
    }                                                           //  End loop over thread blocks
    rangeHost.push_back(1);                                     //  Save size of interaction list
    rangeHost.push_back(0);                                     //  Set begin index of leafs
    rangeHost.push_back(BJN-BJ0);                               //  Set number of leafs
    rangeHost.push_back(Icenter);                               //  Set periodic image flag
    for( B_iter B=BI0; B!=BIN; ++B ) {                          //  Loop over target bodies
      targetHost.push_back(B->X[0]);                            //   Copy x position to GPU buffer
      targetHost.push_back(B->X[1]);                            //   Copy y position to GPU buffer
      targetHost.push_back(B->X[2]);                            //   Copy z position to GPU buffer
      targetHost.push_back(0);                                  //   Initialize target value of GPU buffer
    }                                                           //  End loop over target bodies
    int numPad = blocks * THREADS - (BIN - BI0);                //  Number of elements to pad in target GPU buffer
    for( int i=0; i!=numPad; ++i ) {                            //  Loop over elements to pad
      targetHost.push_back(0);                                  //   Pad x position in GPU buffer
      targetHost.push_back(0);                                  //   Pad y position in GPU buffer
      targetHost.push_back(0);                                  //   Pad z position in GPU buffer
      targetHost.push_back(0);                                  //   Pad target value in GPU buffer
    }                                                           //  End loop over elements to pad
    P2P();                                                      //  Evaluate P2P kernel
    for( B_iter B=BI0; B!=BIN; ++B ) {                          //  Loop over target bodies
      B->pot += targetHost[4*(B-BI0)+0];                        //   Copy potential from GPU buffer
      B->acc[0] += targetHost[4*(B-BI0)+1];                     //   Copy acceleration from GPU buffer
      B->acc[1] += targetHost[4*(B-BI0)+2];                     //   Copy acceleration from GPU buffer
      B->acc[2] += targetHost[4*(B-BI0)+3];                     //   Copy acceleration from GPU buffer
    }                                                           //  End loop over target bodies
    keysHost.clear();                                           //  Clear keys vector
    rangeHost.clear();                                          //  Clear range vector
    constHost.clear();                                          //  Clear const vector
    targetHost.clear();                                         //  Clear target vector
    sourceHost.clear();                                         //  Clear source vector
  }                                                             // Endif for CPU/GPU switch
}
