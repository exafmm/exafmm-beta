#include <sys/time.h>
#include "evaluator.h"

double get_gpu_time(void) {
  cudaThreadSynchronize();
  struct timeval tv;                                            // Time value
  gettimeofday(&tv, NULL);                                      // Get time of day in seconds and microseconds
  return double(tv.tv_sec+tv.tv_usec*1e-6);                     // Combine seconds and microseconds and return
}

void Evaluator::evalP2P(Bodies &ibodies, Bodies &jbodies) {
  B_iter BI0 = ibodies.begin();                                 // Set target bodies begin iterator
  B_iter BIN = ibodies.end();                                   // Set target bodies end iterator
  B_iter BJ0 = jbodies.begin();                                 // Set source bodies begin iterator
  B_iter BJN = jbodies.end();                                   // Set source bodies end iterator
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
    B->pot += targetHost[4*(B-BI0)+3];                          //  Copy target values from GPU buffer
  }                                                             // End loop over target bodies
  keysHost.clear();                                             // Clear keys vector
  rangeHost.clear();                                            // Clear range vector
  targetHost.clear();                                           // Clear target vector
  sourceHost.clear();                                           // Clear source vector
}

void Evaluator::evalP2M(Cells &cells) {                         // Evaluate P2M
  for( CJ=cells.begin(); CJ!=cells.end(); ++CJ ) {              // Loop over cells
    CJ->M = CJ->L = 0;                                          //  Initialize multipole & local coefs
    P2M();                                                      //  Evaluate P2M kernel
  }                                                             // End loop over cells
}

void Evaluator::evalM2M(Cells &cells) {                         // Evaluate M2M
  CJ0 = cells.begin();                                          // Set begin iterator
  for( CJ=cells.begin(); CJ!=cells.end()-1; ++CJ ) {            // Loop over cells bottomup (except root cell)
    CI = CJ0 + CJ->PARENT;                                      //  Set target cell iterator
    M2M();                                                      //  Evaluate M2M kernel
  }                                                             // End loop over cells
}

void Evaluator::evalM2L(Cells &cells) {                         // Evaluate M2L
  CI0 = cells.begin();                                          // Set begin iterator
  Map sourceSize;                                               // Define map for size of source cells in listM2L
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {              // Loop over target cells
    for( L_iter L=listM2L[CI-CI0].begin(); L!=listM2L[CI-CI0].end(); ++L ) {//  Loop over interaction list
      CJ = *L;                                                  //   Set source cell
      sourceSize[CJ] = 2 * NCOEF;                               //   Key : iterator, Value : number of coefs
    }                                                           //  End loop over interaction list
  }                                                             // End loop over target cells
  Map sourceBegin;                                              // Define map for offset of source cells in listM2L
  for( M_iter M=sourceSize.begin(); M!=sourceSize.end(); ++M ) {// Loop over source map
    CJ = M->first;                                              //  Set source cell
    sourceBegin[CJ] = sourceHost.size();                        //  Key : iterator, Value : offset of sources
    sourceHost.push_back(CJ->X[0]);                             //  Copy x position to GPU buffer
    sourceHost.push_back(CJ->X[1]);                             //  Copy y position to GPU buffer
    sourceHost.push_back(CJ->X[2]);                             //  Copy z position to GPU buffer
    for( int i=0; i!=NCOEF; ++i ) {                             //  Loop over coefs in source cell
      sourceHost.push_back((CJ->M[i]).real());                  //   Copy real multipole to GPU buffer
      sourceHost.push_back((CJ->M[i]).imag());                  //   Copy imaginary multipole to GPU buffer
    }                                                           //  End loop over coefs
  }                                                             // End loop over source map
  Map targetBegin;                                              // Define map for offset of target cells
  int key = 0;                                                  // Initialize key to range of coefs in source cells
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {              // Loop over target cells
    if( !listM2L[CI-CI0].empty() ) {                            //  If the interation list is not empty
      keysHost.push_back(key);                                  //   Save key to range of coefs in source cells
      key += 2*listM2L[CI-CI0].size()+1;                        //   Increment key counter
      rangeHost.push_back(listM2L[CI-CI0].size());              //   Save size of interaction list
      for( L_iter L=listM2L[CI-CI0].begin(); L!=listM2L[CI-CI0].end(); ++L ) {//  Loop over interaction list
        CJ = *L;                                                //   Set source cell
        rangeHost.push_back(sourceBegin[CJ]);                   //    Set begin index of coefs in source cell
        rangeHost.push_back(sourceSize[CJ]);                    //    Set number of coefs in source cell
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
  M2L();                                                        // Evaluate M2L kernel
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {              // Loop over target cells
    if( !listM2L[CI-CI0].empty() ) {                            //  If the interation list is not empty
      int begin = targetBegin[CI];                              //   Offset of target coefs
      for( int i=0; i!=NCOEF; ++i ) {                           //   Loop over target coefs
        CI->L[i].real() += targetHost[begin+2*i+0];             //    Copy real target values from GPU buffer
        CI->L[i].imag() += targetHost[begin+2*i+1];             //    Copy imaginary target values from GPU buffer
      }                                                         //   End loop over target coefs
      listM2L[CI-CI0].clear();                                  //   Clear M2L interaction list
    }                                                           //  End if for empty interation list
  }                                                             // End loop over target cells
  keysHost.clear();                                             // Clear keys vector
  rangeHost.clear();                                            // Clear range vector
  targetHost.clear();                                           // Clear target vector
  sourceHost.clear();                                           // Clear source vector
}

void Evaluator::evalM2P(Cells &cells) {                         // Evaluate M2P
  CI0 = cells.begin();                                          // Set begin iterator
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {              // Loop over cells
    while( !listM2P[CI-CI0].empty() ) {                         //  While M2P interaction list is not empty
      CJ = listM2P[CI-CI0].back();                              //   Set source cell iterator
      M2P();                                                    //   Evaluate M2P kernel
      listM2P[CI-CI0].pop_back();                               //   Pop last element from M2P interaction list
    }                                                           //  End while for M2P interaction list
  }                                                             // End loop over cells topdown
}

void Evaluator::evalP2P(Cells &cells) {                         // Evaluate P2P
  CI0 = cells.begin();                                          // Set begin iterator
  Map sourceSize;                                               // Define map for size of source cells in listP2P
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {              // Loop over target cells
    for( L_iter L=listP2P[CI-CI0].begin(); L!=listP2P[CI-CI0].end(); ++L ) {//  Loop over interaction list
      CJ = *L;                                                  //   Set source cell
      sourceSize[CJ] = CJ->NLEAF;                               //   Key : iterator, Value : number of leafs
    }                                                           //  End loop over interaction list
  }                                                             // End loop over target cells
  Map sourceBegin;                                              // Define map for offset of source cells in listP2P
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
  Map targetBegin;                                              // Define map for offset of target cells
  int key = 0;                                                  // Initialize key to range of leafs in source cells
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {              // Loop over target cells
    if( !listP2P[CI-CI0].empty() ) {                            //  If the interation list is not empty
      BI0 = CI->LEAF;                                           //   Set target bodies begin iterator
      BIN = CI->LEAF + CI->NLEAF;                               //   Set target bodies end iterator
      int blocks = (BIN - BI0 - 1) / THREADS + 1;               //   Number of thread blocks needed for this target cell
      for( int i=0; i!=blocks; ++i ) {                          //   Loop over thread blocks
        keysHost.push_back(key);                                //    Save key to range of leafs in source cells
      }                                                         //   End loop over thread blocks
      key += 2*listP2P[CI-CI0].size()+1;                        //   Increment key counter
      rangeHost.push_back(listP2P[CI-CI0].size());              //   Save size of interaction list
      for( L_iter L=listP2P[CI-CI0].begin(); L!=listP2P[CI-CI0].end(); ++L ) {//  Loop over interaction list
        CJ = *L;                                                //   Set source cell
        rangeHost.push_back(sourceBegin[CJ]);                   //    Set begin index of leafs in source cell
        rangeHost.push_back(sourceSize[CJ]);                    //    Set number of leafs in source cell
      }                                                         //  End loop over interaction list
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
  P2P();                                                        // Evaluate P2P kernel
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {              // Loop over target cells
    if( !listP2P[CI-CI0].empty() ) {                            //  If the interation list is not empty
      BI0 = CI->LEAF;                                           //   Set target bodies begin iterator
      BIN = CI->LEAF + CI->NLEAF;                               //   Set target bodies end iterator
      int begin = targetBegin[CI];                              //   Offset of target leafs
      for( B_iter B=BI0; B!=BIN; ++B ) {                        //   Loop over target bodies
        B->pot += targetHost[4*(begin+B-BI0)+3];                //    Copy target values from GPU buffer
      }                                                         //   End loop over target bodies
      listP2P[CI-CI0].clear();                                  //   Clear P2P interaction list
    }                                                           //  End if for empty interation list
  }                                                             // End loop over target cells
  keysHost.clear();                                             // Clear keys vector
  rangeHost.clear();                                            // Clear range vector
  targetHost.clear();                                           // Clear target vector
  sourceHost.clear();                                           // Clear source vector
}

void Evaluator::evalL2L(Cells &cells) {                         // Evaluate L2L
  CI0 = cells.begin();                                          // Set begin iterator
  for( CI=cells.end()-2; CI!=cells.begin()-1; --CI ) {          // Loop over cells topdown (except root cell)
    CJ = CI0 + CI->PARENT;                                      //  Set source cell iterator
    L2L();                                                      //  Evaluate L2L kernel
  }                                                             // End loop over cells topdown
}

void Evaluator::evalL2P(Cells &cells) {                         // Evaluate L2P
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {              // Loop over cells
    if( CI->NCHILD == 0 ) L2P();                                //  If cell is a twig evaluate L2P kernel
  }                                                             // End loop over cells topdown
}
