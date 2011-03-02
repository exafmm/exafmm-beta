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
  Map sourceBegin;                                              // Define map for offset of source cells in listP2M
  Map sourceSize;                                               // Define map for size of source cells in listP2M
  for( CJ=cells.begin(); CJ!=cells.end(); ++CJ ) {              // Loop over target cells
    CJ->M = CJ->L = 0;                                          //  Initialize multipole & local coefficients
    sourceBegin[CJ] = sourceHost.size() / 4;                    //  Key : iterator, Value : offset of source leafs
    sourceSize[CJ] = CJ->NLEAF;                                 //   Key : iterator, Value : number of leafs
    for( B_iter B=CJ->LEAF; B!=CJ->LEAF+CJ->NLEAF; ++B ) {      //  Loop over leafs in source cell
      sourceHost.push_back(B->pos[0]);                          //   Copy x position to GPU buffer
      sourceHost.push_back(B->pos[1]);                          //   Copy y position to GPU buffer
      sourceHost.push_back(B->pos[2]);                          //   Copy z position to GPU buffer
      sourceHost.push_back(B->scal);                            //   Copy mass/charge to GPU buffer
    }                                                           //  End loop over leafs
  }                                                             // End loop over source map
  Map targetBegin;                                              // Define map for offset of target cells
  int key = 0;                                                  // Initialize key to range of leafs in source cells
  for( CJ=cells.begin(); CJ!=cells.end(); ++CJ ) {              // Loop over target cells
    keysHost.push_back(key);                                    //  Save key to range of coefs in source cells
    key += 3;                                                   //  Increment key counter
    rangeHost.push_back(1);                                     //  Save size of interaction list
    rangeHost.push_back(sourceBegin[CJ]);                       //   Set begin index of coefs in source cell
    rangeHost.push_back(sourceSize[CJ]);                        //   Set number of coefs in source cell
    targetBegin[CJ] = targetHost.size();                        //  Key : iterator, Value : offset of target coefs
    targetHost.push_back(CJ->X[0]);                             //  Copy x position to GPU buffer
    targetHost.push_back(CJ->X[1]);                             //  Copy y position to GPU buffer
    targetHost.push_back(CJ->X[2]);                             //  Copy z position to GPU buffer
    for( int i=0; i!=2*NCOEF; ++i ) {                           //  Loop over coefs in target cell
      targetHost.push_back(0);                                  //   Pad GPU buffer
    }                                                           //  End loop over coefs
    int numPad = 2 * THREADS - 2 * NCOEF - 3;                   //  Number of elements to pad in target GPU buffer
    assert(numPad >= 0);                                        //  THREADS must be large enough
    for( int i=0; i!=numPad; ++i ) {                            //  Loop over elements to pad
      targetHost.push_back(0);                                  //   Pad GPU buffer
    }                                                           //  End loop over elements to pad
  }                                                             // End loop over target cells
  P2M();                                                        // Evaluate P2M kernel
  for( CJ=cells.begin(); CJ!=cells.end(); ++CJ ) {              // Loop over target cells
    int begin = targetBegin[CJ];                                //  Offset of target coefs
    for( int i=0; i!=NCOEF; ++i ) {                             //  Loop over target coefs
      CJ->M[i].real() += targetHost[begin+2*i+0];               //   Copy real target values from GPU buffer
      CJ->M[i].imag() += targetHost[begin+2*i+1];               //   Copy imaginary target values from GPU buffer
    }                                                           //  End loop over target coefs
  }                                                             // End loop over target cells
  keysHost.clear();                                             // Clear keys vector
  rangeHost.clear();                                            // Clear range vector
  targetHost.clear();                                           // Clear target vector
  sourceHost.clear();                                           // Clear source vector
}

void Evaluator::evalM2M(Cells &cells) {                         // Evaluate M2M
  CI0 = cells.begin();                                          // Set begin iterator for target
  CJ0 = cells.begin();                                          // Set begin iterator for source
  Lists listM2M(cells.size());                                  // Define M2M interation list vector
  Map sourceBegin;                                              // Define map for offset of source cells in listM2M
  Map sourceSize;                                               // Define map for size of source cells in listM2M
  int level = getLevel(CJ0->I);                                 // Level of twig
  while( level != 0 ) {                                         // While level of source is not root level
    for( CJ=cells.begin(); CJ!=cells.end(); ++CJ ) {            //  Loop over cells bottomup (except root cell)
      if( getLevel(CJ->I) == level ) {                          //   If source cell is at current level
        CI = CJ0 + CJ->PARENT;                                  //    Set target cell iterator
        listM2M[CI-CI0].push_back(CJ);                          //    Push source cell into M2M interaction list
        sourceBegin[CJ] = sourceHost.size();                    //    Key : iterator, Value : offset of sources
        sourceSize[CJ] = 2 * NCOEF;                             //    Key : iterator, Value : number of coefs
        sourceHost.push_back(CJ->X[0]);                         //    Copy x position to GPU buffer
        sourceHost.push_back(CJ->X[1]);                         //    Copy y position to GPU buffer
        sourceHost.push_back(CJ->X[2]);                         //    Copy z position to GPU buffer
        for( int i=0; i!=NCOEF; ++i ) {                         //    Loop over coefs in source cell
          sourceHost.push_back((CJ->M[i]).real());              //     Copy real multipole to GPU buffer
          sourceHost.push_back((CJ->M[i]).imag());              //     Copy imaginary multipole to GPU buffer
        }                                                       //    End loop over coefs
      }                                                         //   Endif for current level
    }                                                           //  End loop over cells
    Map targetBegin;                                            //  Define map for offset of target cells
    int key = 0;                                                //  Initialize key to range of coefs in source cells
    for( CI=cells.begin(); CI!=cells.end(); ++CI ) {            //  Loop over cells bottomup (except root cell)
      if( !listM2M[CI-CI0].empty() ) {                          //   If the interation list is not empty
        keysHost.push_back(key);                                //    Save key to range of coefs in source cells
        key += 2*listM2M[CI-CI0].size()+1;                      //    Increment key counter
        rangeHost.push_back(listM2M[CI-CI0].size());            //    Save size of interaction list
        for( L_iter L=listM2M[CI-CI0].begin(); L!=listM2M[CI-CI0].end(); ++L ) {//  Loop over interaction list
          CJ = *L;                                              //    Set source cell
          rangeHost.push_back(sourceBegin[CJ]);                 //     Set begin index of coefs in source cell
          rangeHost.push_back(sourceSize[CJ]);                  //     Set number of coefs in source cell
        }                                                       //    End loop over interaction list
        targetBegin[CI] = targetHost.size();                    //    Key : iterator, Value : offset of target coefs
        targetHost.push_back(CI->X[0]);                         //    Copy x position to GPU buffer
        targetHost.push_back(CI->X[1]);                         //    Copy y position to GPU buffer
        targetHost.push_back(CI->X[2]);                         //    Copy z position to GPU buffer
        for( int i=0; i!=2*NCOEF; ++i ) {                       //    Loop over coefs in target cell
          targetHost.push_back(0);                              //     Pad GPU buffer
        }                                                       //    End loop over coefs
        int numPad = 2 * THREADS - 2 * NCOEF - 3;               //    Number of elements to pad in target GPU buffer
        assert(numPad >= 0);                                    //    THREADS must be large enough
        for( int i=0; i!=numPad; ++i ) {                        //    Loop over elements to pad
          targetHost.push_back(0);                              //     Pad GPU buffer
        }                                                       //    End loop over elements to pad
      }                                                         //   End if for empty interation list
    }                                                           //  End loop over cells
    M2M();                                                      //  Evaluate M2M kernel
    for( CI=cells.begin(); CI!=cells.end(); ++CI ) {            //  Loop over target cells
      if( !listM2M[CI-CI0].empty() ) {                          //   If the interation list is not empty
        int begin = targetBegin[CI];                            //    Offset of target coefs
        for( int i=0; i!=NCOEF; ++i ) {                         //    Loop over target coefs
          CI->M[i].real() += targetHost[begin+2*i+0];           //     Copy real target values from GPU buffer
          CI->M[i].imag() += targetHost[begin+2*i+1];           //     Copy imaginary target values from GPU buffer
        }                                                       //    End loop over target coefs
        listM2M[CI-CI0].clear();                                //    Clear M2M interaction list
      }                                                         //   End if for empty interation list
    }                                                           //  End loop over target cells
    keysHost.clear();                                           //  Clear keys vector
    rangeHost.clear();                                          //  Clear range vector
    targetHost.clear();                                         //  Clear target vector
    sourceHost.clear();                                         //  Clear source vector
    level--;                                                    //  Decrement level
  }                                                             // End while loop over levels
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
  Map sourceSize;                                               // Define map for size of source cells in listM2P
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {              // Loop over target cells
    for( L_iter L=listM2P[CI-CI0].begin(); L!=listM2P[CI-CI0].end(); ++L ) {//  Loop over interaction list
      CJ = *L;                                                  //   Set source cell
      sourceSize[CJ] = 2 * NCOEF;                               //   Key : iterator, Value : number of coefs
    }                                                           //  End loop over interaction list
  }                                                             // End loop over target cells
  Map sourceBegin;                                              // Define map for offset of source cells in listM2P
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
    if( !listM2P[CI-CI0].empty() ) {                            //  If the interation list is not empty
      BI0 = CI->LEAF;                                           //   Set target bodies begin iterator
      BIN = CI->LEAF + CI->NLEAF;                               //   Set target bodies end iterator
      int blocks = (BIN - BI0 - 1) / THREADS + 1;               //   Number of thread blocks needed for this target cell
      for( int i=0; i!=blocks; ++i ) {                          //   Loop over thread blocks
        keysHost.push_back(key);                                //    Save key to range of leafs in source cells
      }                                                         //   End loop over thread blocks
      key += 2*listM2P[CI-CI0].size()+1;                        //   Increment key counter
      rangeHost.push_back(listM2P[CI-CI0].size());              //   Save size of interaction list
      for( L_iter L=listM2P[CI-CI0].begin(); L!=listM2P[CI-CI0].end(); ++L ) {//  Loop over interaction list
        CJ = *L;                                                //   Set source cell
        rangeHost.push_back(sourceBegin[CJ]);                   //    Set begin index of coefs in source cell
        rangeHost.push_back(sourceSize[CJ]);                    //    Set number of coefs in source cell
      }                                                         //   End loop over interaction list
      targetBegin[CI] = targetHost.size();                      //   Key : iterator, Value : offset of target coefs
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
  M2P();                                                        // Evaluate M2P kernel
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {              // Loop over target cells
    if( !listM2P[CI-CI0].empty() ) {                            //  If the interation list is not empty
      BI0 = CI->LEAF;                                           //   Set target bodies begin iterator
      BIN = CI->LEAF + CI->NLEAF;                               //   Set target bodies end iterator
      int begin = targetBegin[CI];                              //   Offset of target leafs
      for( B_iter B=BI0; B!=BIN; ++B ) {                        //   Loop over target bodies
        B->pot += targetHost[4*(begin+B-BI0)+3];                //    Copy target values from GPU buffer
      }                                                         //   End loop over target bodies
      listM2P[CI-CI0].clear();                                  //   Clear M2P interaction list
    }                                                           //  End if for empty interation list
  }                                                             // End loop over target cells
  keysHost.clear();                                             // Clear keys vector
  rangeHost.clear();                                            // Clear range vector
  targetHost.clear();                                           // Clear target vector
  sourceHost.clear();                                           // Clear source vector
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
  Lists listL2L(cells.size());                                  // Define L2L interation list vector
  Map sourceBegin;                                              // Define map for offset of source cells in listL2L
  Map sourceSize;                                               // Define map for size of source cells in listL2L
  int maxLevel = getLevel(CI0->I);                              // Level of twig
  int level = 1;                                                // Start level from 1
  while( level != maxLevel+1 ) {                                // While level of source is not root level
    for( CI=cells.end()-2; CI!=cells.begin()-1; --CI ) {        //  Loop over cells topdown (except root cell)
      if( getLevel(CI->I) == level ) {                          //   If target cell is at current level
        CJ = CI0 + CI->PARENT;                                  //    Set source cell iterator
        listL2L[CI-CI0].push_back(CJ);                          //    Push source cell into L2L interaction list
        if( sourceBegin[CJ] == 0 ) {                            //    If the source cell has not been stored yet
          sourceBegin[CJ] = sourceHost.size();                  //     Key : iterator, Value : offset of sources
          sourceSize[CJ] = 2 * NCOEF;                           //     Key : iterator, Value : number of coefs
          sourceHost.push_back(CJ->X[0]);                       //     Copy x position to GPU buffer
          sourceHost.push_back(CJ->X[1]);                       //     Copy y position to GPU buffer
          sourceHost.push_back(CJ->X[2]);                       //     Copy z position to GPU buffer
          for( int i=0; i!=NCOEF; ++i ) {                       //     Loop over coefs in source cell
            sourceHost.push_back((CJ->L[i]).real());            //      Copy real local to GPU buffer
            sourceHost.push_back((CJ->L[i]).imag());            //      Copy imaginary local to GPU buffer
          }                                                     //     End loop over coefs
        }                                                       //    Endif for current level
      }                                                         //   Endif for stored source cell
    }                                                           //  End loop over cells topdown
    Map targetBegin;                                            //  Define map for offset of target cells
    int key = 0;                                                //  Initialize key to range of coefs in source cells
    for( CI=cells.end()-2; CI!=cells.begin()-1; --CI ) {        //  Loop over cells topdown (except root cell)
      if( !listL2L[CI-CI0].empty() ) {                          //   If the interation list is not empty
        keysHost.push_back(key);                                //    Save key to range of coefs in source cells
        key += 2*listL2L[CI-CI0].size()+1;                      //    Increment key counter
        rangeHost.push_back(listL2L[CI-CI0].size());            //    Save size of interaction list
        for( L_iter L=listL2L[CI-CI0].begin(); L!=listL2L[CI-CI0].end(); ++L ) {//  Loop over interaction list
          CJ = *L;                                              //    Set source cell
          rangeHost.push_back(sourceBegin[CJ]);                 //     Set begin index of coefs in source cell
          rangeHost.push_back(sourceSize[CJ]);                  //     Set number of coefs in source cell
        }                                                       //    End loop over interaction list
        targetBegin[CI] = targetHost.size();                    //    Key : iterator, Value : offset of target coefs
        targetHost.push_back(CI->X[0]);                         //    Copy x position to GPU buffer
        targetHost.push_back(CI->X[1]);                         //    Copy y position to GPU buffer
        targetHost.push_back(CI->X[2]);                         //    Copy z position to GPU buffer
        for( int i=0; i!=2*NCOEF; ++i ) {                       //    Loop over coefs in target cell
          targetHost.push_back(0);                              //     Pad GPU buffer
        }                                                       //    End loop over coefs
        int numPad = 2 * THREADS - 2 * NCOEF - 3;               //    Number of elements to pad in target GPU buffer
        assert(numPad >= 0);                                    //    THREADS must be large enough
        for( int i=0; i!=numPad; ++i ) {                        //    Loop over elements to pad
          targetHost.push_back(0);                              //     Pad GPU buffer
        }                                                       //    End loop over elements to pad
      }                                                         //   End if for empty interation list
    }                                                           //  End loop over cells topdown
    L2L();                                                      //  Evaluate L2L kernel
    for( CI=cells.end()-2; CI!=cells.begin()-1; --CI ) {        //  Loop over cells topdown (except root cell)
      if( !listL2L[CI-CI0].empty() ) {                          //   If the interation list is not empty
        int begin = targetBegin[CI];                            //    Offset of target coefs
        for( int i=0; i!=NCOEF; ++i ) {                         //    Loop over target coefs
          CI->L[i].real() += targetHost[begin+2*i+0];           //     Copy real target values from GPU buffer
          CI->L[i].imag() += targetHost[begin+2*i+1];           //     Copy imaginary target values from GPU buffer
        }                                                       //    End loop over target coefs
        listL2L[CI-CI0].clear();                                //    Clear L2L interaction list
      }                                                         //   End if for empty interation list
    }                                                           //  End loop over cells topdown
    keysHost.clear();                                           //  Clear keys vector
    rangeHost.clear();                                          //  Clear range vector
    targetHost.clear();                                         //  Clear target vector
    sourceHost.clear();                                         //  Clear source vector
    level++;                                                    //  Increment level
  }                                                             // End while loop over levels
}

void Evaluator::evalL2P(Cells &cells) {                         // Evaluate L2P
  Map sourceBegin;                                              // Define map for offset of source cells in listL2P
  Map sourceSize;                                               // Define map for size of source cells in listL2P
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {              // Loop over cells
    if( CI->NCHILD == 0 ) {                                     //  If cell is a twig evaluate L2P kernel
      sourceBegin[CI] = sourceHost.size();                      //   Key : iterator, Value : offset of sources
      sourceSize[CI] = 2 * NCOEF;                               //   Key : iterator, Value : number of coefs
      sourceHost.push_back(CI->X[0]);                           //   Copy x position to GPU buffer
      sourceHost.push_back(CI->X[1]);                           //   Copy y position to GPU buffer
      sourceHost.push_back(CI->X[2]);                           //   Copy z position to GPU buffer
      for( int i=0; i!=NCOEF; ++i ) {                           //   Loop over coefs in source cell
        sourceHost.push_back((CI->L[i]).real());                //    Copy real multipole to GPU buffer
        sourceHost.push_back((CI->L[i]).imag());                //    Copy imaginary multipole to GPU buffer
      }                                                         //   End loop over coefs
    }                                                           //  Endif for twig cells
  }                                                             // End loop over cells topdown
  Map targetBegin;                                              // Define map for offset of target cells
  int key = 0;                                                  // Initialize key to range of coefs in source cells
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {              // Loop over target cells
    if( CI->NCHILD == 0 ) {                                     //  If cell is a twig evaluate L2P kernel
      BI0 = CI->LEAF;                                           //   Set target bodies begin iterator
      BIN = CI->LEAF + CI->NLEAF;                               //   Set target bodies end iterator
      int blocks = (BIN - BI0 - 1) / THREADS + 1;               //   Number of thread blocks needed for this target cell
      for( int i=0; i!=blocks; ++i ) {                          //   Loop over thread blocks
        keysHost.push_back(key);                                //    Save key to range of leafs in source cells
      }                                                         //   End loop over thread blocks
      key += 3;                                                 //   Increment key counter
      rangeHost.push_back(1);                                   //   Save size of interaction list
      rangeHost.push_back(sourceBegin[CI]);                     //   Set begin index of coefs in source cell
      rangeHost.push_back(sourceSize[CI]);                      //   Set number of coefs in source cell
      targetBegin[CI] = targetHost.size() / 4;                  //   Key : iterator, Value : offset of target coefs
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
    }                                                           //  Endif for twig cells
  }                                                             // End loop over target cells
  L2P();                                                        // Evaluate L2P kernel
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {              // Loop over target cells
    if( CI->NCHILD == 0 ) {                                     //  If cell is a twig evaluate L2P kernel
      BI0 = CI->LEAF;                                           //   Set target bodies begin iterator
      BIN = CI->LEAF + CI->NLEAF;                               //   Set target bodies end iterator
      int begin = targetBegin[CI];                              //   Offset of target leafs
      for( B_iter B=BI0; B!=BIN; ++B ) {                        //   Loop over target bodies
        B->pot += targetHost[4*(begin+B-BI0)+3];                //    Copy target values from GPU buffer
      }                                                         //   End loop over target bodies
    }                                                           //  Endif for twig cells
  }                                                             // End loop over target cells
  keysHost.clear();                                             // Clear keys vector
  rangeHost.clear();                                            // Clear range vector
  targetHost.clear();                                           // Clear target vector
  sourceHost.clear();                                           // Clear source vector
}
