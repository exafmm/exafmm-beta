void Evaluator::setSourceBody() {                               // Set source buffer for bodies
  startTimer("Set sourceB  ");                                  // Start timer
  for( M_iter M=sourceSize.begin(); M!=sourceSize.end(); ++M ) {// Loop over source map
    CJ = M->first;                                              //  Set source cell
    sourceBegin[CJ] = sourceHost.size() / 7;                    //  Key : iterator, Value : offset of source leafs
    for( B_iter B=CJ->LEAF; B!=CJ->LEAF+CJ->NLEAF; ++B ) {      //  Loop over leafs in source cell
      sourceHost.push_back(B->X[0]);                            //   Copy x position to GPU buffer
      sourceHost.push_back(B->X[1]);                            //   Copy y position to GPU buffer
      sourceHost.push_back(B->X[2]);                            //   Copy z position to GPU buffer
      sourceHost.push_back(B->SRC[0]);                          //   Copy 1st source value to GPU buffer
      sourceHost.push_back(B->SRC[1]);                          //   Copy 2nd source value to GPU buffer
      sourceHost.push_back(B->SRC[2]);                          //   Copy 3rd source value to GPU buffer
      sourceHost.push_back(B->SRC[3]);                          //   Copy 4th source value to GPU buffer
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
      }                                                         //   End loop over coefs
    } else {                                                    //  If source is L
      for( int i=0; i!=NCOEF; ++i ) {                           //   Loop over coefs in source cell
        sourceHost.push_back((CJ->L[i]).real());                //    Copy real multipole to GPU buffer
        sourceHost.push_back((CJ->L[i]).imag());                //    Copy imaginary multipole to GPU buffer
      }                                                         //   End loop over coefs
    }                                                           //  Endif for source type
  }                                                             // End loop over source map
  stopTimer("Set sourceC  ");                                   // Stop timer
}

void Evaluator::setTargetBody(Lists lists, Maps flags) {        // Set target buffer for bodies
  startTimer("Set targetB  ");                                  // Start timer
  int key = 0;                                                  // Initialize key to range of coefs in source cells
  for( CI=CIB; CI!=CIE; ++CI ) {                                // Loop over target cells
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
      targetBegin[CI] = targetHost.size() / 6;                  //   Key : iterator, Value : offset of target leafs
      for( B_iter B=BI0; B!=BIN; ++B ) {                        //   Loop over leafs in target cell
        targetHost.push_back(B->X[0]);                          //    Copy x position to GPU buffer
        targetHost.push_back(B->X[1]);                          //    Copy y position to GPU buffer
        targetHost.push_back(B->X[2]);                          //    Copy z position to GPU buffer
        targetHost.push_back(B->SRC[0]);                        //    Copy 1st target value to GPU buffer
        targetHost.push_back(B->SRC[1]);                        //    Copy 2nd target value to GPU buffer
        targetHost.push_back(B->SRC[2]);                        //    Copy 3rd target value to GPU buffer
      }                                                         //   End loop over leafs
      int numPad = blocks * THREADS - (BIN - BI0);              //   Number of elements to pad in target GPU buffer
      for( int i=0; i!=numPad; ++i ) {                          //   Loop over elements to pad
        targetHost.push_back(0);                                //    Pad x position in GPU buffer
        targetHost.push_back(0);                                //    Pad y position in GPU buffer
        targetHost.push_back(0);                                //    Pad z position in GPU buffer
        targetHost.push_back(0);                                //    Pad 1st target value to GPU buffer
        targetHost.push_back(0);                                //    Pad 2nd target value to GPU buffer
        targetHost.push_back(0);                                //    Pad 3rd target value to GPU buffer
      }                                                         //   End loop over elements to pad
    }                                                           //  End if for empty interation list
  }                                                             // End loop over target cells
  stopTimer("Set targetB  ");                                   // Stop timer
}

void Evaluator::setTargetCell(Lists lists, Maps flags) {        // Set target buffer for cells
  startTimer("Set targetC  ");                                  // Start timer
  int key = 0;                                                  // Initialize key to range of coefs in target cells
  for( CI=CIB; CI!=CIE; ++CI ) {                                // Loop over target cells
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
      int numPad = 2 * THREADS * NCOEF / NTERM - 2 * NCOEF - 3; //   Number of elements to pad in target GPU buffer
      assert(numPad >= 0);                                      //   THREADS must be large enough
      for( int i=0; i!=numPad; ++i ) {                          //   Loop over elements to pad
        targetHost.push_back(0);                                //    Pad GPU buffer
      }                                                         //   End loop over elements to pad
    }                                                           //  End if for empty interation list
  }                                                             // End loop over target cells
  stopTimer("Set targetC  ");                                   // Stop timer
}

void Evaluator::getTargetBody(Lists &lists) {                   // Get body values from target buffer
  startTimer("Get targetB  ");                                  // Start timer
  for( CI=CIB; CI!=CIE; ++CI ) {                                // Loop over target cells
    if( !lists[CI-CI0].empty() ) {                              //  If the interation list is not empty
      BI0 = CI->LEAF;                                           //   Set target bodies begin iterator
      BIN = CI->LEAF + CI->NLEAF;                               //   Set target bodies end iterator
      int begin = targetBegin[CI];                              //   Offset of target leafs
      if( kernelName == Gaussian ) {                            //   If Gaussian kernel
        for( B_iter B=BI0; B!=BIN; ++B ) {                      //    Loop over target bodies
          B->TRG[0] += targetHost[6*(begin+B-BI0)+0];           //     Copy 1st target value from GPU buffer
        }                                                       //    End loop over target bodies
      } else {                                                  //   If not Gaussian kernel
        for( B_iter B=BI0; B!=BIN; ++B ) {                      //    Loop over target bodies
          B->TRG[0] += targetHost[6*(begin+B-BI0)+0];           //     Copy 1st target value from GPU buffer
          B->TRG[1] += targetHost[6*(begin+B-BI0)+1];           //     Copy 2nd target value from GPU buffer
          B->TRG[2] += targetHost[6*(begin+B-BI0)+2];           //     Copy 3rd target value from GPU buffer
          B->TRG[3] += targetHost[6*(begin+B-BI0)+3];           //     Copy 4th target value from GPU buffer
        }                                                       //    End loop over target bodies
      }                                                         //   Endif for Gaussian kernel
      lists[CI-CI0].clear();                                    //   Clear interaction list
    }                                                           //  End if for empty interation list
  }                                                             // End loop over target cells
  stopTimer("Get targetB  ");                                   // Stop timer
}

void Evaluator::getTargetCell(Lists &lists, bool isM=true) {    // Get body values from target buffer
  startTimer("Get targetC  ");                                  // Start timer
  for( CI=CIB; CI!=CIE; ++CI ) {                                // Loop over target cells
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

void Evaluator::testMACP2P(C_iter Ci, C_iter Cj) {              // Test multipole acceptance criteria for P2P kernel
  listP2P[Ci-CI0].push_back(Cj);                                // Push source cell into P2P interaction list
  flagP2P[Ci-CI0][Cj] |= Iperiodic;                             // Flip bit of periodic image flag
  NP2P++;                                                       // Count P2P kernel execution
}

void Evaluator::testMACM2L(C_iter Ci, C_iter Cj) {              // Test multipole acceptance criteria for M2L kernel
  vect dist = Ci->X - Cj->X - Xperiodic;                        // Distance vector between cells
  real R = std::sqrt(norm(dist));                               // Distance between cells
  if( Ci->R + Cj->R > THETA*R ) {                               // If cell is too large
    Pair pair(Ci,Cj);                                           //  Form pair of interacting cells
    pairs.push(pair);                                           //  Push interacting pair into stack
  } else {                                                      // If cell is small enough
    listM2L[Ci-CI0].push_back(Cj);                              //  Push source cell into M2L interaction list
    flagM2L[Ci-CI0][Cj] |= Iperiodic;                           //  Flip bit of periodic image flag
    NM2L++;                                                     //  Count M2L kernel execution
  }                                                             // Endif for interaction
}

void Evaluator::testMACM2P(C_iter Ci, C_iter Cj) {              // Test multipole acceptance criteria for M2P kernel
  vect dist = Ci->X - Cj->X - Xperiodic;                        // Distance vector between cells
  real R = std::sqrt(norm(dist));                               // Distance between cells
  if( Ci->NCHILD != 0 || Ci->R + Cj->R > THETA*R ) {            // If target is not twig or cell is too large
    Pair pair(Ci,Cj);                                           //  Form pair of interacting cells
    pairs.push(pair);                                           //  Push interacting pair into stack
  } else {                                                      // If target is twig and cell is small enough
    listM2P[Ci-CI0].push_back(Cj);                              //  Push source cell into M2P interaction list
    flagM2P[Ci-CI0][Cj] |= Iperiodic;                           //  Flip bit of periodic image flag
    NM2P++;                                                     //  Count M2P kernel execution
  }                                                             // Endif for interaction
} 

void Evaluator::timeKernels() {                                 // Time all kernels for auto-tuning
  Bodies ibodies(NCRIT), jbodies(NCRIT);                        // Artificial bodies
  for( B_iter Bi=ibodies.begin(),Bj=jbodies.begin(); Bi!=ibodies.end(); ++Bi, ++Bj ) {// Loop over artificial bodies
    Bi->X = 0;                                                  //  Set coordinates of target body
    Bj->X = 1;                                                  //  Set coordinates of source body
  }                                                             // End loop over artificial bodies
  Cells icells, jcells;                                         // Artificial cells
  icells.resize(100);                                           // 100 artificial target cells
  jcells.resize(100);                                           // 100 artificial source cells
  CI0 = icells.begin();                                         // Set global begin iterator for source
  for( C_iter Ci=icells.begin(); Ci!=icells.end(); ++Ci ) {     // Loop over target cells
    Ci->X = 0;                                                  //  Set coordinates of target cell
    Ci->NLEAF = NCRIT;                                          //  Number of leafs in target cell
    Ci->LEAF = ibodies.begin();                                 //  Leaf iterator in target cell
  }                                                             // End loop over target cells
  for( C_iter Cj=jcells.begin(); Cj!=jcells.end(); ++Cj ) {     // Loop over source cells
    Cj->X = 1;                                                  //  Set coordinates of source cell
    Cj->NLEAF = NCRIT;                                          //  Number of leafs in source cell
    Cj->LEAF = jbodies.begin();                                 //  Leaf iterator in source cell
  }                                                             // End loop over source cells
  listM2L.resize(icells.size());                                // Resize M2L interaction list
  listM2P.resize(icells.size());                                // Resize M2P interaction list
  listP2P.resize(icells.size());                                // Resize P2P interaction list
  flagM2L.resize(icells.size());                                // Resize M2L periodic image flag
  flagM2P.resize(icells.size());                                // Resize M2P periodic image flag
  flagP2P.resize(icells.size());                                // Resize P2P periodic image flag
  for( C_iter Ci=icells.begin(); Ci!=icells.end(); ++Ci ) {     // Loop over target cells
    for( C_iter Cj=jcells.begin(); Cj!=jcells.end(); ++Cj ) {   //  Loop over source cells
      listP2P[Ci-CI0].push_back(Cj);                            //   Push source cell into P2P interaction list
      listM2P[Ci-CI0].push_back(Cj);                            //   Push source cell into P2P interaction list
      listM2L[Ci-CI0].push_back(Cj);                            //   Push source cell into P2P interaction list
    }                                                           //  End loop over source cells
  }                                                             // End loop over target cells
  startTimer("P2P kernel   ");                                  // Start timer
  evalP2P(icells);                                              // Evaluate P2P kernel
  timeP2P = stopTimer("P2P kernel   ") / NCRIT / NCRIT;         // Stop timer
  startTimer("M2L kernel   ");                                  // Start timer
  evalM2L(icells);                                              // Evaluate M2L kernel
  timeM2L = stopTimer("M2L kernel   ");                         // Stop timer
  startTimer("M2P kernel   ");                                  // Start timer
  evalM2P(icells);                                              // Evaluate M2P kernel
  timeM2P = stopTimer("M2P kernel   ") / NCRIT;                 // Stop timer
}

void Evaluator::traversePeriodic(Cells &cells, Cells &jcells, int method) {// Traverse tree for periodic cells
  C_iter Cj = jcells.end()-1;                                   // Initialize iterator for periodic source cell
  for( int level=0; level<IMAGES-1; ++level ) {                 // Loop over sublevels of tree
    for( int I=0; I!=26; ++I, --Cj ) {                          //  Loop over periodic images (exclude center)
      switch (method) {                                         //   Switch between method
      case 0 :                                                  //   0 : treecode
        for( C_iter Ci=cells.begin(); Ci!=cells.end(); ++Ci ) { //   Loop over cells
          if( Ci->NCHILD == 0 ) {                               //     If cell is twig
            listM2P[Ci-CI0].push_back(Cj);                      //      Push source cell into M2P interaction list
            flagM2P[Ci-CI0][Cj] = Icenter;                      //      Flip bit of periodic image flag
          }                                                     //     Endif for twig
        }                                                       //    End loop over cells
        break;                                                  //    Terminate this case
      case 1 :                                                  //   1 : FMM
      case 2 :                                                  //   2 : hybrid
        C_iter Ci = cells.end() - 1;                            //    Set root cell as target
        listM2L[Ci-CI0].push_back(Cj);                          //    Push source cell into M2L interaction list
        flagM2L[Ci-CI0][Cj] = Icenter;                          //    Flip bit of periodic image flag
        break;                                                  //    Terminate this case
      }                                                         //   End switch between methods
    }                                                           //  End loop over x periodic direction
  }                                                             // End loop over sublevels of tree
}

void Evaluator::evalP2P(Bodies &ibodies, Bodies &jbodies, bool onCPU) {// Evaluate P2P
  int numIcall = int(ibodies.size()-1)/MAXBODY+1;               // Number of icall loops
  int numJcall = int(jbodies.size()-1)/MAXBODY+1;               // Number of jcall loops
  int ioffset = 0;                                              // Initialzie offset for icall loops
  for( int icall=0; icall!=numIcall; ++icall ) {                // Loop over icall
    BI0 = ibodies.begin()+ioffset;                              //  Set target bodies begin iterator
    BIN = ibodies.begin()+std::min(ioffset+MAXBODY,int(ibodies.size()));// Set target bodies end iterator
    int joffset = 0;                                            //  Initialize offset for jcall loops
    for( int jcall=0; jcall!=numJcall; ++jcall ) {              //  Loop over jcall
      BJ0 = jbodies.begin()+joffset;                            //  Set source bodies begin iterator
      BJN = jbodies.begin()+std::min(joffset+MAXBODY,int(jbodies.size()));// Set source bodies end iterator
      if( onCPU ) {                                             //  If calculation is to be done on CPU
        Xperiodic = 0;                                          //   Set periodic coordinate offset
        selectP2P_CPU();                                        //   Select P2P_CPU kernel
      } else {                                                  //  If calculation is to be done on GPU
        constHost.push_back(2*R0);                              //   Copy domain size to GPU buffer
        for( B_iter B=BJ0; B!=BJN; ++B ) {                      //   Loop over source bodies
          sourceHost.push_back(B->X[0]);                        //   Copy x position to GPU buffer
          sourceHost.push_back(B->X[1]);                        //   Copy y position to GPU buffer
          sourceHost.push_back(B->X[2]);                        //   Copy z position to GPU buffer
          sourceHost.push_back(B->SRC[0]);                      //   Copy 1st source value to GPU buffer
          sourceHost.push_back(B->SRC[1]);                      //   Copy 2nd source value to GPU buffer
          sourceHost.push_back(B->SRC[2]);                      //   Copy 3rd source value to GPU buffer
          sourceHost.push_back(B->SRC[3]);                      //   Copy 4th source value to GPU buffer
        }                                                       //   End loop over source bodies
        int key = 0;                                            //   Initialize key to range of leafs in source cells
        int blocks = (BIN - BI0 - 1) / THREADS + 1;             //   Number of thread blocks needed for this target cell
        for( int i=0; i!=blocks; ++i ) {                        //   Loop over thread blocks
          keysHost.push_back(key);                              //    Save key to range of leafs in source cells
        }                                                       //   End loop over thread blocks
        rangeHost.push_back(1);                                 //   Save size of interaction list
        rangeHost.push_back(0);                                 //   Set begin index of leafs
        rangeHost.push_back(BJN-BJ0);                           //   Set number of leafs
        rangeHost.push_back(Icenter);                           //   Set periodic image flag
        for( B_iter B=BI0; B!=BIN; ++B ) {                      //   Loop over target bodies
          targetHost.push_back(B->X[0]);                        //    Copy x position to GPU buffer
          targetHost.push_back(B->X[1]);                        //    Copy y position to GPU buffer
          targetHost.push_back(B->X[2]);                        //    Copy z position to GPU buffer
          targetHost.push_back(B->SRC[0]);                      //    Copy 1st target value to GPU buffer
          targetHost.push_back(B->SRC[1]);                      //    Copy 2nd target value to GPU buffer
          targetHost.push_back(B->SRC[2]);                      //    Copy 3rd target value to GPU buffer
        }                                                       //   End loop over target bodies
        int numPad = blocks * THREADS - (BIN - BI0);            //   Number of elements to pad in target GPU buffer
        for( int i=0; i!=numPad; ++i ) {                        //   Loop over elements to pad
          targetHost.push_back(0);                              //    Pad x position in GPU buffer
          targetHost.push_back(0);                              //    Pad y position in GPU buffer
          targetHost.push_back(0);                              //    Pad z position in GPU buffer
          targetHost.push_back(0);                              //    Pad 1st target value to GPU buffer
          targetHost.push_back(0);                              //    Pad 2nd target value to GPU buffer
          targetHost.push_back(0);                              //    Pad 3rd target value to GPU buffer
        }                                                       //   End loop over elements to pad
        selectP2P();                                            //   Select P2P kernel
        if( kernelName == Gaussian ) {                          //   If Gaussian kernel
          for( B_iter B=BI0; B!=BIN; ++B ) {                    //    Loop over target bodies
            B->TRG[0] += targetHost[6*(B-BI0)+0];               //     Copy 1st target value from GPU buffer
          }                                                     //    End loop over target bodies
        } else {                                                //   If not Gaussian kernel
          for( B_iter B=BI0; B!=BIN; ++B ) {                    //    Loop over target bodies
            B->TRG[0] += targetHost[6*(B-BI0)+0];               //     Copy 1st target value from GPU buffer
            B->TRG[1] += targetHost[6*(B-BI0)+1];               //     Copy 2nd target value from GPU buffer
            B->TRG[2] += targetHost[6*(B-BI0)+2];               //     Copy 3rd target value from GPU buffer
            B->TRG[3] += targetHost[6*(B-BI0)+3];               //     Copy 4th target value from GPU buffer
          }                                                     //    End loop over target bodies
        }                                                       //   Endif for Gaussian kernel
        keysHost.clear();                                       //   Clear keys vector
        rangeHost.clear();                                      //   Clear range vector
        constHost.clear();                                      //   Clear const vector
        targetHost.clear();                                     //   Clear target vector
        sourceHost.clear();                                     //   Clear source vector
      }                                                         //  Endif for CPU/GPU switch
      joffset += MAXBODY;                                       //  Increment jcall offset
    }                                                           // End loop over jcall
    ioffset += MAXBODY;                                         // Increment icall offset
  }                                                             // End loop over icall
}

void Evaluator::evalP2M(Cells &cells) {                         // Evaluate P2M
  CI0 = cells.begin();                                          // Set begin iterator for target
  const int numCell = MAXCELL/NCRIT/7;                          // Number of cells per icall
  int numIcall = int(cells.size()-1)/numCell+1;                 // Number of icall loops
  int ioffset = 0;                                              // Initialzie offset for icall loops
  for( int icall=0; icall!=numIcall; ++icall ) {                // Loop over icall
    CIB = cells.begin()+ioffset;                                //  Set begin iterator for target per call
    CIE = cells.begin()+std::min(ioffset+numCell,int(cells.size()));// Set end iterator for target per call
    constHost.push_back(2*R0);                                  //  Copy domain size to GPU buffer
    startTimer("Get list     ");                                //  Start timer
    Lists listP2M(cells.size());                                //  Define P2M interation list vector
    Maps  flagP2M(cells.size());                                //  Define P2M periodic image flag
    for( CI=CIB; CI!=CIE; ++CI ) {                              //  Loop over target cells
      CI->M = CI->L = 0;                                        //   Initialize multipole & local coefficients
      if( CI->NCHILD == 0 ) {                                   //   If cell is a twig
        listP2M[CI-CI0].push_back(CI);                          //    Push source cell into P2M interaction list
        flagP2M[CI-CI0][CI] |= Icenter;                         //    Flip bit of periodic image flag
        sourceSize[CI] = CI->NLEAF;                             //    Key : iterator, Value : number of leafs
      }                                                         //   End loop over cells topdown
    }                                                           //  End loop over source map
    stopTimer("Get list     ");                                 //  Stop timer
    setSourceBody();                                            //  Set source buffer for bodies
    setTargetCell(listP2M,flagP2M);                             //  Set target buffer for cells
    selectP2M();                                                //  Select P2M kernel
    getTargetCell(listP2M);                                     //  Get body values from target buffer
    clearBuffers();                                             //  Clear GPU buffers
    ioffset += numCell;                                         //  Increment ioffset
  }                                                             // End loop over icall
}

void Evaluator::evalM2M(Cells &cells) {                         // Evaluate M2M
  CI0 = cells.begin();                                          // Set begin iterator for target
  const int numCell = MAXCELL/NCOEF/2;                          // Number of cells per icall
  int numIcall = int(cells.size()-1)/numCell+1;                 // Number of icall loops
  int level = getLevel(CI0->ICELL);                             // Level of twig
  while( level != -1 ) {                                        // While level of target is not past root level
    int ioffset = 0;                                            //  Initialzie offset for icall loops
    for( int icall=0; icall!=numIcall; ++icall ) {              //  Loop over icall
      CIB = cells.begin()+ioffset;                              //   Set begin iterator for target per call
      CIE = cells.begin()+std::min(ioffset+numCell,int(cells.size()));// Set end iterator for target per call
      constHost.push_back(2*R0);                                //   Copy domain size to GPU buffer
      startTimer("Get list     ");                              //   Start timer
      Lists listM2M(cells.size());                              //   Define M2M interation list vector
      Maps  flagM2M(cells.size());                              //   Define M2M periodic image flag
      for( CI=CIB; CI!=CIE; ++CI ) {                            //   Loop over cells bottomup (except root cell)
        if( getLevel(CI->ICELL) == level ) {                    //    If target cell is at current level
          for( int i=0; i<CI->NCHILD; ++i ) {                   //     Loop over child cells
            CJ = CI0 + CI->CHILD[i];                            //      Set iterator for source cell
            listM2M[CI-CI0].push_back(CJ);                      //      Push source cell into M2M interaction list
            flagM2M[CI-CI0][CJ] |= Icenter;                     //      Flip bit of periodic image flag
            sourceSize[CJ] = 2 * NCOEF;                         //      Key : iterator, Value : number of coefs
          }                                                     //     End loop over child cells
        }                                                       //    Endif for current level
      }                                                         //   End loop over cells
      stopTimer("Get list     ");                               //   Stop timer
      setSourceCell();                                          //   Set source buffer for cells
      setTargetCell(listM2M,flagM2M);                           //   Set target buffer for cells
      selectM2M();                                              //   Select M2M kernel
      getTargetCell(listM2M);                                   //   Get body values from target buffer
      clearBuffers();                                           //   Clear GPU buffers
      ioffset += numCell;                                       //   Increment ioffset
    }                                                           //  End loop over icall
    level--;                                                    //  Decrement level
  }                                                             // End while loop over levels
}

void Evaluator::evalM2L(Cells &cells, bool kernel) {            // Evaluate M2L
  CI0 = cells.begin();                                          // Set begin iterator
  const int numCell = MAXCELL/NCOEF/2;                          // Number of cells per icall
  int numIcall = int(cells.size()-1)/numCell+1;                 // Number of icall loops
  int ioffset = 0 * kernel;                                     // Initialzie offset for icall loops
  for( int icall=0; icall!=numIcall; ++icall ) {                // Loop over icall
    CIB = cells.begin()+ioffset;                                //  Set begin iterator for target per call
    CIE = cells.begin()+std::min(ioffset+numCell,int(cells.size()));// Set end iterator for target per call
    constHost.push_back(2*R0);                                  //  Copy domain size to GPU buffer
    startTimer("Get list     ");                                //  Start timer
    for( CI=CIB; CI!=CIE; ++CI ) {                              //  Loop over target cells
      for( L_iter L=listM2L[CI-CI0].begin(); L!=listM2L[CI-CI0].end(); ++L ) {//  Loop over interaction list
        CJ = *L;                                                //    Set source cell
        sourceSize[CJ] = 2 * NCOEF;                             //    Key : iterator, Value : number of coefs
      }                                                         //   End loop over interaction list
    }                                                           //  End loop over target cells
    stopTimer("Get list     ");                                 //  Stop timer
    setSourceCell();                                            //  Set source buffer for cells
    setTargetCell(listM2L,flagM2L);                             //  Set target buffer for cells
    selectM2L();                                                //  Select M2L kernel
    getTargetCell(listM2L,false);                               //  Get body values from target buffer
    clearBuffers();                                             //  Clear GPU buffers
    ioffset += numCell;                                         //  Increment ioffset
  }                                                             // End loop over icall
  listM2L.clear();                                              // Clear interaction lists
  flagM2L.clear();                                              // Clear periodic image flags
}

void Evaluator::evalM2P(Cells &cells, bool kernel) {            // Evaluate M2P
  CI0 = cells.begin();                                          // Set begin iterator for target
  const int numCell = MAXCELL/NCRIT/7;                          // Number of cells per icall
  int numIcall = int(cells.size()-1)/numCell+1;                 // Number of icall loops
  int ioffset = 0 * kernel;                                     // Initialzie offset for icall loops
  for( int icall=0; icall!=numIcall; ++icall ) {                // Loop over icall
    CIB = cells.begin()+ioffset;                                //  Set begin iterator for target per call
    CIE = cells.begin()+std::min(ioffset+numCell,int(cells.size()));// Set end iterator for target per call
    constHost.push_back(2*R0);                                  //  Copy domain size to GPU buffer
    startTimer("Get list     ");                                //  Start timer
    for( CI=CIB; CI!=CIE; ++CI ) {                              //  Loop over target cells
      for( L_iter L=listM2P[CI-CI0].begin(); L!=listM2P[CI-CI0].end(); ++L ) {//  Loop over interaction list
        CJ = *L;                                                //    Set source cell
        sourceSize[CJ] = 2 * NCOEF;                             //    Key : iterator, Value : number of coefs
      }                                                         //   End loop over interaction list
    }                                                           //  End loop over target cells
    stopTimer("Get list     ");                                 //  Stop timer
    setSourceCell();                                            //  Set source buffer for cells
    setTargetBody(listM2P,flagM2P);                             //  Set target buffer for bodies
    selectM2P();                                                //  Select M2P kernel
    getTargetBody(listM2P);                                     //  Get body values from target buffer
    clearBuffers();                                             //  Clear GPU buffers
    ioffset += numCell;                                         //  Increment ioffset
  }                                                             // End loop over icall
  listM2P.clear();                                              // Clear interaction lists
  flagM2P.clear();                                              // Clear periodic image flags
}

void Evaluator::evalP2P(Cells &cells, bool kernel) {            // Evaluate P2P
  CI0 = cells.begin();                                          // Set begin iterator
  const int numCell = MAXCELL/NCRIT/7;                          // Number of cells per icall
  int numIcall = int(cells.size()-1)/numCell+1;                 // Number of icall loops
  int ioffset = 0 * kernel;                                     // Initialzie offset for icall loops
  for( int icall=0; icall!=numIcall; ++icall ) {                // Loop over icall
    CIB = cells.begin()+ioffset;                                //  Set begin iterator for target per call
    CIE = cells.begin()+std::min(ioffset+numCell,int(cells.size()));// Set end iterator for target per call
    constHost.push_back(2*R0);                                  //  Copy domain size to GPU buffer
    startTimer("Get list     ");                                //  Start timer
    for( CI=CIB; CI!=CIE; ++CI ) {                              //  Loop over target cells
      for( L_iter L=listP2P[CI-CI0].begin(); L!=listP2P[CI-CI0].end(); ++L ) {//  Loop over interaction list
        CJ = *L;                                                //    Set source cell
        sourceSize[CJ] = CJ->NLEAF;                             //    Key : iterator, Value : number of leafs
      }                                                         //   End loop over interaction list
    }                                                           //  End loop over target cells
    stopTimer("Get list     ");                                 //  Stop timer
    setSourceBody();                                            //  Set source buffer for bodies
    setTargetBody(listP2P,flagP2P);                             //  Set target buffer for bodies
    selectP2P();                                                //  Select P2P kernel
    getTargetBody(listP2P);                                     //  Get body values from target buffer
    clearBuffers();                                             //  Clear GPU buffers
    ioffset += numCell;                                         //  Increment ioffset
  }                                                             // End loop over icall
  listP2P.clear();                                              // Clear interaction lists
  flagP2P.clear();                                              // Clear periodic image flags
}

void Evaluator::evalL2L(Cells &cells) {                         // Evaluate L2L
  CI0 = cells.begin();                                          // Set begin iterator
  const int numCell = MAXCELL/NCOEF/2;                          // Number of cells per icall
  int numIcall = int(cells.size()-1)/numCell+1;                 // Number of icall loops
  int maxLevel = getLevel(CI0->ICELL);                          // Level of twig
  int level = 1;                                                // Start level from 1
  while( level != maxLevel+1 ) {                                // While level of source is not root level
    int ioffset = 0;                                            //  Initialzie offset for icall loops
    for( int icall=0; icall!=numIcall; ++icall ) {              //  Loop over icall
      CIB = cells.begin()+ioffset;                              //   Set begin iterator for target per call
      CIE = cells.begin()+std::min(ioffset+numCell,int(cells.size()));// Set end iterator for target per call
      constHost.push_back(2*R0);                                //   Copy domain size to GPU buffer
      startTimer("Get list     ");                              //   Start timer
      Lists listL2L(cells.size());                              //   Define L2L interation list vector
      Maps  flagL2L(cells.size());                              //   Define L2L periodic image flag
      for( CI=CIE-2; CI!=CIB-1; --CI ) {                        //   Loop over cells topdown (except root cell)
        if( getLevel(CI->ICELL) == level ) {                    //    If target cell is at current level
          CJ = CI0 + CI->PARENT;                                //     Set source cell iterator
          listL2L[CI-CI0].push_back(CJ);                        //     Push source cell into L2L interaction list
          flagL2L[CI-CI0][CJ] |= Icenter;                       //     Flip bit of periodic image flag
          if( sourceSize[CJ] == 0 ) {                           //     If the source cell has not been stored yet
            sourceSize[CJ] = 2 * NCOEF;                         //      Key : iterator, Value : number of coefs
          }                                                     //     Endif for current level
        }                                                       //    Endif for stored source cell
      }                                                         //   End loop over cells topdown
      stopTimer("Get list     ");                               //   Stop timer
      setSourceCell(false);                                     //   Set source buffer for cells
      setTargetCell(listL2L,flagL2L);                           //   Set target buffer for cells
      selectL2L();                                              //   Select L2L kernel
      getTargetCell(listL2L,false);                             //   Get body values from target buffer
      clearBuffers();                                           //   Clear GPU buffers
      ioffset += numCell;                                       //   Increment ioffset
    }                                                           //  End loop over icall
    level++;                                                    //  Increment level
  }                                                             // End while loop over levels
}

void Evaluator::evalL2P(Cells &cells) {                         // Evaluate L2P
  CI0 = cells.begin();                                          // Set begin iterator
  const int numCell = MAXCELL/NCRIT/7;                          // Number of cells per icall
  int numIcall = int(cells.size()-1)/numCell+1;                 // Number of icall loops
  int ioffset = 0;                                              // Initialzie offset for icall loops
  for( int icall=0; icall!=numIcall; ++icall ) {                // Loop over icall
    CIB = cells.begin()+ioffset;                                //  Set begin iterator for target per call
    CIE = cells.begin()+std::min(ioffset+numCell,int(cells.size()));// Set end iterator for target per call
    constHost.push_back(2*R0);                                  //  Copy domain size to GPU buffer
    startTimer("Get list     ");                                //  Start timer
    Lists listL2P(cells.size());                                //  Define L2P interation list vector
    Maps  flagL2P(cells.size());                                //  Define L2P periodic image flag
    for( CI=CIB; CI!=CIE; ++CI ) {                              //  Loop over cells
      if( CI->NCHILD == 0 ) {                                   //   If cell is a twig evaluate L2P kernel
        listL2P[CI-CI0].push_back(CI);                          //    Push source cell into L2P interaction list
        flagL2P[CI-CI0][CI] |= Icenter;                         //    Flip bit of periodic image flag
        sourceSize[CI] = 2 * NCOEF;                             //    Key : iterator, Value : number of coefs
      }                                                         //   Endif for twig cells
    }                                                           //  End loop over cells topdown
    stopTimer("Get list     ");                                 //  Stop timer
    setSourceCell(false);                                       //  Set source buffer for cells
    setTargetBody(listL2P,flagL2P);                             //  Set target buffer for bodies
    selectL2P();                                                //  Select L2P kernel
    getTargetBody(listL2P);                                     //  Get body values from target buffer
    clearBuffers();                                             //  Clear GPU buffers
    ioffset += numCell;                                         //  Increment ioffset
  }                                                             // End loop over icall
}
