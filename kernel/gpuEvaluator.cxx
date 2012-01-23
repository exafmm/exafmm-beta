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
void Evaluator<equation>::setSourceBody() {                     // Set source buffer for bodies
  startTimer("Set sourceB  ");                                  // Start timer
  for( MC_iter M=sourceSize.begin(); M!=sourceSize.end(); ++M ) {// Loop over source map
    C_iter Cj = M->first;                                       //  Set source cell
    sourceBegin[Cj] = sourceHost.size() / 4;                    //  Key : iterator, Value : offset of source leafs
    for( B_iter B=Cj->LEAF; B!=Cj->LEAF+Cj->NDLEAF; ++B ) {     //  Loop over leafs in source cell
      sourceHost.push_back(B->X[0]);                            //   Copy x position to GPU buffer
      sourceHost.push_back(B->X[1]);                            //   Copy y position to GPU buffer
      sourceHost.push_back(B->X[2]);                            //   Copy z position to GPU buffer
      sourceHost.push_back(B->SRC);                             //   Copy source value to GPU buffer
    }                                                           //  End loop over leafs
  }                                                             // End loop over source map
  stopTimer("Set sourceB  ");                                   // Stop timer
}

template<Equation equation>
void Evaluator<equation>::setSourceCell(bool isM) {             // Set source buffer for cells
  startTimer("Set sourceC  ");                                  // Start timer
  for( MC_iter M=sourceSize.begin(); M!=sourceSize.end(); ++M ) {// Loop over source map
    C_iter Cj = M->first;                                       //  Set source cell
    sourceBegin[Cj] = sourceHost.size();                        //  Key : iterator, Value : offset of sources
    sourceHost.push_back(Cj->X[0]);                             //  Copy x position to GPU buffer
    sourceHost.push_back(Cj->X[1]);                             //  Copy y position to GPU buffer
    sourceHost.push_back(Cj->X[2]);                             //  Copy z position to GPU buffer
    if( isM ) {                                                 //  If source is M
      for( int i=0; i!=NTERM; ++i ) {                           //   Loop over coefs in source cell
        sourceHost.push_back(std::real(Cj->M[i]));              //    Copy real multipole to GPU buffer
        sourceHost.push_back(std::imag(Cj->M[i]));              //    Copy imaginary multipole to GPU buffer
      }                                                         //   End loop over coefs
    } else {                                                    //  If source is L
      for( int i=0; i!=NTERM; ++i ) {                           //   Loop over coefs in source cell
        sourceHost.push_back(std::real(Cj->L[i]));              //    Copy real multipole to GPU buffer
        sourceHost.push_back(std::imag(Cj->L[i]));              //    Copy imaginary multipole to GPU buffer
      }                                                         //   End loop over coefs
    }                                                           //  Endif for source type
  }                                                             // End loop over source map
  stopTimer("Set sourceC  ");                                   // Stop timer
}

template<Equation equation>
void Evaluator<equation>::setTargetBody(Lists lists, Maps flags) {// Set target buffer for bodies
  startTimer("Set targetB  ");                                  // Start timer
  int key = 0;                                                  // Initialize key to range of coefs in source cells
  for( C_iter Ci=CiB; Ci!=CiE; ++Ci ) {                         // Loop over target cells
    if( !lists[Ci-Ci0].empty() ) {                              //  If the interation list is not empty
      int blocks = (Ci->NDLEAF - 1) / THREADS + 1;              //   Number of thread blocks needed for this target cell
      for( int i=0; i!=blocks; ++i ) {                          //   Loop over thread blocks
        keysHost.push_back(key);                                //    Save key to range of leafs in source cells
      }                                                         //   End loop over thread blocks
      key += 3*lists[Ci-Ci0].size()+1;                          //   Increment key counter
      rangeHost.push_back(lists[Ci-Ci0].size());                //   Save size of interaction list
      for( LC_iter L=lists[Ci-Ci0].begin(); L!=lists[Ci-Ci0].end(); ++L ) {//  Loop over interaction list
        C_iter Cj = *L;                                         //   Set source cell
        rangeHost.push_back(sourceBegin[Cj]);                   //    Set begin index of coefs in source cell
        rangeHost.push_back(sourceSize[Cj]);                    //    Set number of coefs in source cell
        rangeHost.push_back(flags[Ci-Ci0][Cj]);                 //    Set periodic image flag of source cell
      }                                                         //   End loop over interaction list
      targetBegin[Ci] = targetHost.size() / 4;                  //   Key : iterator, Value : offset of target leafs
      for( B_iter B=Ci->LEAF; B!=Ci->LEAF+Ci->NDLEAF; ++B ) {   //   Loop over leafs in target cell
        targetHost.push_back(B->X[0]);                          //    Copy x position to GPU buffer
        targetHost.push_back(B->X[1]);                          //    Copy y position to GPU buffer
        targetHost.push_back(B->X[2]);                          //    Copy z position to GPU buffer
        targetHost.push_back(B->SRC);                           //    Copy target value to GPU buffer
      }                                                         //   End loop over leafs
      int numPad = blocks * THREADS - Ci->NDLEAF;               //   Number of elements to pad in target GPU buffer
      for( int i=0; i!=numPad; ++i ) {                          //   Loop over elements to pad
        targetHost.push_back(0);                                //    Pad x position in GPU buffer
        targetHost.push_back(0);                                //    Pad y position in GPU buffer
        targetHost.push_back(0);                                //    Pad z position in GPU buffer
        targetHost.push_back(0);                                //    Pad target value to GPU buffer
      }                                                         //   End loop over elements to pad
    }                                                           //  End if for empty interation list
  }                                                             // End loop over target cells
  stopTimer("Set targetB  ");                                   // Stop timer
}

template<Equation equation>
void Evaluator<equation>::setTargetCell(Lists lists, Maps flags) {// Set target buffer for cells
  startTimer("Set targetC  ");                                  // Start timer
  int key = 0;                                                  // Initialize key to range of coefs in target cells
  for( C_iter Ci=CiB; Ci!=CiE; ++Ci ) {                         // Loop over target cells
    if( !lists[Ci-Ci0].empty() ) {                              //  If the interation list is not empty
      keysHost.push_back(key);                                  //   Save key to range of coefs in target cells
      key += 3*lists[Ci-Ci0].size()+1;                          //   Increment key counter
      rangeHost.push_back(lists[Ci-Ci0].size());                //   Save size of interaction list
      for( LC_iter L=lists[Ci-Ci0].begin(); L!=lists[Ci-Ci0].end(); ++L ) {//  Loop over interaction list
        C_iter Cj = *L;                                         //    Set target cell
        int begin = sourceBegin[Cj];                            //    Get begin index of coefs in source cell
        int size = sourceSize[Cj];                              //    Get number of coefs in source cell
        rangeHost.push_back(begin);                             //    Set begin index of coefs in source cell
        rangeHost.push_back(size);                              //    Set number of coefs in source cell
        rangeHost.push_back(flags[Ci-Ci0][Cj]);                 //    Set periodic image flag of source cell
      }                                                         //   End loop over interaction list
      targetBegin[Ci] = targetHost.size();                      //   Key : iterator, Value : offset of target coefs
      targetHost.push_back(Ci->X[0]);                           //   Copy x position to GPU buffer
      targetHost.push_back(Ci->X[1]);                           //   Copy y position to GPU buffer
      targetHost.push_back(Ci->X[2]);                           //   Copy z position to GPU buffer
      for( int i=0; i!=2*NTERM; ++i ) {                         //   Loop over coefs in target cell
        targetHost.push_back(0);                                //    Pad GPU buffer
      }                                                         //   End loop over coefs
      int numPad = 2 * THREADS - 2 * NTERM - 3;                 //   Number of elements to pad in target GPU buffer
      assert(numPad >= 0);                                      //   THREADS must be large enough
      for( int i=0; i!=numPad; ++i ) {                          //   Loop over elements to pad
        targetHost.push_back(0);                                //    Pad GPU buffer
      }                                                         //   End loop over elements to pad
    }                                                           //  End if for empty interation list
  }                                                             // End loop over target cells
  stopTimer("Set targetC  ");                                   // Stop timer
}

template<Equation equation>
void Evaluator<equation>::getTargetBody(Lists &lists) {         // Get body values from target buffer
  startTimer("Get targetB  ");                                  // Start timer
  for( C_iter Ci=CiB; Ci!=CiE; ++Ci ) {                         // Loop over target cells
    if( !lists[Ci-Ci0].empty() ) {                              //  If the interation list is not empty
      int begin = targetBegin[Ci];                              //   Offset of target leafs
        for( B_iter B=Ci->LEAF; B!=Ci->LEAF+Ci->NDLEAF; ++B ) { //    Loop over target bodies
          B->TRG[0] += targetHost[4*(begin+B-Ci->LEAF)+0];      //     Copy 1st target value from GPU buffer
          B->TRG[1] += targetHost[4*(begin+B-Ci->LEAF)+1];      //     Copy 2nd target value from GPU buffer
          B->TRG[2] += targetHost[4*(begin+B-Ci->LEAF)+2];      //     Copy 3rd target value from GPU buffer
          B->TRG[3] += targetHost[4*(begin+B-Ci->LEAF)+3];      //     Copy 4th target value from GPU buffer
        }                                                       //    End loop over target bodies
      lists[Ci-Ci0].clear();                                    //   Clear interaction list
    }                                                           //  End if for empty interation list
  }                                                             // End loop over target cells
  stopTimer("Get targetB  ");                                   // Stop timer
}

template<Equation equation>
void Evaluator<equation>::getTargetCell(Lists &lists, bool isM) {// Get body values from target buffer
  startTimer("Get targetC  ");                                  // Start timer
  for( C_iter Ci=CiB; Ci!=CiE; ++Ci ) {                         // Loop over target cells
    if( !lists[Ci-Ci0].empty() ) {                              //  If the interation list is not empty
      int begin = targetBegin[Ci];                              //   Offset of target coefs
      if( isM ) {                                               //   If target is M
        for( int i=0; i!=NTERM; ++i ) {                         //    Loop over coefs in target cell
          Ci->M[i].real() += targetHost[begin+2*i+0];           //     Copy real target values from GPU buffer
          Ci->M[i].imag() += targetHost[begin+2*i+1];           //     Copy imaginary target values from GPU buffer
        }                                                       //    End loop over coefs
      } else {                                                  //   If target is L
        for( int i=0; i!=NTERM; ++i ) {                         //    Loop over coefs in target cell
          Ci->L[i].real() += targetHost[begin+2*i+0];           //     Copy real target values from GPU buffer
          Ci->L[i].imag() += targetHost[begin+2*i+1];           //     Copy imaginary target values from GPU buffer
        }                                                       //    End loop over coefs
      }                                                         //   Endif for target type
      lists[Ci-Ci0].clear();                                    //   Clear interaction list
    }                                                           //  End if for empty interation list
  }                                                             // End loop over target cells
  stopTimer("Get targetC  ");                                   // Stop timer
}

template<Equation equation>
void Evaluator<equation>::clearBuffers() {                      // Clear GPU buffers
  startTimer("Clear buffer ");                                  // Start timer
  constHost.clear();                                            // Clear const vector
  keysHost.clear();                                             // Clear keys vector
  rangeHost.clear();                                            // Clear range vector
  sourceHost.clear();                                           // Clear source vector
  targetHost.clear();                                           // Clear target vector
  sourceBegin.clear();                                          // Clear map for offset of source cells
  sourceSize.clear();                                           // Clear map for size of source cells
  targetBegin.clear();                                          // Clear map for offset of target cells
  stopTimer("Clear buffer ");                                   // Stop timer
}

template<Equation equation>
void Evaluator<equation>::evalP2P(Bodies &ibodies, Bodies &jbodies, bool onCPU) {// Evaluate all P2P kernels
  int numIcall = int(ibodies.size()-1)/MAXBODY+1;               // Number of icall loops
  int numJcall = int(jbodies.size()-1)/MAXBODY+1;               // Number of jcall loops
  int ioffset = 0;                                              // Initialzie offset for icall loops
  Cells cells(2);                                               // Two cells to put target and source bodies
  for( int icall=0; icall!=numIcall; ++icall ) {                // Loop over icall
    B_iter Bi0 = ibodies.begin()+ioffset;                       //  Set target bodies begin iterator
    B_iter BiN = ibodies.begin()+std::min(ioffset+MAXBODY,int(ibodies.size()));// Set target bodies end iterator
    cells[0].LEAF = Bi0;                                        //  Iterator of first target leaf
    cells[0].NDLEAF = BiN-Bi0;                                  //  Number of target leafs
    int joffset = 0;                                            //  Initialize offset for jcall loops
    for( int jcall=0; jcall!=numJcall; ++jcall ) {              //  Loop over jcall
      B_iter Bj0 = jbodies.begin()+joffset;                     //  Set source bodies begin iterator
      B_iter BjN = jbodies.begin()+std::min(joffset+MAXBODY,int(jbodies.size()));// Set source bodies end iterator
      cells[1].LEAF = Bj0;                                      //  Iterator of first source leaf
      cells[1].NDLEAF = BjN-Bj0;                                //  Number of source leafs
      C_iter Ci = cells.begin(), Cj = cells.begin()+1;          //  Iterator of target and source cells
      if( onCPU ) {                                             //  If calculation is to be done on CPU
        Xperiodic = 0;                                          //   Set periodic coordinate offset
        P2P(Ci,Cj);                                             //   Perform P2P kernel on CPU
      } else {                                                  //  If calculation is to be done on GPU
        constHost.push_back(2*R0);                              //   Copy domain size to GPU buffer
        for( B_iter B=Bj0; B!=BjN; ++B ) {                      //   Loop over source bodies
          sourceHost.push_back(B->X[0]);                        //   Copy x position to GPU buffer
          sourceHost.push_back(B->X[1]);                        //   Copy y position to GPU buffer
          sourceHost.push_back(B->X[2]);                        //   Copy z position to GPU buffer
          sourceHost.push_back(B->SRC);                         //   Copy source value to GPU buffer
        }                                                       //   End loop over source bodies
        int key = 0;                                            //   Initialize key to range of leafs in source cells
        int blocks = (BiN - Bi0 - 1) / THREADS + 1;             //   Number of thread blocks needed for this target cell
        for( int i=0; i!=blocks; ++i ) {                        //   Loop over thread blocks
          keysHost.push_back(key);                              //    Save key to range of leafs in source cells
        }                                                       //   End loop over thread blocks
        rangeHost.push_back(1);                                 //   Save size of interaction list
        rangeHost.push_back(0);                                 //   Set begin index of leafs
        rangeHost.push_back(BjN-Bj0);                           //   Set number of leafs
        rangeHost.push_back(Icenter);                           //   Set periodic image flag
        for( B_iter B=Bi0; B!=BiN; ++B ) {                      //   Loop over target bodies
          targetHost.push_back(B->X[0]);                        //    Copy x position to GPU buffer
          targetHost.push_back(B->X[1]);                        //    Copy y position to GPU buffer
          targetHost.push_back(B->X[2]);                        //    Copy z position to GPU buffer
          targetHost.push_back(B->SRC);                         //    Copy target value to GPU buffer
        }                                                       //   End loop over target bodies
        int numPad = blocks * THREADS - (BiN - Bi0);            //   Number of elements to pad in target GPU buffer
        for( int i=0; i!=numPad; ++i ) {                        //   Loop over elements to pad
          targetHost.push_back(0);                              //    Pad x position in GPU buffer
          targetHost.push_back(0);                              //    Pad y position in GPU buffer
          targetHost.push_back(0);                              //    Pad z position in GPU buffer
          targetHost.push_back(0);                              //    Pad target value to GPU buffer
        }                                                       //   End loop over elements to pad
        allocate();                                             //   Allocate GPU memory
        hostToDevice();                                         //   Copy from host to device
        P2P();                                                  //   Perform P2P kernel
        deviceToHost();                                         //   Copy from device to host
        for( B_iter B=Bi0; B!=BiN; ++B ) {                      //   Loop over target bodies
          B->TRG[0] += targetHost[4*(B-Bi0)+0];                 //    Copy 1st target value from GPU buffer
          B->TRG[1] += targetHost[4*(B-Bi0)+1];                 //    Copy 2nd target value from GPU buffer
          B->TRG[2] += targetHost[4*(B-Bi0)+2];                 //    Copy 3rd target value from GPU buffer
          B->TRG[3] += targetHost[4*(B-Bi0)+3];                 //    Copy 4th target value from GPU buffer
        }                                                       //   End loop over target bodies
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

template<Equation equation>
void Evaluator<equation>::evalP2M(Cells &cells) {               // Evaluate all P2M kernels
  Ci0 = cells.begin();                                          // Set begin iterator for target
  const int numCell = MAXCELL/NCRIT/4;                          // Number of cells per icall
  int numIcall = int(cells.size()-1)/numCell+1;                 // Number of icall loops
  int ioffset = 0;                                              // Initialzie offset for icall loops
  for( int icall=0; icall!=numIcall; ++icall ) {                // Loop over icall
    CiB = cells.begin()+ioffset;                                //  Set begin iterator for target per call
    CiE = cells.begin()+std::min(ioffset+numCell,int(cells.size()));// Set end iterator for target per call
    constHost.push_back(2*R0);                                  //  Copy domain size to GPU buffer
    startTimer("Get list     ");                                //  Start timer
    Lists listP2M(cells.size());                                //  Define P2M interation list vector
    Maps  flagP2M(cells.size());                                //  Define P2M periodic image flag
    for( C_iter Ci=CiB; Ci!=CiE; ++Ci ) {                       //  Loop over target cells
      Ci->M = Ci->L = 0;                                        //   Initialize multipole & local coefficients
      if( Ci->NCHILD == 0 ) {                                   //   If cell is a twig
        listP2M[Ci-Ci0].push_back(Ci);                          //    Push source cell into P2M interaction list
        flagP2M[Ci-Ci0][Ci] |= Icenter;                         //    Flip bit of periodic image flag
        sourceSize[Ci] = Ci->NDLEAF;                            //    Key : iterator, Value : number of leafs
      }                                                         //   End loop over cells topdown
    }                                                           //  End loop over source map
    stopTimer("Get list     ");                                 //  Stop timer
    setSourceBody();                                            //  Set source buffer for bodies
    setTargetCell(listP2M,flagP2M);                             //  Set target buffer for cells
    allocate();                                                 //  Allocate GPU memory
    hostToDevice();                                             //  Copy from host to device
    P2M();                                                      //  Perform P2M kernel
    deviceToHost();                                             //  Copy from device to host
    getTargetCell(listP2M,true);                                //  Get body values from target buffer
    clearBuffers();                                             //  Clear GPU buffers
    ioffset += numCell;                                         //  Increment ioffset
  }                                                             // End loop over icall
}

template<Equation equation>
void Evaluator<equation>::evalM2M(Cells &cells, Cells &jcells) {// Evaluate all M2M kernels
  Ci0 = cells.begin();                                          // Set begin iterator for target
  Cj0 = jcells.begin();                                         // Set begin iterator for source
  const int numCell = MAXCELL / NTERM / 2;                      // Number of cells per icall
  int numIcall = int(cells.size()-1) / numCell + 1;             // Number of icall loops
  int level = getLevel(cells.front().ICELL);                    // Level of twig
  while( level != -1 ) {                                        // While level of target is not past root level
    int ioffset = 0;                                            //  Initialzie offset for icall loops
    for( int icall=0; icall!=numIcall; ++icall ) {              //  Loop over icall
      CiB = cells.begin()+ioffset;                              //   Set begin iterator for target per call
      CiE = cells.begin()+std::min(ioffset+numCell,int(cells.size()));// Set end iterator for target per call
      constHost.push_back(2*R0);                                //   Copy domain size to GPU buffer
      startTimer("Get list     ");                              //   Start timer
      Lists listM2M(cells.size());                              //   Define M2M interation list vector
      Maps  flagM2M(cells.size());                              //   Define M2M periodic image flag
      for( C_iter Ci=CiB; Ci!=CiE; ++Ci ) {                     //   Loop over cells bottomup (except root cell)
        if( getLevel(Ci->ICELL) == level ) {                    //    If target cell is at current level
          for( int i=0; i<Ci->NCHILD; ++i ) {                   //     Loop over child cells
            C_iter Cj = Cj0 + Ci->CHILD+i;                      //      Set iterator for source cell
            listM2M[Ci-Ci0].push_back(Cj);                      //      Push source cell into M2M interaction list
            flagM2M[Ci-Ci0][Cj] |= Icenter;                     //      Flip bit of periodic image flag
            sourceSize[Cj] = 2 * NTERM;                         //      Key : iterator, Value : number of coefs
          }                                                     //     End loop over child cells
        }                                                       //    Endif for current level
      }                                                         //   End loop over cells
      stopTimer("Get list     ");                               //   Stop timer
      setSourceCell(true);                                      //   Set source buffer for cells
      setTargetCell(listM2M,flagM2M);                           //   Set target buffer for cells
      allocate();                                               //   Allocate GPU memory
      hostToDevice();                                           //   Copy from host to device
      M2M();                                                    //   Perform M2M kernel
      deviceToHost();                                           //   Copy from device to host
      getTargetCell(listM2M,true);                              //   Get body values from target buffer
      clearBuffers();                                           //   Clear GPU buffers
      ioffset += numCell;                                       //   Increment ioffset
    }                                                           //  End loop over icall
    level--;                                                    //  Decrement level
  }                                                             // End while loop over levels
}

template<Equation equation>
void Evaluator<equation>::evalM2L(C_iter Ci, C_iter Cj) {       // Queue single M2L kernel
  listM2L[Ci-Ci0].push_back(Cj);                                // Push source cell into M2L interaction list
  flagM2L[Ci-Ci0][Cj] |= Iperiodic;                             // Flip bit of periodic image flag
  NM2L++;                                                       // Count M2L kernel execution
}

template<Equation equation>
void Evaluator<equation>::evalM2L(Cells &cells) {               // Evaluate queued M2L kernels
  Ci0 = cells.begin();                                          // Set begin iterator
  const int numCell = MAXCELL / NTERM / 2;                      // Number of cells per icall
  int numIcall = int(cells.size()-1)/numCell+1;                 // Number of icall loops
  int ioffset = 0;                                              // Initialzie offset for icall loops
  for( int icall=0; icall!=numIcall; ++icall ) {                // Loop over icall
    CiB = cells.begin()+ioffset;                                //  Set begin iterator for target per call
    CiE = cells.begin()+std::min(ioffset+numCell,int(cells.size()));// Set end iterator for target per call
    constHost.push_back(2*R0);                                  //  Copy domain size to GPU buffer
    startTimer("Get list     ");                                //  Start timer
    for( C_iter Ci=CiB; Ci!=CiE; ++Ci ) {                       //  Loop over target cells
      for( LC_iter L=listM2L[Ci-Ci0].begin(); L!=listM2L[Ci-Ci0].end(); ++L ) {//  Loop over interaction list
        C_iter Cj = *L;                                         //    Set source cell
        sourceSize[Cj] = 2 * NTERM;                             //    Key : iterator, Value : number of coefs
      }                                                         //   End loop over interaction list
    }                                                           //  End loop over target cells
    stopTimer("Get list     ");                                 //  Stop timer
    setSourceCell(true);                                        //  Set source buffer for cells
    setTargetCell(listM2L,flagM2L);                             //  Set target buffer for cells
    allocate();                                                 //  Allocate GPU memory
    hostToDevice();                                             //  Copy from host to device
    M2L();                                                      //  Perform M2L kernel
    deviceToHost();                                             //  Copy from device to host
    getTargetCell(listM2L,false);                               //  Get body values from target buffer
    clearBuffers();                                             //  Clear GPU buffers
    ioffset += numCell;                                         //  Increment ioffset
  }                                                             // End loop over icall
  listM2L.clear();                                              // Clear interaction lists
  flagM2L.clear();                                              // Clear periodic image flags
}

template<Equation equation>
void Evaluator<equation>::evalM2P(C_iter Ci, C_iter Cj) {       // Queue single M2P kernel
  listM2P[Ci-Ci0].push_back(Cj);                                // Push source cell into M2P interaction list
  flagM2P[Ci-Ci0][Cj] |= Iperiodic;                             // Flip bit of periodic image flag
  NM2P++;                                                       // Count M2P kernel execution
}

template<Equation equation>
void Evaluator<equation>::evalM2P(Cells &cells) {               // Evaluate queued M2P kernels
  Ci0 = cells.begin();                                          // Set begin iterator for target
  const int numCell = MAXCELL/NCRIT/4;                          // Number of cells per icall
  int numIcall = int(cells.size()-1)/numCell+1;                 // Number of icall loops
  int ioffset = 0;                                              // Initialzie offset for icall loops
  for( int icall=0; icall!=numIcall; ++icall ) {                // Loop over icall
    CiB = cells.begin()+ioffset;                                //  Set begin iterator for target per call
    CiE = cells.begin()+std::min(ioffset+numCell,int(cells.size()));// Set end iterator for target per call
    constHost.push_back(2*R0);                                  //  Copy domain size to GPU buffer
    startTimer("Get list     ");                                //  Start timer
    for( C_iter Ci=CiB; Ci!=CiE; ++Ci ) {                       //  Loop over target cells
      for( LC_iter L=listM2P[Ci-Ci0].begin(); L!=listM2P[Ci-Ci0].end(); ++L ) {//  Loop over interaction list
        C_iter Cj = *L;                                         //    Set source cell
        sourceSize[Cj] = 2 * NTERM;                             //    Key : iterator, Value : number of coefs
      }                                                         //   End loop over interaction list
    }                                                           //  End loop over target cells
    stopTimer("Get list     ");                                 //  Stop timer
    setSourceCell(true);                                        //  Set source buffer for cells
    setTargetBody(listM2P,flagM2P);                             //  Set target buffer for bodies
    allocate();                                                 //  Allocate GPU memory
    hostToDevice();                                             //  Copy from host to device
    M2P();                                                      //  Perform M2P kernel
    deviceToHost();                                             //  Copy from device to host
    getTargetBody(listM2P);                                     //  Get body values from target buffer
    clearBuffers();                                             //  Clear GPU buffers
    ioffset += numCell;                                         //  Increment ioffset
  }                                                             // End loop over icall
  listM2P.clear();                                              // Clear interaction lists
  flagM2P.clear();                                              // Clear periodic image flags
}

template<Equation equation>
void Evaluator<equation>::evalP2P(C_iter Ci, C_iter Cj) {       // Queue single P2P kernel
  listP2P[Ci-Ci0].push_back(Cj);                                // Push source cell into P2P interaction list
  flagP2P[Ci-Ci0][Cj] |= Iperiodic;                             // Flip bit of periodic image flag
  NP2P++;                                                       // Count P2P kernel execution
}

template<Equation equation>
void Evaluator<equation>::evalP2P(Cells &cells) {               // Evaluate queued P2P kernels
  Ci0 = cells.begin();                                          // Set begin iterator
  const int numCell = MAXCELL/NCRIT/4;                          // Number of cells per icall
  int numIcall = int(cells.size()-1)/numCell+1;                 // Number of icall loops
  int ioffset = 0;                                              // Initialzie offset for icall loops
  for( int icall=0; icall!=numIcall; ++icall ) {                // Loop over icall
    CiB = cells.begin()+ioffset;                                //  Set begin iterator for target per call
    CiE = cells.begin()+std::min(ioffset+numCell,int(cells.size()));// Set end iterator for target per call
    constHost.push_back(2*R0);                                  //  Copy domain size to GPU buffer
    startTimer("Get list     ");                                //  Start timer
    for( C_iter Ci=CiB; Ci!=CiE; ++Ci ) {                       //  Loop over target cells
      for( LC_iter L=listP2P[Ci-Ci0].begin(); L!=listP2P[Ci-Ci0].end(); ++L ) {//  Loop over interaction list
        C_iter Cj = *L;                                         //    Set source cell
        sourceSize[Cj] = Cj->NDLEAF;                            //    Key : iterator, Value : number of leafs
      }                                                         //   End loop over interaction list
    }                                                           //  End loop over target cells
    stopTimer("Get list     ");                                 //  Stop timer
    setSourceBody();                                            //  Set source buffer for bodies
    setTargetBody(listP2P,flagP2P);                             //  Set target buffer for bodies
    allocate();                                                 //  Allocate GPU memory
    hostToDevice();                                             //  Copy from host to device
    P2P();                                                      //  Perform P2P kernel
    deviceToHost();                                             //  Copy from device to host
    getTargetBody(listP2P);                                     //  Get body values from target buffer
    clearBuffers();                                             //  Clear GPU buffers
    ioffset += numCell;                                         //  Increment ioffset
  }                                                             // End loop over icall
  listP2P.clear();                                              // Clear interaction lists
  flagP2P.clear();                                              // Clear periodic image flags
}

template<Equation equation>
void Evaluator<equation>::evalL2L(Cells &cells) {               // Evaluate all L2L kenrels
  Ci0 = cells.begin();                                          // Set begin iterator
  const int numCell = MAXCELL / NTERM / 2;                      // Number of cells per icall
  int numIcall = int(cells.size()-1)/numCell+1;                 // Number of icall loops
  int maxLevel = getLevel(Ci0->ICELL);                          // Level of twig
  int level = 1;                                                // Start level from 1
  while( level != maxLevel+1 ) {                                // While level of source is not root level
    int ioffset = 0;                                            //  Initialzie offset for icall loops
    for( int icall=0; icall!=numIcall; ++icall ) {              //  Loop over icall
      CiB = cells.begin()+ioffset;                              //   Set begin iterator for target per call
      CiE = cells.begin()+std::min(ioffset+numCell,int(cells.size()));// Set end iterator for target per call
      constHost.push_back(2*R0);                                //   Copy domain size to GPU buffer
      startTimer("Get list     ");                              //   Start timer
      Lists listL2L(cells.size());                              //   Define L2L interation list vector
      Maps  flagL2L(cells.size());                              //   Define L2L periodic image flag
      for( C_iter Ci=CiE-2; Ci!=CiB-1; --Ci ) {                 //   Loop over cells topdown (except root cell)
        if( getLevel(Ci->ICELL) == level ) {                    //    If target cell is at current level
          C_iter Cj = Ci0 + Ci->PARENT;                         //     Set source cell iterator
          listL2L[Ci-Ci0].push_back(Cj);                        //     Push source cell into L2L interaction list
          flagL2L[Ci-Ci0][Cj] |= Icenter;                       //     Flip bit of periodic image flag
          if( sourceSize[Cj] == 0 ) {                           //     If the source cell has not been stored yet
            sourceSize[Cj] = 2 * NTERM;                         //      Key : iterator, Value : number of coefs
          }                                                     //     Endif for current level
        }                                                       //    Endif for stored source cell
      }                                                         //   End loop over cells topdown
      stopTimer("Get list     ");                               //   Stop timer
      setSourceCell(false);                                     //   Set source buffer for cells
      setTargetCell(listL2L,flagL2L);                           //   Set target buffer for cells
      allocate();                                               //   Allocate GPU memory
      hostToDevice();                                           //   Copy from host to device
      L2L();                                                    //   Perform L2L kernel
      deviceToHost();                                           //   Copy from device to host
      getTargetCell(listL2L,false);                             //   Get body values from target buffer
      clearBuffers();                                           //   Clear GPU buffers
      ioffset += numCell;                                       //   Increment ioffset
    }                                                           //  End loop over icall
    level++;                                                    //  Increment level
  }                                                             // End while loop over levels
}

template<Equation equation>
void Evaluator<equation>::evalL2P(Cells &cells) {               // Evaluate all L2P kernels
  Ci0 = cells.begin();                                          // Set begin iterator
  const int numCell = MAXCELL/NCRIT/4;                          // Number of cells per icall
  int numIcall = int(cells.size()-1)/numCell+1;                 // Number of icall loops
  int ioffset = 0;                                              // Initialzie offset for icall loops
  for( int icall=0; icall!=numIcall; ++icall ) {                // Loop over icall
    CiB = cells.begin()+ioffset;                                //  Set begin iterator for target per call
    CiE = cells.begin()+std::min(ioffset+numCell,int(cells.size()));// Set end iterator for target per call
    constHost.push_back(2*R0);                                  //  Copy domain size to GPU buffer
    startTimer("Get list     ");                                //  Start timer
    Lists listL2P(cells.size());                                //  Define L2P interation list vector
    Maps  flagL2P(cells.size());                                //  Define L2P periodic image flag
    for( C_iter Ci=CiB; Ci!=CiE; ++Ci ) {                       //  Loop over cells
      if( Ci->NCHILD == 0 ) {                                   //   If cell is a twig evaluate L2P kernel
        listL2P[Ci-Ci0].push_back(Ci);                          //    Push source cell into L2P interaction list
        flagL2P[Ci-Ci0][Ci] |= Icenter;                         //    Flip bit of periodic image flag
        sourceSize[Ci] = 2 * NTERM;                             //    Key : iterator, Value : number of coefs
      }                                                         //   Endif for twig cells
    }                                                           //  End loop over cells topdown
    stopTimer("Get list     ");                                 //  Stop timer
    setSourceCell(false);                                       //  Set source buffer for cells
    setTargetBody(listL2P,flagL2P);                             //  Set target buffer for bodies
    allocate();                                                 //  Allocate GPU memory
    hostToDevice();                                             //  Copy from host to device
    L2P();                                                      //  Perform L2P kernel
    deviceToHost();                                             //  Copy from device to host
    getTargetBody(listL2P);                                     //  Get body values from target buffer
    clearBuffers();                                             //  Clear GPU buffers
    ioffset += numCell;                                         //  Increment ioffset
  }                                                             // End loop over icall
}

template<Equation equation>
void Evaluator<equation>::evalEwaldReal(C_iter Ci, C_iter Cj) { // Queue single Ewald real kernel
  listP2P[Ci-Ci0].push_back(Cj);                                // Push source cell into P2P interaction list
  flagP2P[Ci-Ci0][Cj] |= Iperiodic;                             // Flip bit of periodic image flag
}

template<Equation equation>
void Evaluator<equation>::evalEwaldReal(Cells &cells) {         // Evaluate queued Ewald real kernels
  Ci0 = cells.begin();                                          // Set begin iterator
  const int numCell = MAXCELL/NCRIT/4;                          // Number of cells per icall
  int numIcall = int(cells.size()-1)/numCell+1;                 // Number of icall loops
  int ioffset = 0;                                              // Initialzie offset for icall loops
  for( int icall=0; icall!=numIcall; ++icall ) {                // Loop over icall
    CiB = cells.begin()+ioffset;                                //  Set begin iterator for target per call
    CiE = cells.begin()+std::min(ioffset+numCell,int(cells.size()));// Set end iterator for target per call
    constHost.push_back(2*R0);                                  //  Copy domain size to GPU buffer
    constHost.push_back(ALPHA);                                 //  Copy Ewald scaling to GPU buffer
    startTimer("Get list     ");                                //  Start timer
    for( C_iter Ci=CiB; Ci!=CiE; ++Ci ) {                       //  Loop over target cells
      for( LC_iter L=listP2P[Ci-Ci0].begin(); L!=listP2P[Ci-Ci0].end(); ++L ) {//  Loop over interaction list
        C_iter Cj = *L;                                         //    Set source cell
        sourceSize[Cj] = Cj->NDLEAF;                            //    Key : iterator, Value : number of leafs
      }                                                         //   End loop over interaction list
    }                                                           //  End loop over target cells
    stopTimer("Get list     ");                                 //  Stop timer
    setSourceBody();                                            //  Set source buffer for bodies
    setTargetBody(listP2P,flagP2P);                             //  Set target buffer for bodies
    allocate();                                                 //  Allocate GPU memory
    hostToDevice();                                             //  Copy from host to device
    EwaldReal();                                                //  Perform Ewald real kernel
    deviceToHost();                                             //  Copy from device to host
    getTargetBody(listP2P);                                     //  Get body values from target buffer
    clearBuffers();                                             //  Clear GPU buffers
    ioffset += numCell;                                         //  Increment ioffset
  }                                                             // End loop over icall
  listP2P.clear();                                              // Clear interaction lists
  flagP2P.clear();                                              // Clear periodic image flags
}

template<Equation equation>
void Evaluator<equation>::timeKernels() {                       // Time all kernels for auto-tuning
  Bodies ibodies(100), jbodies(100);                            // Artificial bodies
  for( B_iter Bi=ibodies.begin(),Bj=jbodies.begin(); Bi!=ibodies.end(); ++Bi, ++Bj ) {// Loop over artificial bodies
    Bi->X = 0;                                                  //  Set coordinates of target body
    Bj->X = 1;                                                  //  Set coordinates of source body
  }                                                             // End loop over artificial bodies
  Cells icells, jcells;                                         // Artificial cells
  icells.resize(10);                                            // 100 artificial target cells
  jcells.resize(100);                                           // 100 artificial source cells
  Ci0 = icells.begin();                                         // Set global begin iterator for source
  for( C_iter Ci=icells.begin(); Ci!=icells.end(); ++Ci ) {     // Loop over target cells
    Ci->X = 0;                                                  //  Set coordinates of target cell
    Ci->NDLEAF = 100;                                           //  Number of leafs in target cell
    Ci->LEAF = ibodies.begin();                                 //  Leaf iterator in target cell
  }                                                             // End loop over target cells
  for( C_iter Cj=jcells.begin(); Cj!=jcells.end(); ++Cj ) {     // Loop over source cells
    Cj->X = 1;                                                  //  Set coordinates of source cell
    Cj->NDLEAF = 100;                                           //  Number of leafs in source cell
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
      listP2P[Ci-Ci0].push_back(Cj);                            //   Push source cell into P2P interaction list
      listM2P[Ci-Ci0].push_back(Cj);                            //   Push source cell into P2P interaction list
      listM2L[Ci-Ci0].push_back(Cj);                            //   Push source cell into P2P interaction list
    }                                                           //  End loop over source cells
  }                                                             // End loop over target cells
  startTimer("P2P kernel   ");                                  // Start timer
  evalP2P(icells);                                              // Evaluate queued P2P kernels
  timeP2P = stopTimer("P2P kernel   ") / 100 / 100;             // Stop timer
  startTimer("M2L kernel   ");                                  // Start timer
  evalM2L(icells);                                              // Evaluate queued M2L kernels
  timeM2L = stopTimer("M2L kernel   ");                         // Stop timer
  startTimer("M2P kernel   ");                                  // Start timer
  evalM2P(icells);                                              // Evaluate queued M2P kernels
  timeM2P = stopTimer("M2P kernel   ") / 100;                   // Stop timer
}

#if QUARK
template<Equation equation>
void Evaluator<equation>::interact(C_iter CI, C_iter CJ, Quark*) {
  PairQueue privateQueue;
  Pair pair(CI,CJ);
  privateQueue.push(pair);
  while( !privateQueue.empty() ) {
    Pair Cij = privateQueue.front();
    privateQueue.pop();
    if(splitFirst(Cij.first,Cij.second)) {
      C_iter C = Cij.first;
      for( C_iter Ci=Ci0+C->CHILD; Ci!=Ci0+C->CHILD+C->NCHILD; ++Ci ) {
        interact(Ci,Cij.second,privateQueue);
      }
    } else {
      C_iter C = Cij.second;
      for( C_iter Cj=Cj0+C->CHILD; Cj!=Cj0+C->CHILD+C->NCHILD; ++Cj ) {
        interact(Cij.first,Cj,privateQueue);
      }
    }
  }
}
#endif
