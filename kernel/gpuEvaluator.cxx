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
      if( kernelName == "Gaussian" ) {                          //  If Gaussian kernel
        for( B_iter B=BI0; B!=BIN; ++B ) {                      //   Loop over target bodies
          B->TRG[0] += targetHost[6*(begin+B-BI0)+0];           //    Copy 1st target value from GPU buffer
        }                                                       //   End loop over target bodies
      } else {                                                  //  If not Gaussian kernel
        for( B_iter B=BI0; B!=BIN; ++B ) {                      //   Loop over target bodies
          B->TRG[0] += targetHost[6*(begin+B-BI0)+0];           //    Copy 1st target value from GPU buffer
          B->TRG[1] += targetHost[6*(begin+B-BI0)+1];           //    Copy 2nd target value from GPU buffer
          B->TRG[2] += targetHost[6*(begin+B-BI0)+2];           //    Copy 3rd target value from GPU buffer
          B->TRG[3] += targetHost[6*(begin+B-BI0)+3];           //    Copy 4th target value from GPU buffer
        }                                                       //   End loop over target bodies
      }                                                         //  Endif for Gaussian kernel
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
        if( kernelName == "Laplace" ) {                         //   If Laplace kernel
          LaplaceP2P_CPU();                                     //    Evaluate P2P kernel
        } else if ( kernelName == "BiotSavart" ) {              //   If Biot Savart kernel
          BiotSavartP2P_CPU();                                  //    Evaluate P2P kernel
        } else if ( kernelName == "Stretching" ) {              //   If Stretching kernel
          StretchingP2P_CPU();                                  //    Evaluate P2P kernel
        } else if ( kernelName == "Gaussian" ) {                //   If Gaussian kernel
          GaussianP2P_CPU();                                    //    Evaluate P2P kernel
        } else {                                                //   If kernel is none of the above
          if(MPIRANK == 0) std::cout << "Invalid kernel type" << std::endl;// Invalid kernel type
          abort();                                              //   Abort execution
        }                                                       //   Endif for kernel type
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
        if( kernelName == "Laplace" ) {                         //   If Laplace kernel
          LaplaceP2P();                                         //    Evaluate P2P kernel
        } else if ( kernelName == "BiotSavart" ) {              //   If Biot Savart kernel
          BiotSavartP2P();                                      //    Evaluate P2P kernel
        } else if ( kernelName == "Stretching" ) {              //   If Stretching kernel
          StretchingP2P();                                      //    Evaluate P2P kernel
        } else if ( kernelName == "Gaussian" ) {                //   If Gaussian kernel
          GaussianP2P();                                        //    Evaluate P2P kernel
        } else {                                                //   If kernel is none of the above
          if(MPIRANK == 0) std::cout << "Invalid kernel type" << std::endl;// Invalid kernel type
          abort();                                              //   Abort execution
        }                                                       //   Endif for kernel type
        if( kernelName == "Gaussian" ) {                        //   If Gaussian kernel
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
    if( kernelName == "Laplace" ) {                             //  If Laplace kernel
      LaplaceP2M();                                             //   Evaluate P2M kernel
    } else if ( kernelName == "BiotSavart" ) {                  //  If Biot Savart kernel
      BiotSavartP2M();                                          //   Evaluate P2M kernel
    } else if ( kernelName == "Stretching" ) {                  //  If Stretching kernel
      StretchingP2M();                                          //   Evaluate P2M kernel
    } else if ( kernelName == "Gaussian" ) {                    //  If Gaussian kernel
      GaussianP2M();                                            //   Evaluate P2M kernel
    } else {                                                    //  If kernel is none of the above
      if(MPIRANK == 0) std::cout << "Invalid kernel type" << std::endl;// Invalid kernel type
      abort();                                                  //   Abort execution
    }                                                           //  Endif for kernel type
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
      if( kernelName == "Laplace" ) {                           //   If Laplace kernel
        LaplaceM2M();                                           //    Evaluate M2M kernel
      } else if ( kernelName == "BiotSavart" ) {                //   If Biot Savart kernel
        BiotSavartM2M();                                        //    Evaluate M2M kernel
      } else if ( kernelName == "Stretching" ) {                //   If Stretching kernel
        StretchingM2M();                                        //    Evaluate M2M kernel
      } else if ( kernelName == "Gaussian" ) {                  //   If Gaussian kernel
        GaussianM2M();                                          //    Evaluate M2M kernel
      } else {                                                  //   If kernel is none of the above
        if(MPIRANK == 0) std::cout << "Invalid kernel type" << std::endl;// Invalid kernel type
        abort();                                                //   Abort execution
      }                                                         //   Endif for kernel type
      getTargetCell(listM2M);                                   //   Get body values from target buffer
      clearBuffers();                                           //   Clear GPU buffers
      ioffset += numCell;                                       //   Increment ioffset
    }                                                           //  End loop over icall
    level--;                                                    //  Decrement level
  }                                                             // End while loop over levels
}

void Evaluator::evalM2L(Cells &cells) {                         // Evaluate M2L
  CI0 = cells.begin();                                          // Set begin iterator
  const int numCell = MAXCELL/NCOEF/2;                          // Number of cells per icall
  int numIcall = int(cells.size()-1)/numCell+1;                 // Number of icall loops
  int ioffset = 0;                                              // Initialzie offset for icall loops
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
    if( kernelName == "Laplace" ) {                             //  If Laplace kernel
      LaplaceM2L();                                             //   Evaluate M2L kernel
    } else if ( kernelName == "BiotSavart" ) {                  //  If Biot Savart kernel
      BiotSavartM2L();                                          //   Evaluate M2L kernel
    } else if ( kernelName == "Stretching" ) {                  //  If Stretching kernel
      StretchingM2L();                                          //   Evaluate M2L kernel
    } else if ( kernelName == "Gaussian" ) {                    //  If Gaussian kernel
      GaussianM2L();                                            //   Evaluate M2L kernel
    } else {                                                    //  If kernel is none of the above
      if(MPIRANK == 0) std::cout << "Invalid kernel type" << std::endl;// Invalid kernel type
      abort();                                                  //   Abort execution
    }                                                           //  Endif for kernel type
    getTargetCell(listM2L,false);                               //  Get body values from target buffer
    clearBuffers();                                             //  Clear GPU buffers
    ioffset += numCell;                                         //  Increment ioffset
  }                                                             // End loop over icall
  listM2L.clear();                                              // Clear interaction lists
  flagM2L.clear();                                              // Clear periodic image flags
}

void Evaluator::evalM2P(Cells &cells) {                         // Evaluate M2P
  CI0 = cells.begin();                                          // Set begin iterator for target
  const int numCell = MAXCELL/NCRIT/7;                          // Number of cells per icall
  int numIcall = int(cells.size()-1)/numCell+1;                 // Number of icall loops
  int ioffset = 0;                                              // Initialzie offset for icall loops
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
    if( kernelName == "Laplace" ) {                             //  If Laplace kernel
      LaplaceM2P();                                             //   Evaluate M2P kernel
    } else if ( kernelName == "BiotSavart" ) {                  //  If Biot Savart kernel
      BiotSavartM2P();                                          //   Evaluate M2P kernel
    } else if ( kernelName == "Stretching" ) {                  //  If Stretching kernel
      StretchingM2P();                                          //   Evaluate M2P kernel
    } else if ( kernelName == "Gaussian" ) {                    //  If Gaussian kernel
      GaussianM2P();                                            //   Evaluate M2P kernel
    } else {                                                    //  If kernel is none of the above
      if(MPIRANK == 0) std::cout << "Invalid kernel type" << std::endl;// Invalid kernel type
      abort();                                                  //   Abort execution
    }                                                           //  Endif for kernel type
    getTargetBody(listM2P);                                     //  Get body values from target buffer
    clearBuffers();                                             //  Clear GPU buffers
    ioffset += numCell;                                         //  Increment ioffset
  }                                                             // End loop over icall
  listM2P.clear();                                              // Clear interaction lists
  flagM2P.clear();                                              // Clear periodic image flags
}

void Evaluator::evalP2P(Cells &cells) {                         // Evaluate P2P
  CI0 = cells.begin();                                          // Set begin iterator
  const int numCell = MAXCELL/NCRIT/7;                          // Number of cells per icall
  int numIcall = int(cells.size()-1)/numCell+1;                 // Number of icall loops
  int ioffset = 0;                                              // Initialzie offset for icall loops
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
    if( kernelName == "Laplace" ) {                             //  If Laplace kernel
      LaplaceP2P();                                             //   Evaluate P2P kernel
    } else if ( kernelName == "BiotSavart" ) {                  //  If Biot Savart kernel
      BiotSavartP2P();                                          //   Evaluate P2P kernel
    } else if ( kernelName == "Stretching" ) {                  //  If Stretching kernel
      StretchingP2P();                                          //   Evaluate P2P kernel
    } else if ( kernelName == "Gaussian" ) {                    //  If Gaussian kernel
      GaussianP2P();                                            //   Evaluate P2P kernel
    } else {                                                    //  If kernel is none of the above
      if(MPIRANK == 0) std::cout << "Invalid kernel type" << std::endl;// Invalid kernel type
      abort();                                                  //   Abort execution
    }                                                           //  Endif for kernel type
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
      if( kernelName == "Laplace" ) {                           //   If Laplace kernel
        LaplaceL2L();                                           //    Evaluate L2L kernel
      } else if ( kernelName == "BiotSavart" ) {                //   If Biot Savart kernel
        BiotSavartL2L();                                        //    Evaluate L2L kernel
      } else if ( kernelName == "Stretching" ) {                //   If Stretching kernel
        StretchingL2L();                                        //    Evaluate L2L kernel
      } else if ( kernelName == "Gaussian" ) {                  //   If Gaussian kernel
        GaussianL2L();                                          //    Evaluate L2L kernel
      } else {                                                  //   If kernel is none of the above
        if(MPIRANK == 0) std::cout << "Invalid kernel type" << std::endl;// Invalid kernel type
        abort();                                                //    Abort execution
      }                                                         //   Endif for kernel type
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
    if( kernelName == "Laplace" ) {                             //  If Laplace kernel
      LaplaceL2P();                                             //   Evaluate L2P kernel
    } else if ( kernelName == "BiotSavart" ) {                  //  If Biot Savart kernel
      BiotSavartL2P();                                          //   Evaluate L2P kernel
    } else if ( kernelName == "Stretching" ) {                  //  If Stretching kernel
      StretchingL2P();                                          //   Evaluate L2P kernel
    } else if ( kernelName == "Gaussian" ) {                    //  If Gaussian kernel
      GaussianL2P();                                            //   Evaluate L2P kernel
    } else {                                                    //  If kernel is none of the above
      if(MPIRANK == 0) std::cout << "Invalid kernel type" << std::endl;// Invalid kernel type
      abort();                                                  //   Abort execution
    }                                                           //  Endif for kernel type
    getTargetBody(listL2P);                                     //  Get body values from target buffer
    clearBuffers();                                             //  Clear GPU buffers
    ioffset += numCell;                                         //  Increment ioffset
  }                                                             // End loop over icall
}
