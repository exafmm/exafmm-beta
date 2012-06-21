template<Equation equation>
void Evaluator<equation>::timeKernels() {
  Bodies ibodies(1000), jbodies(1000);
  for( B_iter Bi=ibodies.begin(),Bj=jbodies.begin(); Bi!=ibodies.end(); ++Bi, ++Bj ) {
    Bi->X = 0;
    Bj->X = 1;
  }
  Cells cells;
  cells.resize(2);
  C_iter Ci = cells.begin(), Cj = cells.begin()+1;
  Ci->X = 0;
  Ci->NDLEAF = 10;
  Ci->LEAF = ibodies.begin();
  Ci->M = 0;
  Ci->L = 0;
  Cj->X = 1;
  Cj->NDLEAF = 1000;
  Cj->LEAF = jbodies.begin();
  Cj->M = 0;
  startTimer("P2P kernel");
  P2P(Ci,Cj);
  timeP2P = stopTimer("P2P kernel") / 10000;
  startTimer("M2L kernel");
  for( int i=0; i!=1000; ++i ) M2L(Ci,Cj);
  timeM2L = stopTimer("M2L kernel") / 1000;
  startTimer("M2P kernel");
  for( int i=0; i!=100; ++i ) M2P(Ci,Cj);
  timeM2P = stopTimer("M2P kernel") / 1000;
}

template<Equation equation>
inline void Evaluator<equation>::direct(Bodies &ibodies, Bodies &jbodies) {// Evaluate direct summation
  Cells cells;
  cells.resize(2);
  C_iter Ci = cells.begin(), Cj = cells.begin()+1;
  Ci->LEAF = ibodies.begin();
  Ci->NDLEAF = ibodies.size();
  Cj->LEAF = jbodies.begin();
  Cj->NDLEAF = jbodies.size();
  P2P(Ci,Cj);
}

template<Equation equation>
inline void Evaluator<equation>::evalP2M(Cells &cells) {        // Evaluate all P2M kernels
  if( TOPDOWN ) {                                               // If tree was constructed top down
    for( C_iter C=cells.end()-1; C!=cells.begin()-1; --C ) {    //  Loop over cells forwards
      if( C->NCHILD == 0 ) P2M(C);                              //   If cell is a twig do P2M
    }                                                           //  End loop over cells
  } else {                                                      // If tree was constructed bottom up
    for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {        //  Loop over cells backwards
      if( C->NCHILD == 0 ) P2M(C);                              //   If cell is a twig do P2M
    }                                                           //  End loop over cells
  }                                                             // Endif for tree construction
}

template<Equation equation>
inline void Evaluator<equation>::evalM2M(Cells &cells, Cells &jcells) {// Evaluate all M2M kernels
  Cj0 = jcells.begin();                                         // Set begin iterator
  if( TOPDOWN ) {                                               // If tree was constructed top down
    for( C_iter C=cells.end()-1; C!=cells.begin()-1; --C ) {    //  Loop over cells forwards
      if( C->NCHILD != 0 ) M2M(C);                              //   If cell is not a twig do M2M
    }                                                           //  End loop over cells
    cells.front().X = X0;                                       //  Set root center
    cells.front().R = R0;                                       //  Set root radius
  } else {                                                      // If tree was constructed bottom up
    for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {        //  Loop over cells backwards
      if( C->NCHILD != 0 ) M2M(C);                              //   If cell is not a twig do M2M
    }                                                           //  End loop over cells
    cells.back().X = X0;                                        //  Set root center
    cells.back().R = R0;                                        //  Set root radius
  }                                                             // Endif for tree construction
}

template<Equation equation>
void Evaluator<equation>::evalM2L(C_iter Ci, C_iter Cj) {       // Evaluate single M2L kernel
  M2L(Ci,Cj);                                                   // Call M2L kernel
//  NM2L++;                                                       // Count M2L kernel execution
}

template<Equation equation>
void Evaluator<equation>::evalM2P(C_iter Ci, C_iter Cj) {       // Evaluate single M2P kernel
  M2P(Ci,Cj);                                                   // Call M2P kernel
//  NM2P++;                                                       // Count M2P kernel execution
}

template<Equation equation>
void Evaluator<equation>::evalP2P(C_iter Ci, C_iter Cj) {       // Evaluate single P2P kernel
  P2P(Ci,Cj);                                                   // Call P2P kernel
//  NP2P++;                                                       // Count P2P kernel execution
}

template<Equation equation>
inline void Evaluator<equation>::evalL2L(Cells &cells) {        // Evaluate all L2L kernels
  if( TOPDOWN ) {                                               // If tree was constructed top down
    for( C_iter C=cells.begin()+1; C!=cells.end(); ++C ) {      //  Loop over cells forwards
      L2L(C);                                                   //   Do L2L
    }                                                           //  End loop over cells
  } else {                                                      // If tree was constructed bottom up
    for( C_iter C=cells.end()-2; C!=cells.begin()-1; --C ) {    //  Loop over cells backwards
      L2L(C);                                                   //   Do L2L
    }                                                           //  End loop over cells
  }                                                             // Endif for tree construction
}

template<Equation equation>
inline void Evaluator<equation>::evalL2P(Cells &cells) {        // Evaluate all L2P kernels
  if( TOPDOWN ) {                                               // If tree was constructed top down
    for( C_iter C=cells.begin()+1; C!=cells.end(); ++C ) {      //  Loop over cells forwards
      if( C->NCHILD == 0 ) L2P(C);                              //   If cell is a twig do L2P
    }                                                           //  End loop over cells
  } else {                                                      // If tree was constructed bottom up
    for( C_iter C=cells.end()-2; C!=cells.begin()-1; --C ) {    //  Loop over cells backwards
      L2P(C);                                                   //   Do L2P
    }                                                           //  End loop over cells
  }                                                             // Endif for tree construction
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
  privateQueue.push_back(pair);
  while( !privateQueue.empty() ) {
    Pair Cij = privateQueue.front();
    privateQueue.pop_front();
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
