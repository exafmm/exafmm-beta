#include "evaluator.h"

void Evaluator::evalP2P(Bodies &ibodies, Bodies &jbodies) {
  B_iter BI0 = ibodies.begin();
  B_iter BIN = ibodies.end();
  B_iter BJ0 = jbodies.begin();
  B_iter BJN = jbodies.end();

  hostConstant[0] = BJN - BJ0;
  int NI = ((BIN - BI0) / THREADS + 1) * THREADS;
  int NJ = ((BJN - BJ0) / THREADS + 1) * THREADS;
  float4 *sourceHost, *sourceDevc;
  float  *targetHost, *targetDevc;

  targetHost = (float *)     malloc( NI*sizeof(float ) );
  sourceHost = (float4*)     malloc( NJ*sizeof(float4) );
  cudaMalloc(  (void**) &targetDevc, NI*sizeof(float ) );
  cudaMalloc(  (void**) &sourceDevc, NJ*sizeof(float4) );

  for( B_iter B=BJ0; B!=BJN; ++B ) {
    sourceHost[B-BJ0].x = B->pos[0];
    sourceHost[B-BJ0].y = B->pos[1];
    sourceHost[B-BJ0].z = B->pos[2];
    sourceHost[B-BJ0].w = B->scal;
  }

  cudaMemcpyToSymbol(deviceConstant,hostConstant,4*sizeof(int));
  cudaMemcpy(sourceDevc,sourceHost,Nround*sizeof(float4),cudaMemcpyHostToDevice);
  P2P(sourceDevc,targetDevc);
  cudaMemcpy(targetHost,targetDevc,Nround*sizeof(float ),cudaMemcpyDeviceToHost);

  for( B_iter B=BI0; B!=BIN; ++B ) {
    B->pot = targetHost[B-BI0];
  }
}

void Evaluator::evalP2M(Cells &cells) {                       // Evaluate P2M
  for( CJ=cells.begin(); CJ!=cells.end(); ++CJ ) {            // Loop over cells
    CJ->M = CJ->L = 0;                                        //  Initialize multipole & local coefficients
    P2M();                                                    //  Evaluate P2M kernel
  }                                                           // End loop over cells
}

void Evaluator::evalM2M(Cells &cells) {                       // Evaluate M2M
  CJ0 = cells.begin();                                        // Set begin iterator
  for( CJ=cells.begin(); CJ!=cells.end()-1; ++CJ ) {          // Loop over cells bottomup (except root cell)
    CI = CJ0+CJ->PARENT;                                      //  Set target cell iterator
    M2M();                                                    //  Evaluate M2M kernel
  }                                                           // End loop over cells
}

void Evaluator::evalM2L(Cells &cells) {                       // Evaluate M2L
  CI0 = cells.begin();                                        // Set begin iterator
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {            // Loop over cells
    while( !listM2L[CI-CI0].empty() ) {                       //  While M2L interaction list is not empty
      CJ = listM2L[CI-CI0].back();                            //   Set source cell iterator
      M2L();                                                  //   Evaluate M2L kernel
      listM2L[CI-CI0].pop_back();                             //   Pop last element from M2L interaction list
    }                                                         //  End while for M2L interaction list
  }                                                           // End loop over cells topdown
}

void Evaluator::evalM2P(Cells &cells) {                       // Evaluate M2P
  CI0 = cells.begin();                                        // Set begin iterator
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {            // Loop over cells
    while( !listM2P[CI-CI0].empty() ) {                       //  While M2P interaction list is not empty
      CJ = listM2P[CI-CI0].back();                            //   Set source cell iterator
      M2P();                                                  //   Evaluate M2P kernel
      listM2P[CI-CI0].pop_back();                             //   Pop last element from M2P interaction list
    }                                                         //  End while for M2P interaction list
  }                                                           // End loop over cells topdown
}

void Evaluator::evalP2P(Cells &cells) {                       // Evaluate P2P
  CI0 = cells.begin();                                        // Set begin iterator
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {            // Loop over cells
    while( !listP2P[CI-CI0].empty() ) {                       //  While M2P interaction list is not empty
      CJ = listP2P[CI-CI0].back();                            //   Set source cell iterator
      BI0 = CI->LEAF;                                         //   Set target bodies begin iterator
      BIN = CI->LEAF+CI->NLEAF;                               //   Set target bodies end iterator
      BJ0 = CJ->LEAF;                                         //   Set source bodies begin iterator
      BJN = CJ->LEAF+CJ->NLEAF;                               //   Set source bodies end iterator
      P2P();                                                  //   Evaluate P2P kernel
      listP2P[CI-CI0].pop_back();                             //   Pop last element from M2P interaction list
    }                                                         //  End while for M2P interaction list
  }                                                           // End loop over cells topdown
}

void Evaluator::evalL2L(Cells &cells) {                       // Evaluate L2L
  CI0 = cells.begin();                                        // Set begin iterator
  for( CI=cells.end()-2; CI!=cells.begin()-1; --CI ) {        // Loop over cells topdown (except root cell)
    CJ = CI0+CI->PARENT;                                      //  Set source cell iterator
    L2L();                                                    //  Evaluate L2L kernel
  }                                                           // End loop over cells topdown
}

void Evaluator::evalL2P(Cells &cells) {                       // Evaluate L2P
  for( CI=cells.begin(); CI!=cells.end(); ++CI ) {            // Loop over cells
    if( CI->NCHILD == 0 ) L2P();                              //  If cell is a twig evaluate L2P kernel
  }                                                           // End loop over cells topdown
}
