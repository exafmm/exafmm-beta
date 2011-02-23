#include "evaluator.h"

void Evaluator::evalP2P(Bodies &ibodies, Bodies &jbodies) {
  P2P(ibodies.begin(),ibodies.end(),jbodies.begin(),jbodies.end());
}

void Evaluator::evalP2M(Cells &cells) {                       // Evaluate P2M kernel
  for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {        // Loop over cells
    C->M = C->L = 0;                                          //  Initialize multipole & local coefficients
    P2M(C);                                                   //  Evaluate P2M kernel
  }                                                           // End loop over cells
}

void Evaluator::evalM2M(Cells &cells) {                       // Evaluate M2M kernel
  CI0 = cells.begin();                                        // Set begin iterator
  for( C_iter C=cells.begin(); C!=cells.end()-1; ++C ) {      // Loop over cells bottomup (except root cell)
    M2M(CI0+C->PARENT,C);                                     //  Evaluate M2M kernel
  }                                                           // End loop over cells
}

void Evaluator::evalM2L(Cells &cells) {                       // Evaluate M2L kernel
  CI0 = cells.begin();                                        // Set begin iterator
  for( C_iter CI=cells.begin(); CI!=cells.end(); ++CI ) {     // Loop over cells
    while( !listM2L[CI-CI0].empty() ) {                       //  While M2L interaction list is not empty
      C_iter CJ = listM2L[CI-CI0].back();                     //   Set source cell iterator
      M2L(CI,CJ);                                             //   Evaluate M2L kernel
      listM2L[CI-CI0].pop_back();                             //   Pop last element from M2L interaction list
    }                                                         //  End while for M2L interaction list
  }                                                           // End loop over cells topdown
}

void Evaluator::evalM2P(Cells &cells) {                       // Evaluate M2P kernel
  CI0 = cells.begin();                                        // Set begin iterator
  for( C_iter CI=cells.begin(); CI!=cells.end(); ++CI ) {     // Loop over cells
    while( !listM2P[CI-CI0].empty() ) {                       //  While M2P interaction list is not empty
      C_iter CJ = listM2P[CI-CI0].back();                     //   Set source cell iterator
      M2P(CI,CJ);                                             //   Evaluate M2P kernel
      listM2P[CI-CI0].pop_back();                             //   Pop last element from M2P interaction list
    }                                                         //  End while for M2P interaction list
  }                                                           // End loop over cells topdown
}

void Evaluator::evalP2P(Cells &cells) {                       // Evaluate P2P kernel
  CI0 = cells.begin();                                        // Set begin iterator
  for( C_iter CI=cells.begin(); CI!=cells.end(); ++CI ) {     // Loop over cells
    while( !listP2P[CI-CI0].empty() ) {                       //  While M2P interaction list is not empty
      C_iter CJ = listP2P[CI-CI0].back();                     //   Set source cell iterator
      P2P(CI->LEAF,CI->LEAF+CI->NLEAF,CJ->LEAF,CJ->LEAF+CJ->NLEAF);// Evaluate P2P kernel
      listP2P[CI-CI0].pop_back();                             //   Pop last element from M2P interaction list
    }                                                         //  End while for M2P interaction list
  }                                                           // End loop over cells topdown
}

void Evaluator::evalL2L(Cells &cells) {                       // Evaluate L2L kernel
  CI0 = cells.begin();                                        // Set begin iterator
  for( C_iter C=cells.end()-2; C!=cells.begin()-1; --C ) {    // Loop over cells topdown (except root cell)
    L2L(C,CI0+C->PARENT);                                     //  Evaluate L2L kernel
  }                                                           // End loop over cells topdown
}

void Evaluator::evalL2P(Cells &cells) {                       // Evaluate L2P kernel
  for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {        // Loop over cells
    if( C->NCHILD == 0 ) L2P(C);                              //  If cell is a twig evaluate L2P kernel
  }                                                           // End loop over cells topdown
}
