#ifndef evaluator_h
#define evaluator_h
#include "kernel.h"

class Evaluator {
private:
  C_iter CI0;                                                   // icells.begin()
  C_iter CJ0;                                                   // jcells.begin()
  Lists  listM2L;                                               // M2L interaction list
  Lists  listM2P;                                               // M2P interaction list
  Lists  listP2P;                                               // P2P interaction list
  Pairs  pairs;
  Kernel K;

  void tryM2L(C_iter CI, C_iter CJ) {                           // Interface for M2L kernel
    vect dist = CI->X - CJ->X;                                  // Distance vector between cells
    real R = std::sqrt(norm(dist));                             // Distance between cells
    if( CI->R + CJ->R > THETA*R ) {                             // If cell is too large
      Pair pair(CI,CJ);                                         //  Form pair of interacting cells
      pairs.push(pair);                                         //  Push interacting pair into stack
    } else {                                                    // If cell is small enough
      listM2L[CI-CI0].push_back(CJ);                            // Push source cell into M2L interaction list
    }                                                           // Endif for interaction
  }

  void tryM2P(C_iter CI, C_iter CJ) {                           // Interface for M2P kernel
    vect dist = CI->X - CJ->X;                                  // Distance vector between cells
    real R = std::sqrt(norm(dist));                             // Distance between cells
    if( CI->NCHILD != 0 || CI->R + CJ->R > THETA*R ) {          // If target is not twig or cell is too large
      Pair pair(CI,CJ);                                         //  Form pair of interacting cells
      pairs.push(pair);                                         //  Push interacting pair into stack
    } else {                                                    // If target is twig and cell is small enough
      listM2P[CI-CI0].push_back(CJ);                            // Push source cell into M2P interaction list
    }                                                           // Endif for interaction
  }

  void treecode(C_iter CI, C_iter CJ) {                         // Tree walk for treecode
    if( CI->NCHILD == 0 && CJ->NCHILD == 0) {                   // If both cells are twigs
      if( CJ->NLEAF != 0 ) {                                    // If the twig has leafs
        listP2P[CI-CI0].push_back(CJ);                          // Push source cell into P2P interaction list
      } else {                                                  // If the twig has no leafs
#ifdef DEBUG
        std::cout << "CJ->I=" << CJ->I << " has no leaf. Doing M2P instead." << std::endl;
#endif
        listM2P[CI-CI0].push_back(CJ);                          // Push source cell into M2P interaction list
      }                                                         // Endif for twigs with leafs
    } else if ( CI->NCHILD != 0 ) {                             // If target is not twig
      for( int i=0; i<CI->NCHILD; i++ ) {                       //  Loop over child cells of target
        tryM2P(CI0+CI->CHILD[i],CJ);                            //   Try to evaluate M2P kernel
      }                                                         //  End loop over child cells of target
    } else {                                                    // If target is twig
      for( int i=0; i<CJ->NCHILD; i++ ) {                       //  Loop over child cells of source
        tryM2P(CI,CJ0+CJ->CHILD[i]);                            //   Try to evaluate M2P kernel
      }                                                         //  End loop over child cells of source
    }                                                           // Endif for type of interaction
  }

  void FMM(C_iter CI, C_iter CJ) {                              // Tree walk for FMM
    if( CI->NCHILD == 0 && CJ->NCHILD == 0 ) {                  // If both cells are twigs
      if( CJ->NLEAF != 0 ) {                                    // If the twig has leafs
        listP2P[CI-CI0].push_back(CJ);                          // Push source cell into P2P interaction list
      } else {                                                  // If the twig has no leafs
#ifdef DEBUG
        std::cout << "CJ->I=" << CJ->I << " has no leaf. Doing M2P instead." << std::endl;
#endif
        listM2P[CI-CI0].push_back(CJ);                          // Push source cell into M2P interaction list
      }                                                         // Endif for twigs with leafs
    } else if ( CJ->NCHILD == 0 || (CI->NCHILD != 0 && CI->R > CJ->R) ) {// If source is twig or target is larger
      for( int i=0; i<CI->NCHILD; i++ ) {                       //  Loop over child cells of target
        tryM2L(CI0+CI->CHILD[i],CJ);                            //   Try to evaluate M2L kernel
      }                                                         //  End loop over child cells of target
    } else {                                                    // If target is twig or source is larger
      for( int i=0; i<CJ->NCHILD; i++ ) {                       //  Loop over child cells of source
        tryM2L(CI,CJ0+CJ->CHILD[i]);                            //   Try to evaluate M2L kernel
      }                                                         //  End loop over child cells of source
    }                                                           // Endif for type of interaction
  }

public:
  void initialize() {
    K.initialize();
  }

  void traverse(Cells &cells, Cells &jcells, int method) {      // Traverse tree to get interaction list
    C_iter root = cells.end()-1;                                // Iterator for root cell
    C_iter jroot = jcells.end()-1;                              // Iterator for root cell
    CI0 = cells.begin();                                        // Set begin iterator for icells
    CJ0 = jcells.begin();                                       // Set begin iterator for jcells
    Pair pair(root,jroot);                                      // Form pair of root cells
    pairs.push(pair);                                           // Push pair to stack
    listM2L.resize(cells.size());                               // Resize M2L interaction list
    listM2P.resize(cells.size());                               // Resize M2P interaction list
    listP2P.resize(cells.size());                               // Resize P2P interaction list
    while( !pairs.empty() ) {                                   // While interaction stack is not empty
      pair = pairs.top();                                       //  Get interaction pair from top of stack
      pairs.pop();                                              //  Pop interaction stack
      switch (method) {                                         //  Swtich between methods
        case 0 : treecode(pair.CI,pair.CJ); break;              //   0 : treecode
        case 1 : FMM(pair.CI,pair.CJ);      break;              //   1 : FMM
      }                                                         //  End switch between methods
    }                                                           // End while loop for interaction stack
  }

  void P2P(Cells &cells) {                                      // Evaluate P2P kernel
    for( C_iter CI=cells.begin(); CI!=cells.end(); ++CI ) {     // Loop over cells
      while( !listP2P[CI-CI0].empty() ) {                       //  While M2P interaction list is not empty
        C_iter CJ = listP2P[CI-CI0].back();                     //   Set source cell iterator
        K.P2P(CI->LEAF,CI->LEAF+CI->NLEAF,CJ->LEAF,CJ->LEAF+CJ->NLEAF);// Evaluate P2P kernel
        listP2P[CI-CI0].pop_back();                             //   Pop last element from M2P interaction list
      }                                                         //  End while for M2P interaction list
    }                                                           // End loop over cells topdown
  }

  void P2M(Cells &twigs) {                                      // Evaluate P2M kernel
    for( C_iter C=twigs.begin(); C!=twigs.end(); ++C ) {        // Loop over twigs
      C->M = C->L = 0;                                          //  Initialize multipole/local coefficients
      K.P2M(C);                                                 //  Evaluate P2M kernel
    }                                                           // End loop over cells
  }

  void M2M(Cells &cells) {                                      // Evaluate M2M kernel
    for( C_iter C=cells.begin(); C!=cells.end()-1; ++C ) {      // Loop over cells bottomup (except root cell)
      K.M2M(cells.begin()+C->PARENT,C);                         //  Evaluate M2M kernel
    }                                                           // End loop over cells
  }

  void M2L(Cells &cells) {                                      // Evaluate M2L kernel
    for( C_iter CI=cells.begin(); CI!=cells.end(); ++CI ) {     // Loop over cells
      while( !listM2L[CI-CI0].empty() ) {                       //  While M2L interaction list is not empty
        C_iter CJ = listM2L[CI-CI0].back();                     //   Set source cell iterator
        K.M2L(CI,CJ);                                           //   Evaluate M2L kernel
        listM2L[CI-CI0].pop_back();                             //   Pop last element from M2L interaction list
      }                                                         //  End while for M2L interaction list
    }                                                           // End loop over cells topdown
  }

  void L2L(Cells &cells) {                                      // Evaluate L2L kernel
    for( C_iter C=cells.end()-2; C!=cells.begin()-1; --C ) {    // Loop over cells topdown (except root cell)
      K.L2L(C,CI0+C->PARENT);                                   //  Evaluate L2L kernel
    }                                                           // End loop over cells topdown
  }

  void L2P(Cells &cells) {                                      // Evaluate L2P kernel
    for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {        // Loop over cells
      if( C->NCHILD == 0 ) K.L2P(C);                            //  If cell is a twig evaluate L2P kernel
    }                                                           // End loop over cells topdown
  }

  void M2P(Cells &cells) {                                      // Evaluate M2P kernel
    for( C_iter CI=cells.begin(); CI!=cells.end(); ++CI ) {     // Loop over cells
      while( !listM2P[CI-CI0].empty() ) {                       //  While M2P interaction list is not empty
        C_iter CJ = listM2P[CI-CI0].back();                     //   Set source cell iterator
        K.M2P(CI,CJ);                                           //   Evaluate M2P kernel
        listM2P[CI-CI0].pop_back();                             //   Pop last element from M2P interaction list
      }                                                         //  End while for M2P interaction list
    }                                                           // End loop over cells topdown
  }

  void finalize() {
    K.finalize();
  }
};

#endif
