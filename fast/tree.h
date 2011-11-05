#ifndef tree_h
#define tree_h
#include "bottomup.h"

class TreeConstructor : public BottomUp {
public:
  void topdown(Bodies &bodies, Cells &cells) {
    TOPDOWN = true;
    TopDown::setDomain(bodies);
    TopDown::buildTree();
    TopDown::linkTree(bodies,cells);
    TopDown::upwardPass(cells);
  }

  void bottomup(Bodies &bodies, Cells &cells) {
    TOPDOWN = false;
    BottomUp::setDomain(bodies);
    BottomUp::buildTree(bodies,cells);
    BottomUp::linkTree(cells);
    BottomUp::upwardPass(cells);
  }
};

class FastMultipoleMethod : public TreeConstructor {
public:
  void direct(Bodies &bodies, Cells &cells) {
    setRootCell(cells);
    P2P(ROOT);
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      B->TRG /= B->SRC[0];
    }
  }

  void direct2(Bodies &bodies, Cells &cells) {
    setRootCell(cells);
    P2P(ROOT,ROOT,false);
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      B->TRG /= B->SRC[0];
    }
  }

  void approximate(Cells &cells) {
    setRootCell(cells);
    startTimer("Traverse     ");
    traverse();
    stopTimer("Traverse     ",printNow);
    startTimer("Downward pass");
    if( TOPDOWN ) {
      TopDown::downwardPass(cells);
    } else {
      BottomUp::downwardPass(cells);
    }
    stopTimer("Downward pass",printNow);
    if(printNow) printTreeData(cells);
  }

  void approximate2(Cells &cells) {
    setRootCell(cells);
    startTimer("Traverse     ");
    traverse(false);
    stopTimer("Traverse     ",printNow);
    startTimer("Downward pass");
    if( TOPDOWN ) {
      TopDown::downwardPass(cells);
    } else {
      BottomUp::downwardPass(cells);
    }
    stopTimer("Downward pass",printNow);
    if(printNow) printTreeData(cells);
  }

};

#endif
