#ifndef tree_h
#define tree_h
#include "bottomup.h"

class TreeConstructor : public BottomUp {
private:
  bool TOPDOWN;

  void printTreeData() const {
    std::cout << "------------------------" << std::endl;
    std::cout << "Root center   : " << ROOT->X              << std::endl;
    std::cout << "Root radius   : " << R0                   << std::endl;
    std::cout << "Bodies        : " << ROOT->NDLEAF         << std::endl;
    std::cout << "Cells         : " << NCELL                << std::endl;
    std::cout << "Tree depth    : " << MAXLEVEL             << std::endl;
    std::cout << "Total charge  : " << std::abs(ROOT->M[0]) << std::endl;
    std::cout << "------------------------" << std::endl;
  }

public:
  void topdown(Bodies &bodies) {
    TOPDOWN = true;
    TopDown::setDomain(bodies);
    TopDown::build();
    TopDown::link(bodies);
    TopDown::upward();
  }

  void bottomup(Bodies &bodies) {
    TOPDOWN = false;
    BottomUp::setDomain(bodies);
    BottomUp::build(bodies);
    BottomUp::link();
    BottomUp::upward();
  }

  void exact(Bodies &bodies, bool IeqJ=true) const {
    if( IeqJ ) {
      P2P(ROOT);
    } else {
//      P2P(ROOT,ROOT,false);
    }
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      B->TRG /= B->SRC[0];
    }
  }

  void approximate(bool IeqJ=true) {
    startTimer("Traverse     ");
    if( IeqJ ) {
      traverse();
    } else {
//      traverse(false);
    }
    stopTimer("Traverse     ",printNow);
    startTimer("Downward     ");
    if( TOPDOWN ) {
      TopDown::downward();
    } else {
      BottomUp::downward();
    }
    stopTimer("Downward     ",printNow);
    if(printNow) printTreeData();
    cells.clear();
  }
};

#endif
