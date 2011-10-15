#ifndef tree_h
#define tree_h
#include "bottomup.h"

class TreeConstructor : public BottomUp {
private:
  void write() const {
    std::cout<<" root center:           "<<CN->X            <<'\n';
    std::cout<<" root radius:           "<<R0               <<'\n';
    std::cout<<" bodies loaded:         "<<CN->NDLEAF       <<'\n';
    std::cout<<" total scal:            "<<CN->M[0]         <<'\n';
    std::cout<<" cells used:            "<<NCELL            <<'\n';
    std::cout<<" maximum level:         "<<MAXLEVEL         <<'\n';
  }

public:
  void topdown(Bodies &bodies) {
    TopDown::setDomain(bodies);
    TopDown::build();
    TopDown::link(bodies);
    startTimer("Upward       ");
    TopDown::upward();
    stopTimer("Upward       ",printNow);
  }

  void bottomup(Bodies &bodies) {
    BottomUp::setDomain(bodies);
    BottomUp::build(bodies);
    BottomUp::link();
    startTimer("Upward       ");
    BottomUp::upward();
    stopTimer("Upward       ",printNow);
  }

  void exact(Bodies &bodies) {
#if IEQJ
    P2P(CN);
#else
    P2P(CN,CN,false);
#endif
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      B->TRG /= B->SRC[0];
    }
  }

  void approximate() {
    startTimer("Traverse     ");
#if IEQJ
    traverse();
#else
    traverse(false);
#endif
    stopTimer("Traverse     ",printNow);
    startTimer("Downward     ");
    for( C_iter C=C0+CN->CHILD; C!=C0+CN->CHILD+CN->NCHILD; ++C ) {
      downward(C);
    }
#ifndef MANY
    write();
#endif
    cells.clear();
    stopTimer("Downward     ",printNow);
  }
};

#endif
