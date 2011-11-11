#ifndef serialfmm_h
#define serialfmm_h
#include "topdown.h"
#include "bottomup.h"

//! Serial FMM interface
class SerialFMM : public TopDown, public BottomUp {
public:
//! Constructor
  SerialFMM() : TopDown(), BottomUp() {
    preCalculation();
  }
//! Destructor
  ~SerialFMM() {
    postCalculation();
  }

//! Topdown tree constructor interface. Input: bodies, Output: cells
  void topdown(Bodies &bodies, Cells &cells) {
    TopDown::grow(bodies);                                      // Grow tree structure topdown

    TopDown::setIndex();                                        // Set index of cells

    buffer.resize(bodies.size());                               // Resize sort buffer
    sortBodies(bodies,buffer,false);                            // Sort bodies in descending order

    Cells twigs;                                                // Twigs are cells at the bottom of tree
    bodies2twigs(bodies,twigs);                                 // Turn bodies to twigs

    Cells sticks;                                               // Sticks are twigs from other processes that are not twigs in the current process
    twigs2cells(twigs,cells,sticks);                            // Turn twigs to cells
  }

//! Bottomup tree constructor interface. Input: bodies, Output: cells
  void bottomup(Bodies &bodies, Cells &cells) {
    BottomUp::setIndex(bodies);                                 // Set index of cells

    buffer.resize(bodies.size());                               // Resize sort buffer
    sortBodies(bodies,buffer,false);                            // Sort bodies in descending order

/*
    prune(bodies);                                              // Prune tree structure bottomup

    BottomUp::grow(bodies);                                     // Grow tree structure at bottom if necessary

    sortBodies(bodies,buffer,false);                            // Sort bodies in descending order
*/

    Cells twigs;                                                // Twigs are cells at the bottom of tree
    bodies2twigs(bodies,twigs);                                 // Turn bodies to twigs

    Cells sticks;                                               // Sticks are twigs from other processes not twigs here
    twigs2cells(twigs,cells,sticks);                            // Turn twigs to cells
  }
};

#endif
