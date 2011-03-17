#ifndef kernel_h
#define kernel_h
#include "logger.h"

class Kernel : public Logger {                                  // Unified CPU/GPU kernel class
protected:
  B_iter BI0;                                                   // Target bodies begin iterator
  B_iter BIN;                                                   // Target bodies end iterator
  B_iter BJ0;                                                   // Source bodies begin iterator
  B_iter BJN;                                                   // Source bodies end iterator
  C_iter CI;                                                    // Target cell iterator
  C_iter CJ;                                                    // Source cell iterator
  vect   Xperiodic;                                             // Coordinate offset of periodic image

  std::vector<int>   keysHost;                                  // Offsets for rangeHost
  std::vector<int>   rangeHost;                                 // Offsets for sourceHost
  std::vector<float> constHost;                                 // Constants on host
  std::vector<float> sourceHost;                                // Sources on host
  std::vector<float> targetHost;                                // Targets on host
  Map                sourceBegin;                               // Define map for offset of source cells
  Map                sourceSize;                                // Define map for size of source cells
  Map                targetBegin;                               // Define map for offset of target cells

  int getLevel(bigint index) {                                  // Get level from cell index
    int level = -1;                                             // Initialize level counter
    while( index >= 0 ) {                                       // While cell index is non-negative
      level++;                                                  //  Increment level
      index -= 1 << 3*level;                                    //  Subtract number of cells in that level
    }                                                           // End while loop for cell index
    return level;                                               // Return the level
  }

public:
  void initialize();

  void allocGPU();

  void deallocGPU();

  void P2M();

  void M2M();

  void M2L();

  void M2P();

  void P2P();

  void L2L();

  void L2P();

  void finalize();
};

#endif
