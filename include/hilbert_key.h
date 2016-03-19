#ifndef HILBERT_KEY_H
#define HILBERT_KEY_H
//+++++++++++++++++++++++++++ PUBLIC-DOMAIN SOFTWARE ++++++++++++++++++++++++++
// Functions: TransposetoAxes AxestoTranspose
// Purpose: Transform in-place between Hilbert transpose and geometrical axes
// Example: b=5 bits for each of n=3 coordinates.
// 15-bit Hilbert integer = A B C D E F G H I J K L M N O is stored
// as its Transpose
// X[0] = A D G J M X[2]|
// X[1] = B E H K N <-------> | /X[1]
// X[2] = C F I L O axes |/
// high low 0------ X[0]
// Axes are stored conventially as b-bit integers.
// Author: John Skilling 20 Apr 2001 to 11 Oct 2003
//-----------------------------------------------------------------------------
#define cast_coord(V) static_cast<int>(V)

//! Transform in-place between geometrical axes and Hilbert transpose
inline void axesToTranspose(int X[], int order) {
  int M = 1 << (order - 1), P, Q, t;
  // Inverse undo
  for (Q = M; Q > 1; Q >>= 1) {
    P = Q - 1;
    for (int i = 0; i < 3; i++)
      if (X[i] & Q)
        X[0] ^= P;
      else {
        t = (X[0] ^ X[i]) & P;
        X[0] ^= t;
        X[i] ^= t;
      }
  }
  // Gray encode
  for (int i = 1; i < 3; i++) X[i] ^= X[i - 1];
  t = 0;
  for (Q = M; Q > 1; Q >>= 1) {
    if (X[2] & Q)
      t ^= Q - 1;
  }
  for (int i = 0; i < 3; i++) X[i] ^= t;
}

//! Output one 64-bit Hilbert order from the 3D transposed key
inline int64_t flattenTransposedKey(int X[], int order) {
  int64_t key = 0;
  int shifts = order - 1;
  for (int i = shifts; i >= 0; --i) {
    for (int j = 0; j < 3; ++j) {
      if (X[j] >> i & 1)
        key |= 1ull << ((3 * i) + (3 - j - 1));
    }
  }
  return key;
}

//! swap a & b values
void swap(int & a, int & b) {
  int c(a); a = b; b = c;
};

//! get linear Hilbert address given the x,y,z coordinates
//! migrated from grouptargets.h, but not tested yet
uint64_t getHilbert(int iX[3], int nbits) {
  const int octantMap[8] = {0, 1, 7, 6, 3, 2, 4, 5};
  int mask = 1 << (nbits - 1);
  uint64_t key = 0;
#pragma unroll
  for (int i = 0; i < nbits; i++) {
    const int ix = (iX[0] & mask) ? 1 : 0;
    const int iy = (iX[1] & mask) ? 1 : 0;
    const int iz = (iX[2] & mask) ? 1 : 0;
    const int octant = (ix << 2) + (iy << 1) + iz;
    if (octant == 0) {
      swap(iX[1], iX[2]);
    } else if (octant == 1 || octant == 5) {
      swap(iX[0], iX[1]);
    } else if (octant == 4 || octant == 6) {
      iX[0] = (iX[0]) ^ (-1);
      iX[2] = (iX[2]) ^ (-1);
    } else if (octant == 3 || octant == 7) {
      iX[0] = (iX[0]) ^ (-1);
      iX[1] = (iX[1]) ^ (-1);
      swap(iX[0], iX[1]);
    } else {
      iX[1] = (iX[1]) ^ (-1);
      iX[2] = (iX[2]) ^ (-1);
      swap(iX[1], iX[2]);
    }
    key = (key << 3) + octantMap[octant];
    mask >>= 1;
  }
  return key;
}

#endif
