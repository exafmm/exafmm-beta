#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <algorithm>
#include <vector>
#include "pocsort.h"


/* test program for pocsort.h 

  try to sort a data of "item" below.

  pocsort can sort vector<T> of any type
  T that implements 
  (i) operator<   // >
  (ii) sort_key()

  sort_key() should be an unsinged 
  integer. 
*/

typedef unsigned long key_type;
template<int item_size>
struct item {
  union {
    key_type x;
    char pad[item_size];
  };
  item(key_type x_) : x(x_) {}
  key_type sort_key() { return x; }
};

template<int item_size>
bool operator<(const item<item_size>& a, const item<item_size>& b) {
  return a.x < b.x;
}

inline double cur_time() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return double(tv.tv_sec+tv.tv_usec*1e-6);
}

template<int item_size>
int test_sort(int key_size, long n) {
  printf("item_size: %lu\n", sizeof(item<item_size>));
  printf("number of items: %lu\n", n);
  printf("data_size: %lu\n", sizeof(item<item_size>) * n);
  /* make an array to sort */
  unsigned short xsub[3] = { 10, 20, 30 };
  std::vector<item<item_size> > A; // array to sort by stl
  std::vector<item<item_size> > S; // array to sort by pocsort
  std::vector<item<item_size> > T; // tmp array

  for (long i = 0; i < n; i++) {
    unsigned x = nrand48(xsub) & ((1L << n) - 1);
    A.push_back(item<item_size>(x));
  }
  S = A;
  T = A;
  pocsort<item<item_size> > ps;
  double t0 = cur_time();
  ps.sort_vector(S, T, key_size);
  double t1 = cur_time();
  std::sort(A.begin(), A.end());
  double t2 = cur_time();
  for (size_t i = 0; i < A.size(); i++) {
    if (A[i].sort_key() != S[i].sort_key()) {
      fprintf(stderr, "error: A[%ld].sort_key = %lu != S[%ld].sort_key() = %lu\n",
	      i, A[i].sort_key(), i, S[i].sort_key());
      exit(1);
    }
  }
  printf("OK %f sec by stl sort, %f sec by poc sort\n", t2 - t1, t1 - t0);
  return 0;
}

int main(int argc, char ** argv) {
  long n = (argc > 1 ? atol(argv[1]) : 10000);
  test_sort<64>(32, n);
  return 0;
}
