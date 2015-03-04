#ifdef CILK

#if defined TBB || defined MTHREAD || defined QTHREAD || defined SERIAL
#include "build_tree_cilk.h"
#else
#include "build_tree_tbb.h"
#endif

#elif defined TBB || defined MTHREAD || defined QTHREAD || defined SERIAL
#include "build_tree_tbb.h"

#else
#include "build_tree_omp.h"

#endif
