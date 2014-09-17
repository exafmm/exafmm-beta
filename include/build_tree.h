#ifdef CILK

#if defined TBB || defined MTHREAD || defined QTHREAD || defined SERIAL
#pragma message "Using build_tree_cilk.h"
#include "build_tree_cilk.h"
#else
#pragma message "Using build_tree_tbb.h"
#include "build_tree_tbb.h"
#endif

#elif defined TBB || defined MTHREAD || defined QTHREAD || defined SERIAL
#pragma message "Using build_tree_tbb.h"
#include "build_tree_tbb.h"

#else
#pragma message "Using build_tree_omp.h"
#include "build_tree_omp.h"

#endif
