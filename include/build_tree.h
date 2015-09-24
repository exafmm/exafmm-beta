#ifndef build_tree_h
#define build_tree_h

#include "config.h"

#ifdef WITH_CILK

#if defined WITH_TBB || defined WITH_MTHREAD || defined WITH_QTHREAD || defined DISABLE_THREAD
#include "build_tree_cilk.h"
#else
#include "build_tree_tbb.h"
#endif

#elif defined WITH_TBB || defined WITH_MTHREAD || defined WITH_QTHREAD || defined DISABLE_THREAD
#include "build_tree_tbb.h"

#else
#include "build_tree_omp.h"

#endif

#endif
