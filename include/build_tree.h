#ifndef build_tree_h
#define build_tree_h

#include "config.h"

#ifdef EXAFMM_WITH_CILK

#if defined EXAFMM_WITH_TBB || defined EXAFMM_WITH_MTHREAD || defined EXAFMM_WITH_QTHREAD
#include "build_tree_cilk.h"
#else
#include "build_tree_tbb.h"
#endif

#elif defined EXAFMM_WITH_TBB || defined EXAFMM_WITH_MTHREAD || defined EXAFMM_WITH_QTHREAD
#include "build_tree_tbb.h"

#else
#include "build_tree_omp.h"

#endif

#endif
