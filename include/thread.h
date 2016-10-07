#ifndef thread_h
#define thread_h
#include "config.h"

#if EXAFMM_WITH_CILK
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#endif

#if EXAFMM_WITH_TBB
#define num_threads(E)                tbb::task_scheduler_init init(E)
#if DAG_RECORDER == 2  /* TBB with DAG Recorder */
#define TO_TBB 1 
#include <tpswitch/tpswitch.h>
#include <tbb/task_scheduler_init.h>
using namespace mtbb;
#else  /* TBB without DAG Recorder; original TBB */
#include <tbb/task_group.h>
#include <tbb/task_scheduler_init.h>
using namespace tbb;
#define mk_task_group                 task_group tg;
#define wait_tasks                    tg.wait()
#define create_taskc(E)               tg.run(E)
#define create_taskc_if(x, E)         if(x) { create_taskc(E); } else { E(); }
#endif

#elif EXAFMM_WITH_MTHREAD
/* MassiveThreads (TBB-like interface on top of MassiveThreads)  */
//#define num_threads(E)		      myth_init_ex(E, 1 << 16)
#define num_threads(E)		      do { myth_globalattr_set_stacksize(0, 1 << 16); myth_globalattr_set_n_workers(0, E); myth_init(); } while (0)
#define TO_MTHREAD_NATIVE 1
#include <tpswitch/tpswitch.h>

#elif EXAFMM_WITH_CILK
#define num_threads(E)                char nworkers[32]; sprintf(nworkers,"%d",E); __cilkrts_set_param("nworkers",nworkers)
#define mk_task_group
#define wait_tasks                    cilk_sync
template<class Call>
void call(Call C) { C(); }
#define create_taskc(E)               cilk_spawn call(E)
#define create_taskc_if(x, E)         if(x) { create_taskc(E); } else { E(); }

#elif EXAFMM_WITH_OPENMP
#include <omp.h>
#define num_threads(E)                omp_set_num_threads(E);
#if DAG_RECORDER == 2		/* OpenMP with DAG Recorder */
#define TO_OMP 1
#include <tpswitch/tpswitch.h>
#else  /* OpenMP without DAG Recorder; original OpenMP */
#define PRAGMA_OMP(x)                 _Pragma( #x )
#define mk_task_group
#define wait_tasks                    PRAGMA_OMP(omp taskwait)
#define create_taskc(E)               PRAGMA_OMP(omp task) E()
#define create_taskc_if(x, E)         if(x) { create_taskc(E); } else { E(); }
#endif

#else
#define mk_task_group
#define wait_tasks
#define create_taskc(E)               E()
#define create_taskc_if(x, E)         E()
#define num_threads(E)

#endif

#endif
