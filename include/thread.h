#ifndef thread_h
#define thread_h

#if TBB
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

#elif OPENMP
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

#elif MTHREAD
/* MassiveThreads (TBB-like interface on top of MassiveThreads)  */
#define num_threads(E)		myth_init_withparam(E, 1 << 16)
#define TO_MTHREAD_NATIVE 1
#include <tpswitch/tpswitch.h>

#elif QTHREAD
/* Qthreads (TBB-like interface on top of MassiveThreads).
   make sure you set QTHREAD_STACK_SIZE large enough (e.g., 131072) */
#define num_threads(E)		do { char n[30]; sprintf(n, "%d", E); setenv("QTHREAD_STACK_SIZE", "65536", 0); setenv("QTHREAD_NUM_SHEPHERDS", n, 0); setenv("QTHREAD_NUM_WORKERS_PER_SHEPHERD", "1", 0); qthread_initialize(); } while(0)
#define TO_QTHREAD 1
#include <tpswitch/tpswitch.h>

#elif SERIAL
/* Qthreads (TBB-like interface on top of MassiveThreads)  */
#define num_threads(E)
#define TO_SERIAL 1
#include <tpswitch/tpswitch.h>

#else
#define mk_task_group
#define wait_tasks
#define create_taskc(E)               E()
#define create_taskc_if(x, E)         E()
#define num_threads(E)

#endif

#endif
