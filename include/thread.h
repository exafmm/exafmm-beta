#ifndef thread_h
#define thread_h

#if TBB
#include <tbb/task_group.h>
using namespace tbb;
#elif MTHREAD
#define TO_MTHREAD_NATIVE 1
#include <tpswitch/tpswitch.h>
#elif QTHREAD
#define TO_QTHREAD 1
#include <tpswitch/tpswitch.h>
#endif

#if TBB
#define mk_task_group                 task_group tg;
#define wait_tasks                    tg.wait()
#define create_taskc(E)               tg.run(E)
#define create_taskc_if(x, E)         if(x) { create_taskc(E); } else { E(); }

#elif OPENMP
#define PRAGMA_OMP(x)                 _Pragma( #x )
#define mk_task_group
#define wait_tasks                    PRAGMA_OMP(omp taskwait)
#define create_taskc(E)               PRAGMA_OMP(omp task) E()
#define create_taskc_if(x, E)         if(x) { create_taskc(E); } else { E(); }

#elif MTHREAD || QTHREAD

#else
#define mk_task_group
#define wait_tasks
#define create_taskc(E)               E()
#define create_taskc_if(x, E)         E()

#endif

#endif
