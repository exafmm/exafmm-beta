#ifndef thread_h
#define thread_h

#if TBB

/* if you want to use TBB and do not want to
   install MassiveThreads */

#include <tbb/task_group.h>
#define mk_task_group        tbb::task_group __tg__
#define wait_tasks           __tg__.wait()
#define create_taskc(E)      __tg__.run(E)
#define create_taskc_if(x, E) if (x) { create_taskc(E); } else { E(); }

#else

/* this works for TBB, MassiveThreads, Nanos++, Qthreads,
   and OpenMP, provided you have installed MassiveThreads.
   switch between them by giving one of 
   -DTO_TBB, -DTO_MTHREAD, -DTO_QTHREAD, -DTO_NANOX, or
   -DTO_OMP */

#include <mtbb/task_group.h>
#include <tpswitch/tpswitch.h>

#endif


#if 0

// You can erase everything below


#if COMMON_CILKH
#include <common.cilkh>

#else  /* ! COMMON_CILKH */

#if TBB
#include <tbb/task_group.h>
typedef tbb::task_group task_group_t;
#elif MTHREAD || QTHREAD || NANOX
#include <mtbb/task_group.h>
typedef mtbb::task_group task_group_t;
#endif

#if TBB || MTHREAD || QTHREAD || NANOX
#define task_group           task_group_t __tg__
#define wait_tasks           __tg__.wait()
#define create_task(E)       __tg__.run(E)
#define create_task_if(x, E) if (x) { create_task(E); } else { E(); }
#else /* not TBB || MTHREAD || QTHREAD || NANOX */
#define task_group
#define wait_tasks
#define create_task(E)       E();
#define create_task_if(x, E) E();
#endif /* TBB || MTHREAD || QTHREAD || NANOX */

#endif /* COMMON_CLIKH */

#endif	/* the end of #if 0 */

#endif
