#ifndef thread_h
#define thread_h

/* 

   With new MassiveThreads, we can just include
   <mtbb/task_group.h> instead of
   <tbb/task_group.h> and you can use the same
   task_group class for BOTH TBB and
   MassiveThreads (and Qthreads and Nanos++, for
   that matter).  Switch them by defining one of
   TO_TBB, TO_MTHREAD, TO_QTHREAD, and TO_NANOX in
   the command line.

   tpswitch/tpswitch.h provides create_taskc, 
   create_taskc_if, etc.
 */


#if 1
// These two lines are all we need
#include <mtbb/task_group.h>
#include <tpswitch/tpswitch.h>


#else

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

#endif

#endif
