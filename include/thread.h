#ifndef thread_h
#define thread_h

#if COMMON_CILKH
#include <common.cilkh>
#else

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
