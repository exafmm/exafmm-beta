#ifndef thread_h
#define thread_h

// Thread model
#if TBB
#include <tbb/task_group.h>
#include <tbb/task_scheduler_init.h>
using namespace tbb;

#elif MTHREAD
#include <mtbb/task_group.h>

#endif

/* pragma_omp macro */

#define DO_PRAGMA(x)                  _Pragma( #x )

#if _OPENMP
#define pragma_omp(x)              DO_PRAGMA(omp x)
#else
#define pragma_omp(x)              
#endif

// Task based threading macros.
// Usage example
// spawn_tasks {
//   spawn_task0(foo());
//   spawn_task1(x, x = bar());
//   sync_tasks;
//   spawn_task0_if(n > 10, baz(n));
//   spawn_task2(n > 10, a, x, x = bar(a));
//   sync_tasks;
// }
#if _OPENMP

#define spawn_tasks
#define sync_tasks                    pragma_omp(taskwait)

#define spawn_task0(E)                pragma_omp(task) E
#define spawn_task1(s0, E)            pragma_omp(task shared(s0)) E
#define spawn_task2(s0, s1, E)        pragma_omp(task shared(s0,s1)) E
#define spawn_task0_if(x, E)          if(x) { spawn_task0(E); } else { E; }
#define spawn_task1_if(x, s0, E)      if(x) { spawn_task1(s0,E); } else { E; }
#define spawn_task2_if(x, s0, s1, E)  if(x) { spawn_task2(s0,s1,E); } else { E; }

// tau: I think the following isn't quite right .
// I think 
// #pragma omp parallel {
// #pragma omp single {
//   }
// }
//  must be done only once in the toplevel
//
//#define DO_PRAGMA(x)                  _Pragma( #x )
//#define __init_tasks__                DO_PRAGMA(omp parallel) DO_PRAGMA(omp single)
//#define __sync_tasks__                DO_PRAGMA(omp taskwait)
//#define spawn_task0(E)                DO_PRAGMA(omp task) { E; }
//#define spawn_task1(s0, E)            DO_PRAGMA(omp task shared (s0)) { E; }  
//#define spawn_task2(s0, s1, E)        DO_PRAGMA(omp task shared(s0,s1)) { E; }
//#define spawn_task0_if(x, E)          if(x) { DO_PRAGMA(omp task) { E; } } else { E; }

#elif TBB || MTHREAD

#define spawn_tasks                   task_group tg;
#define sync_tasks                    tg.wait()

#define spawn_task0(E)                tg.run([=] { E; })
#define spawn_task1(s0, E)            tg.run([=,&s0] { E; })
#define spawn_task2(s0, s1, E)        tg.run([=,&s0,&s1] { E; })

#define spawn_task0_if(x, E)          if(x) { spawn_task0(E); } else { E; }
#define spawn_task1_if(x, s0, E)      if(x) { spawn_task1(s0,E); } else { E; }
#define spawn_task2_if(x, s0, s1, E)  if(x) { spawn_task2(s0,s1,E); } else { E; }

#else  /* not _OPENMP, TBB, or MTHREAD */

#define spawn_tasks
#define sync_tasks

#define spawn_task0(E)                E
#define spawn_task1(s0, E)            E
#define spawn_task2(s0, s1, E)        E

#define spawn_task0_if(x, E)          E
#define spawn_task1_if(x, s0, E)      E
#define spawn_task2_if(x, s0, s1, E)  E

#endif

#endif	/* thread_h */
