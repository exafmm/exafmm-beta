#ifndef thread_h
#define thread_h

// Thread model
#if TBB
#include <tbb/task_group.h>
using namespace tbb;
#define __init_tasks__                task_group tg
#define __sync_tasks__                tg.wait()
#define spawn_task0(E)                tg.run([=] { E; })
#define spawn_task1(s0, E)            tg.run([=,&s0] { E; })
#define spawn_task2(s0, s1, E)        tg.run([=,&s0,&s1] { E; })
#define spawn_task0_if(x, E)          if (x) { tg.run([=] { E; }); } else { E; }

#elif MTHREAD
#include <mtbb/task_group.h>
#define __init_tasks__                task_group tg
#define __sync_tasks__                tg.wait()
#define spawn_task0(E)                tg.run([=] { E; })
#define spawn_task1(s0, E)            tg.run([=,&s0] { E; })
#define spawn_task2(s0, s1, E)        tg.run([=,&s0,&s1] { E; })
#define spawn_task0_if(x, E)          if (x) { tg.run([=] { E; }); } else { E; }

#elif OPENMP
#define DO_PRAGMA(x)                  _Pragma( #x )
#define __init_tasks__                DO_PRAGMA(omp parallel) DO_PRAGMA(omp single)
#define __sync_tasks__                DO_PRAGMA(omp taskwait)
#define spawn_task0(E)                DO_PRAGMA(omp task) { E; }
#define spawn_task1(s0, E)            DO_PRAGMA(omp task shared (s0)) { E; }  
#define spawn_task2(s0, s1, E)        DO_PRAGMA(omp task shared(s0,s1)) { E; }
#define spawn_task0_if(x, E)          if(x) { DO_PRAGMA(omp task) { E; } } else { E; }

#else
#define __init_tasks__
#define __sync_tasks__
#define spawn_task0(E)                E
#define spawn_task1(s0, E)            E
#define spawn_task2(s0, s1, E)        E
#define spawn_task0_if(x, E)          E

#endif

#endif
