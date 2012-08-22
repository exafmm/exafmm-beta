#ifndef thread_h
#define thread_h
#if TBB
#include <tbb/task_group.h>
#define TASKS 1
using namespace tbb;

#elif MTHREAD
#include <mtbb/task_group.h>
#define TASKS 1

#endif

#if TASKS

#define __init_tasks__                task_group tg
#define __sync_tasks__                tg.wait()
#define spawn_task0(E)                tg.run([=] { E; })
#define spawn_task1(s0, E)            tg.run([=,&s0] { E; })
#define spawn_task2(s0, s1, E)        tg.run([=,&s0,&s1] { E; })
#define spawn_task0_if(x, E)          if (x) { tg.run([=] { E; }); } else { E; }

#else

#define __init_tasks__
#define __sync_tasks__
#define spawn_task0(E)                E
#define spawn_task1(s0, E)            E
#define spawn_task2(s0, s1, E)        E
#define spawn_task0_if(x, E)          E

#endif

#endif
