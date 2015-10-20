#ifndef thread_h
#define thread_h

#if TBB
#include <tbb/task_group.h>
using namespace tbb;
#elif MTHREAD
#include <mtbb/task_group.h>
#endif

#if TBB || MTHREAD
#define spawn_tasks                   task_group tg;
#define sync_tasks                    tg.wait()
#define spawn_task0(E)                tg.run([=] { E; })
#define spawn_task1(s0, E)            tg.run([=,&s0] { E; })
#define spawn_task2(s0, s1, E)        tg.run([=,&s0,&s1] { E; })
#define spawn_task0_if(x, E)          if(x) { spawn_task0(E); } else { E; }
#define spawn_task1_if(x, s0, E)      if(x) { spawn_task1(s0,E); } else { E; }
#define spawn_task2_if(x, s0, s1, E)  if(x) { spawn_task2(s0,s1,E); } else { E; }

#else  /* not TBB, or MTHREAD */
#define spawn_tasks
#define sync_tasks
#define spawn_task0(E)                E
#define spawn_task1(s0, E)            E
#define spawn_task2(s0, s1, E)        E
#define spawn_task0_if(x, E)          E
#define spawn_task1_if(x, s0, E)      E
#define spawn_task2_if(x, s0, s1, E)  E

#endif

#endif
