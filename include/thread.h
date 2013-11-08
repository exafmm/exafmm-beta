#ifndef thread_h
#define thread_h

#if COMMON_CILKH
#include <common.cilkh>
#else

#if _OPENMP
#elif TBB
#include <tbb/task_group.h>
typedef tbb::task_group task_group_t;
#elif MTHREAD || QTHREAD || NANOX
#include <mtbb/task_group.h>
typedef mtbb::task_group task_group_t;
#endif

#if _OPENMP
#define PRAGMA_OMP(x)                 _Pragma( #x )
#define task_group
#define wait_tasks                    PRAGMA_OMP(omp taskwait)
#if CXX_LAMBDA
#define create_task0(E)               PRAGMA_OMP(omp task) E;
#define create_task1(s0, E)           PRAGMA_OMP(omp task shared(s0)) E;
#define create_task2(s0, s1, E)       PRAGMA_OMP(omp task shared(s0,s1)) E;
#define create_taskA(E)               PRAGMA_OMP(omp task default(shared)) E;
#define create_task0_if(x, E)         PRAGMA_OMP(omp task if (x)) E;
#define create_task1_if(x, s0, E)     PRAGMA_OMP(omp task shared(s0) if (x)) E;
#define create_task2_if(x, s0, s1, E) PRAGMA_OMP(omp task shared(s0,s1) if (x)) E;
#define create_taskA_if(x, E)         PRAGMA_OMP(omp task default(shared) if (x)) E;
#else
#define create_taskc(E)               PRAGMA_OMP(omp task) E();
#define create_taskc_if(x, E)         PRAGMA_OMP(omp task if (x)) E();
#endif

#elif TBB || MTHREAD || QTHREAD || NANOX
#define task_group                    task_group_t __tg__
#define wait_tasks                    __tg__.wait()
#if CXX_LAMBDA
#define create_task0(E)               __tg__.run([=] { E; })
#define create_task1(s0, E)           __tg__.run([=,&s0] { E; })
#define create_task2(s0, s1, E)       __tg__.run([=,&s0,&s1] { E; })
#define create_taskA(E)               __tg__.run([&] { E; })
#define create_task0_if(x, E)         if (x) { create_task0(E); } else { E; }
#define create_task1_if(x, s0, E)     if (x) { create_task1(s0, E); } else { E; }
#define create_task2_if(x, s0, s1, E) if (x) { create_task2(s0, s1, E); } else { E; }
#define create_taskA_if(x, E)         if (x) { create_taskA(E); } else { E; }
#else
#define create_taskc(E)               __tg__.run(E)
#define create_taskc_if(x, E)         if (x) { create_taskc(E); } else { E(); }
#endif

#else  /* not _OPENMP, TBB, MTHREAD, QTHREAD, or NANOX */
#define task_group
#define wait_tasks
#if CXX_LAMBDA
#define create_task0(E)                E;
#define create_task1(s0, E)            E;
#define create_task2(s0, s1, E)        E;
#define create_taskA(E)                E;
#define create_task0_if(x, E)          E;
#define create_task1_if(x, s0, E)      E;
#define create_task2_if(x, s0, s1, E)  E;
#define create_taskA_if(x, E)          E;
#else
#define create_taskc(E)                E();
#define create_taskc_if(x, E)          E();
#endif
#endif
#endif	/* COMMON_CLIKH */

#endif
