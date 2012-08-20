/* 
 * task_parallel.h
 */

#pragma once

#ifndef task_parallel_h
#define task_parallel_h

#if TBB
/* on TBB, task_group is already in */
#include <tbb/task_group.h>
#define task_group_defined 1
using namespace tbb;

#elif MTHREAD
/* on MassiveThreads, it's already in the library too */
#include <mtbb/task_group.h>
#define task_group_defined 1

#elif PTHREAD || QTHREAD || NANOX
/* roll our own in the way similar to MassiveThreads on
   Pthreads/Qthread/Nano++ */

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <functional>

#if PTHREAD
#define th_func_ret_type void *
#endif
#if QTHREAD
#define th_func_ret_type aligned_t
#endif
#if NANOX
#define th_func_ret_type void
#endif

struct task {
  pthread_t tid;
  std::function<void ()> f;
#if QTHREAD
  aligned_t ret;
#endif
};

th_func_ret_type invoke_task(void * arg_) {
  task * arg = (task *)arg_;
  std::function<void()> f = arg->f;
  f();
#if PTHREAD || QTHREAD
  return 0;
#endif
}

#if NANOX
nanos_smp_args_t invoke_task_arg={invoke_task};
#endif

#if !defined(TASK_GROUP_INIT_SZ)
#define TASK_GROUP_INIT_SZ 10
#endif

#if !defined(TASK_GROUP_NULL_CREATE)
#define TASK_GROUP_NULL_CREATE 0
#endif

struct task_list_node {
  task_list_node * next;
  int capacity;
  int n;
  task a[TASK_GROUP_INIT_SZ];
};

struct task_group {
  task_list_node first_chunk_[1];
  task_list_node * head;
  task_list_node * tail;
  task_group() {
    head = first_chunk_;
    tail = first_chunk_;
    head->next = NULL;
    head->capacity = TASK_GROUP_INIT_SZ;
    head->n = 0;
  }
  ~task_group() {
    task_list_node * q = NULL;
    for (task_list_node * p = head; p; p = q) {
      q = p->next;
      if (p != first_chunk_) delete p;
    }
  }

  void extend() {
    task_list_node * new_node = new task_list_node();
    new_node->next = NULL;
    new_node->n = 0;
    new_node->capacity = TASK_GROUP_INIT_SZ;
    tail->next = new_node;
    tail = new_node;
  }
  void run(std::function<void ()> f) {
    if (tail->n == tail->capacity) {
      if (tail->next == NULL) extend();
      else tail = tail->next;
      assert(tail->n == 0);
    }
    task * t = &tail->a[tail->n];
    t->f = f;
    if (TASK_GROUP_NULL_CREATE) {
      invoke_task((void *)t);
    } else {
      tail->n++;
#if PTHREAD
      pthread_create(&t->tid, NULL, invoke_task, (void*)t);
#elif QTHREAD
      qthread_fork(invoke_task, (void*)t, &t->ret);
#elif NANOX
      nanos_wd_t wd=NULL;
      nanos_device_t dev[1] = {NANOS_SMP_DESC(invoke_task_arg)};
      nanos_wd_props_t props;	// originally, ={true,false,false};
      props.mandatory_creation = true;
      props.tied = false;
      props.reserved0 = false;
      NANOS_SAFE(nanos_create_wd(&wd,1,dev,sizeof(struct task),
				 __alignof__(struct task),
				 (void**)&t,nanos_current_wd(),&props,0,NULL));
      NANOS_SAFE(nanos_submit(wd,0,0,0));
#else
#error "neither PTHREAD, QTHREAD, nor NANOX defined"
#endif
    }
  }

  void wait() {
    int n_joined = 0;
    for (task_list_node * p = head; p && p->n; p = p->next) {
      for (int i = 0; i < p->n; i++) {
#if TASK_GROUP_NULL_CREATE 
	/* noop */
#elif PTHREAD
	void * ret;
	pthread_join(p->a[i].tid, &ret);
#elif QTHREAD
	aligned_t ret;
	qthread_readFF(&ret,&p->a[i].ret);
#elif NANOX
	if (n_joined == 0) {
	  NANOS_SAFE(nanos_wg_wait_completion(nanos_current_wd()));
	}
#else
#error "neither PTHREAD, QTHREAD, nor NANOX defined"
#endif
	n_joined++;
      }
      p->n = 0;
    }
  }

};
#define task_group_defined 1

#else	/* ! TBB || MTHREAD || PTHREAD || QTHREAD || NANOX || ... */

#define task_group_defined 0

#endif


/* 
 * macros that absorb syntactic difference between
 * Cilk's spawn and sync and task_group's run and wait.
 * they are also made syntactically valid on other platforms.
 * 
 * An example:
 *
 * {
 *    int a[3];
 *    int b[1000];
 *    __spawn_tasks__;  // declare a block that creates tasks
 *    spawn_task0(spawn f(x));
 *    spawn_task1(a, a[0] = spawn f(x));
 *    spawn_task2(a, b, a[1] = spawn g(x, b));
 *    __sync__;
 * }
 *
 * CAUTION(1): with Cilk, you need 'spawn' keyword inside the macro.
 * spawn_task0(f(x))  just calls f(x) serially.  This is to unify
 * spawning void functions and spawning value-returning functions
 * into a single macro;
 * 
 * With Cilk, the above gets expanded to:
 *
 *    int a[3];
 *    int b[1000];
 *    
 *    spawn f(x);
 *    a[0] = spawn f(x);
 *    a[1] = spawn g(x, b);
 *    sync;
 *
 * With TBB and other platforms supporting task_group, the above gets
 * expanded to:
 *    int a[3];
 *    int b[1000];
 *    task_group tg;
 *    tg.run([=] { f(x) });
 *    tg.run([=,&a] { a[0] = f(x); });
 *    tg.run([=,&a,&b] { a[1] = spawn g(x, b); });
 *    tg.wait();

 * spawn_task{0,1,2} are expanded using tg.run.  
 * spawn_task{0,1,2} takes {0,1,2} variables
 * that are shared between the parent and child, respectively.
 *
 * CAUTION(2): for value-returning spawn, you typically need to specify
 * at least one shared variable---the variable via which you receive the
 * result.  For example, the code below is almost certainly not what you want:
 *    spawn_task0(y = spawn f(x));
 * as it does not share the variable y between the parent and child.
 * you probably want to do the following instead.
 *    spawn_task1(y, y = spawn f(x));
 * To accommodate cases in which the lefthand side of the assignment
 * is not a variable (e.g., a[3], a.f.g, a->f), we must owe the programmer
 * to specify what is the variable name you actually want to share. e.g.,
 *    int a[3]; 
 *    spawn_task1(a, a[1] = spawn f(x));
 * 
 */

#if task_group_defined

#define __spawn_tasks__               task_group tg
#define __sync__                      tg.wait()

#define spawn_task0(E)                tg.run([=] { E; })
#define spawn_task1(s0, E)            tg.run([=,&s0] { E; })
#define spawn_task2(s0, s1, E)        tg.run([=,&s0,&s1] { E; })

#define spawn_task0_if(x, E)          if (x) { tg.run([=] { E; }); } else { E; }
#define spawn_task1_if(x, s0, E)      if (x) { tg.run([=,&s0] { E; }); } else { E; }
#define spawn_task2_if(x, s0, s1, E)  if (x) { tg.run([=,&s0,&s1] { E; }); } else { E; }

#define call_task(E)                  E
#define spawn 

#else  /* task_group_defined */

#define __spawn_tasks__
#define __sync__

#define spawn_task0(E)                E
#define spawn_task1(s0, E)            E
#define spawn_task2(s0, s1, E)        E

#define spawn_task0_if(x, E)          E
#define spawn_task1_if(x, s0, E)      E
#define spawn_task2_if(x, s0, s1, E)  E

#define call_task(E)                  E
#define spawn 

#endif	/* task_group_defined */

#endif
