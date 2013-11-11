/* 
 * pocsort  (recursive counting sort)
 *
 * vector<int> v;
 * pocsort<int> a;
 * a.sort_vector(v, 20);
 * 
 */
#pragma once
#ifndef __pocsort_h__
#define __pocsort_h__

#include <assert.h>
#include <vector>
#include "thread.h"

#define create_task0(E)       __tg__.run([=] { E; })
#define create_task0_if(x, E) if (x) { create_task0(E); } else { E; }

/* 
recursive binsort with processor-oblivious parallel counting and shuffling

it is a hybrid of quicksort and bin sort. 

i don't know whether the algorithm has been published or given a
name.

assume we are given an array to sort A, each value of which has a key
= 0 and < M.  we assume we are also given a temporary array of the
same size as A. the algorithm works as follows.

(1) divide the key range into F equal bins (i.e.,
[0,M) = [0,M/F) + [M/F,2M/F) + ... + [(F-1)M/F,M))
      = B0      + B1         + ... + B_{F-1}

(2) count phase: count the number of elements that fall into each
bin, in parallel; let C_k be the number of elements that fall into
bin B_k.

(3) shuffle phase: as a result of (2), we know where values in each
bin B_k should go in the result array; i.e., values in bin B_0 should
go to [0,C0), values in bin B_1 to [C0,C0+C1), etc.  so we move
values to an appropriate region, again in parallel.

(4) after we are done, we recursively sort each sub-array in parallel

special case 1)

when F = 2, it's almost equivalent to quicksort, except that we do
not carefully choose pivot (although we could, with a slightly
complex partitioning code), and our method partitions values into a
separate array rather than in place.  our advantage and interest are,
of course, we can make F larger than 2 (e.g., 8 or 16).

special case 2)

when we make F = the number of workers, and do not recursively
parallelize step (4) beyond the first level (i.e., worker k sorts bin
B_k sequentially), then it is equivalent to the usual binsort.  In
that case, load balancing may not be good if key distributions are
skewed.  Our method overcomes this by recursively parallelizing sorting
each bin.

one difficulty associated with our method is how to parallelize
count/shuffle phase (step (2) and (3) above), in a
"processor-oblivious" manner, without knowing the number of workers.
our practical interest here is how to make parallel sort in
task-parallel languages as fast as parallel binsort/radix sort; a
theoretical interest here is binsort type method that has a critical
path length in o(n) when we have an extremely large (O(n^{1/2}))
processors.

to this end, the counting phase creates a tree describing how many
values fall into each bin from recursively divided subsections; the
root node describes how many values of the entire array fall into
each bin; the root node has F children, each representing the
consecutive 1/F of the original array, etc. such a tree can be
trivially built by parallel recursion; given this tree, we shuffle
values again using parallel recursion.

*/

const int masklen = 5;
const int F = 1 << masklen;/* number of recursions */


/* a vector of F long values */
struct long_vec {
  long c[F];

  /* v.prefix_sum returns a long_vec x
    s.t. 
  x.c[0] = begin;
    x.c[1] = begin + v.c[0]
    x.c[2] = begin + v.c[0] + v.c[1]
      ...
    x.c[F-1] = begin + v.c[0] + v.c[1] + ... + v.c[F-2]
 */
  long_vec prefix_sum() {
    long_vec s;
    long p = 0;
    for (long i = 0; i < F; i++) {
      s.c[i] = p;
      p += c[i];
    }
    return s;
  }

  long_vec operator+=(const long_vec& b) {
    for (int k = 0; k < F; k++) {
      this->c[k] += b.c[k];
    }
    return *this;
  }

  long_vec operator+(const long_vec& b) {
    long_vec c;
    for (int k = 0; k < F; k++) {
      c.c[k] = this->c[k] + b.c[k];
    }
    return c;
  }

};

/* a tree of long_vec */
struct long_vec_tree {
  long_vec counts;
  long_vec_tree * children[F];
};

template <class T>
class pocsort {
  typedef typename std::vector<T>::iterator TI;

private:
  /* return the maximum number of nodes necessary to count */
  long max_nodes_to_count(long n, long pcnt_th) {
    if (n <= pcnt_th) return 1;
    else return 1 + F * max_nodes_to_count((n + F - 1) / F, pcnt_th);
  }

  long max_nodes_to_sort(long n, long pcnt_th) {
    long d = (F - 1) * pcnt_th;
    return (n * F * F + d - 1) / d;
  }

  /* count how many values between [a_beg,a_end) fall into each bin 
    and put the result in root; if [a_beg,a_end) is large, 
    divide the region into subregions, count in parallel, and 
    aggregate the result; we keep the result for each divided
    subregions to later shuffle values in parallel. so the result
    is a tree each node describes how many elements of a particular
    region of the array fall into each bin */
  long_vec_tree * count(TI a_beg, TI a_end, int keylen,
			long_vec_tree * root,
			long_vec_tree * c_beg,
			__attribute__ ((unused)) long_vec_tree * c_end,
			long prec_th, long pcnt_th) {
    for (int k = 0; k < F; k++) {
      root->counts.c[k] = 0;
      root->children[k] = NULL;
    }
    if (a_end - a_beg <= pcnt_th) {
      /* small region; count elements sequentially and return
	 a singleton tree */
      int shift = (keylen >= masklen ? keylen - masklen : 0);
      for (TI p = a_beg; p < a_end; p++) {
	int k = (p->x >> shift) & ((1 << masklen) - 1);
	root->counts.c[k]++;
      }
    } else {
      /* divide the array, count each sub array, and make a tree */
      long_vec_tree * c = c_beg;
      task_group;
      for (int k = 0; k < F; k++) {
	TI b = a_beg + ((a_end - a_beg) * k) / F;
	TI e = a_beg + ((a_end - a_beg) * (k + 1)) / F;
	long n_nodes = max_nodes_to_count(e - b, pcnt_th);
	assert(c + n_nodes <= c_end);
	create_task0(root->children[k] 
		     = this->count(b, e, keylen, c, c + 1, 
				   c + n_nodes, prec_th, pcnt_th));
	c += n_nodes;
      }
      wait_tasks;
      for (int k = 0; k < F; k++) {
	root->counts += root->children[k]->counts;
      }
    }
    return root;
  }

  /* move values in [a_beg,a_end) into appropriate
    positions in [t_beg,...); specifically values
    in bin k should be put from t_beg[offsets[k] */
  void move(TI a_beg, TI a_end, TI t_beg,
	    long_vec offsets, long_vec_tree * r, int keylen) {
    if (r->children[0] == NULL) {
      /* leaf. shuffle values sequentially */
      int shift = (keylen >= masklen ? keylen - masklen : 0);
      for (TI p = a_beg; p < a_end; p++) {
	unsigned long k = (p->sort_key() >> shift) & ((1UL << masklen) - 1);
	t_beg[offsets.c[k]] = *p;
	offsets.c[k]++;
      }
    } else {
      /* otherwise recursively divide the region and 
	 shuffle values of each subregion in parallel */
#if 1
      task_group;
      for (int k = 0; k < F; k++) {
	assert(r->children[k]);
	TI b = a_beg + ((a_end - a_beg) * k) / F;
	TI e = a_beg + ((a_end - a_beg) * (k + 1)) / F;
	create_task0(this->move(b, e, t_beg, offsets, r->children[k], keylen));
	offsets += r->children[k]->counts;
      }
      wait_tasks;
#else
      mtbb::parallel_for(0, F, [=] (int k) {
	  assert(r->children[k]);
	  long_vec offs = offsets;
	  for (int kk = 1; kk < k; kk++) 
	    offs += r->children[kk]->counts;
	  TI b = a_beg + ((a_end - a_beg) * k) / F;
	  TI e = a_beg + ((a_end - a_beg) * (k + 1)) / F;
	  this->move(b, e, t_beg, offs, r->children[k], keylen);
	});
#endif
    }
  }

  inline long choose_min_idx(TI l, long begin, long end) {
    long m = begin;
    for (long i = begin + 1; i < end; i++) {
      if (l[i] < l[m]) {/* operator<  */
	m = i;
      }
    }
    return m;
  }

  inline void selection_sort_range(TI l, TI r, TI d) {
    long n = r - l;
    for (long i = 0; i < n; i++) {
      long j = choose_min_idx(l, i, n);
      T t = l[i];
      d[i] = l[j];
      l[j] = t;
    }
  }

  /* sort values in [a_beg,a_end) according to the lowermost keylen
    bits of each element.  depending on whethere dest == 0 or 1, the
    result goes either to [a_beg,...) (when dest == 0) or into
    [t_beg,...)  (when dest == 1).


  */
  void count_and_move(TI a_beg, TI a_end, TI t_beg, int dest,
		      int keylen, 
		      long_vec_tree * v_beg, long_vec_tree * v_end, 
		      long rec_th, long prec_th, long pcnt_th) {
    assert(keylen > 0);
    assert(max_nodes_to_sort(a_end - a_beg, pcnt_th) <= v_end - v_beg + 1);
    if (a_end - a_beg <= rec_th) {
      this->selection_sort_range(a_beg, a_end, (dest == 0 ? a_beg : t_beg));
    } else {
      long_vec_tree root[1]; 
      long_vec_tree * r = count(a_beg, a_end, keylen, root, v_beg, v_end, prec_th, pcnt_th);
      long_vec counts = r->counts;
      long_vec offsets = counts.prefix_sum();
      move(a_beg, a_end, t_beg, offsets, r, keylen);
      if (keylen - masklen > 0) {
	long_vec_tree * v = v_beg;
	task_group;
	for (int k = 0; k < F; k++) {
	  long n_nodes = max_nodes_to_sort(counts.c[k], pcnt_th);
	  assert(v + n_nodes - 1 <= v_end);
	  create_task0_if(a_end - a_beg > prec_th,
			  this->count_and_move(t_beg + offsets.c[k], 
					       t_beg + offsets.c[k] + counts.c[k], 
					       a_beg + offsets.c[k], 
					       1 - dest, keylen - masklen, 
					       v, v + n_nodes - 1,
					       rec_th, prec_th, pcnt_th));
	  v += n_nodes - 1;
	}
	wait_tasks;
      }
    }
  }

public:

  /* sort a into s */
  void sort_vector(std::vector<T>& a, std::vector<T>& s, int key_size) {
    long rec_th = 10;/* use this algorithm until n <= this value */
    long prec_th = 100;/* use parallel recursion until n <= this value */
    long pcnt_th = 1000;/* use parallel counting until n <= this value */
    long n_nodes = max_nodes_to_sort(a.size(), pcnt_th);

    long_vec_tree * v_beg = new long_vec_tree[n_nodes];
    long_vec_tree * v_end = v_beg + n_nodes;
    count_and_move(a.begin(), a.end(), s.begin(), 0, key_size, v_beg, v_end,
		   rec_th, prec_th, pcnt_th);
    delete[] v_beg;
  }
};

#endif
