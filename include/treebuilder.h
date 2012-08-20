#ifndef topdown_h
#define topdown_h
#include "evaluator.h"

#if PARALLEL_EVERYTHING

/* 8 elements counter to record the number of bodies
   that should go to each child */
typedef vec<8,int> ivec8;

/* tree of 8 elements counters.
   a node corresponds to a contiguous region of 
   the bodies vector. */
struct ivec8Tree {
  ivec8 counts;
  ivec8Tree * children[2];
};

struct XNode {
  bool  NOCHILD;                                                //!< Flag for twig nodes
  int   LEVEL;                                                  //!< Level in the tree structure
  int BODY;
  int   NLEAF;                                                  //!< Number of descendant leafs
  XNode * CHILD[8];                                             //!< Pointer to child node
  vec3  X;                                                      //!< Coordinate at center
  int NNODE;
};

#endif

class TreeBuilder : public Evaluator {
private:
  Nodes    nodes;
  Leafs    leafs;
  B_iter   BN;
  C_iter   CN;
  unsigned NLEAF;
  unsigned NCELL;
  int      MAXLEVEL;

protected:
  real_t localRadius;                                           //!< Radius of local root cell
  vec3 localCenter;                                             //!< Center of local root cell
  vec3 localXmin;                                               //!< Local Xmin for a given rank
  vec3 localXmax;                                               //!< Local Xmax for a given rank

private:
  void init(Node &node) {
    node.NOCHILD = true;
    node.NLEAF = 0;
    node.LEAF = NULL;
    for (int b=0; b<8; b++) node.CHILD[b] = -1;
  }

  inline long long getIndex(vec3 X, int level) {
    float d = 2 * localRadius / (1 << level);
    long long ix = int((X[0] - localXmin[0]) / d);
    long long iy = int((X[1] - localXmin[1]) / d);
    long long iz = int((X[2] - localXmin[2]) / d);
    long long id = ((1 << 3 * level) - 1) / 7;
    assert( level < 22 );
    for (int l=0; l<level; l++) {
      id += (ix & 1) << (3 * l);
      id += (iy & 1) << (3 * l + 1);
      id += (iz & 1) << (3 * l + 2);
      ix >>= 1;
      iy >>= 1;
      iz >>= 1;
    }
    return id;
  }

  inline void addChild(int octant, int n) {
    N_iter N = nodes.begin()+n;
    Node child;
    init(child);
    child.LEVEL = N->LEVEL+1;
    child.X = N->X;
    real_t r = localRadius / (1 << child.LEVEL);
    for (int d=0; d<3; d++) {                                   // Loop over dimensions
      child.X[d] += r * (((octant & 1 << d) >> d) * 2 - 1);     //  Calculate new center position
    }                                                           // End loop over dimensions
    N->NOCHILD = false;
    N->CHILD[octant] = nodes.size();
    nodes.push_back(child);
    NCELL++;
  }

  void splitNode(int n) {
    while (nodes[n].NLEAF > NCRIT) {
      int c = 0;
      Leaf *Ln;
      for (Leaf *L=nodes[n].LEAF; L; L=Ln) {
        N_iter N = nodes.begin()+n;
        Ln = L->NEXT;
        int octant = (L->X[0] > N->X[0]) + ((L->X[1] > N->X[1]) << 1) + ((L->X[2] > N->X[2]) << 2);
        if (N->CHILD[octant] == -1) {
          addChild(octant,n);
        }
        c = nodes[n].CHILD[octant];
        Node *child = &nodes[c];
        L->NEXT = child->LEAF;
        child->LEAF = &*L;
        child->NLEAF++;
      }
      n = c;
    }
  }

  void nodes2cells(int i, C_iter C) {
    C->R      = localRadius / (1 << nodes[i].LEVEL);
    C->X      = nodes[i].X;
    C->NDLEAF = nodes[i].NLEAF;
    C->LEAF   = BN;
    C->ICELL  = getIndex(C->X,nodes[i].LEVEL);
    if (nodes[i].NOCHILD) {
      C->CHILD = 0;
      C->NCHILD = 0;
      C->NCLEAF = nodes[i].NLEAF;
      for (Leaf *L=nodes[i].LEAF; L; L=L->NEXT) {
        BN->IBODY = L->I;
        BN++;
      }
    } else {
      C->NCLEAF = 0;
      int nsub=0;
      for (int octant=0; octant<8; octant++) {
        if (nodes[i].CHILD[octant] != -1) {
          ++nsub;
        }
      }
      C_iter Ci = CN;
      C->CHILD = Ci - Ci0;
      C->NCHILD = nsub;
      CN += nsub;
      for (int octant=0; octant<8; octant++) {
        if (nodes[i].CHILD[octant] != -1) {
          Ci->PARENT = C - Ci0;
          nodes2cells(nodes[i].CHILD[octant], Ci++);
        }
      }
    }
  }

  void permuteBodies(Bodies &bodies) {
    Bodies buffer = bodies;
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
      B->ICELL = buffer[B->IBODY].ICELL;
      B->X     = buffer[B->IBODY].X;
      B->SRC   = buffer[B->IBODY].SRC;
      B->TRG   = buffer[B->IBODY].TRG;
    }
  }

protected:
  void setLeafs(Bodies &bodies) {
    startTimer("Set leafs");
    NLEAF = bodies.size();
    leafs.reserve(NLEAF);
    for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
      Leaf leaf;
      leaf.I = B-bodies.begin();
      leaf.X = B->X;
      leafs.push_back(leaf);
    }
    stopTimer("Set leafs",printNow);
  }

  void growTree() {
    startTimer("Grow tree");
    NCELL = 1;
    nodes.reserve(NLEAF);
    Node node;
    init(node);
    node.LEVEL = 0;
    node.X = localCenter;
    nodes.push_back(node);
    MAXLEVEL = 0;
    for (L_iter L=leafs.begin(); L!=leafs.end(); L++) {
      int n = 0;
      while (!nodes[n].NOCHILD) {
        N_iter N = nodes.begin() + n;
        int octant = (L->X[0] > N->X[0]) + ((L->X[1] > N->X[1]) << 1) + ((L->X[2] > N->X[2]) << 2);
        N->NLEAF++;
        if (nodes[n].CHILD[octant] == -1) addChild(octant,n);
        n = nodes[n].CHILD[octant];
      }
      L->NEXT = nodes[n].LEAF;
      nodes[n].LEAF = &*L;
      nodes[n].NLEAF++;
      if (nodes[n].NLEAF > NCRIT) splitNode(n);
      if (MAXLEVEL < nodes[n].LEVEL) MAXLEVEL = nodes[n].LEVEL;
    }
    MAXLEVEL++;
    stopTimer("Grow tree",printNow);
  }
  
  void linkTree(Bodies &bodies, Cells &cells) {
    startTimer("Link tree");
    cells.clear();
    cells.resize(NCELL);
    Ci0 = cells.begin();
    BN = bodies.begin();
    CN = Ci0 + 1;
    nodes2cells(0,Ci0);
    nodes.clear();
    leafs.clear();
    permuteBodies(bodies);
    stopTimer("Link tree",printNow);
  }

#if PARALLEL_EVERYTHING
 private:
  ivec8 prefixSum(ivec8 a, int beg) {
    ivec8 s;
    int p = beg;
    for (int i=0; i<8; i++) {
      s[i] = p;
      p += a[i];
    }
    return s;
  }

  /* make XNode instance that covers to bodies[beg] ... bodies[end] */
  XNode * makeNode(int level, int beg, int end, vec3 X, bool nochild) {
    XNode * n = new XNode();
    n->LEVEL = level; 
    n->BODY = beg; 
    n->NLEAF = end - beg; 
    n->X = X;
    n->NNODE = 1;
    n->NOCHILD = nochild;
    if (nochild) {
      for (int k=0; k<8; k++) n->CHILD[k] = NULL;
    }
    return n;
  }

  int max_ivec8_nodes_to_count(int n, int leaf_len) {
    if (n <= leaf_len) return 1;
    else return 4 * ((n - 1) / leaf_len) - 1;
  }

  /* maximum ivec8 nodes for n bodies */
  int max_ivec8_nodes_to_build(int n, int leaf_len) {
    return (4 * n) / leaf_len;
  }

  /* given 
       bodies[beg:end] : a set of bodies
       X               : the geometric center of a cubic region,
     count how many bodies in bodies[beg:end] as well as its subsections
     are in each of the eight child cubes.
     the result is returned as a tree of counters.  its root node describes
     how many bodies in bodies[beg:end] are in ach of the eight child cubes,
     each of its eight children describes how many bodies in each of its
     1/8 sections, and so on. the section is divided until the number of
     bodies <= leaf_len.
     the root node is allocated at t_root, its descendents allocated between
     t_beg and t_end.
*/
  ivec8Tree * countBodies(Bodies& bodies, int beg, int end, vec3 X, 
                          ivec8Tree * t_root, ivec8Tree * t_beg, ivec8Tree * t_end, 
                          int leaf_len) {
    assert(max_ivec8_nodes_to_count(end - beg, leaf_len) <= t_end - t_beg + 1);
    if (end - beg <= leaf_len) {
      /* the section is small enough -> count sequentially */
      for (int k=0; k<8; k++) t_root->counts[k] = 0;
      t_root->children[0] = t_root->children[1] = NULL;
      for (int i=beg; i<end; i++) {
        vec3 x = bodies[i].X;
        int oct = (x[0] > X[0]) + ((x[1] > X[1]) << 1) + ((x[2] > X[2]) << 2);
        t_root->counts[oct]++;
      } 
    } else {
      /* divide the section into two subsections and count them
         in parallel */
      int mid = (beg + end) / 2;
      int n0 = max_ivec8_nodes_to_count(mid - beg, leaf_len);
      int n1 = max_ivec8_nodes_to_count(end - mid, leaf_len);
      ivec8Tree * t_mid = t_beg + n0;
      assert(t_end - t_beg >= n0 + n1);
      __spawn_tasks__;
#if _OPENMP
#pragma omp task shared(bodies)
#endif
      spawn_task1(bodies,
                  spawn t_root->children[0] 
                  = countBodies(bodies, beg, mid, X, t_beg, t_beg + 1, t_beg + n0, leaf_len));
      call_task(spawn t_root->children[1] 
                = countBodies(bodies, mid, end, X, t_mid, t_mid + 1, t_mid + n1, leaf_len));
#if _OPENMP
#pragma omp taskwait
#endif
      __sync__;
      t_root->counts = t_root->children[0]->counts + t_root->children[1]->counts;
    }
    return t_root;
  }

  /* move bodies in bodies[beg:end] into t_bodies[beg:end], so each
     particle will be in the right subcube. positions are described
     in offsets. */
  void moveBodies(Bodies& bodies, Bodies& t_bodies, int beg, int end, 
                  ivec8Tree * t, ivec8 offsets, vec3 X) {
    if (t->children[0] == NULL) {
      /* it is leaf, so we move sequentially */
      for (int i=beg; i<end; i++) {
        vec3 x = bodies[i].X;
        int oct = (x[0] > X[0]) + ((x[1] > X[1]) << 1) + ((x[2] > X[2]) << 2);
        t_bodies[offsets[oct]] = bodies[i];
        offsets[oct]++;
      }
    } else {
      /* divide the section into two subsections
         and work on each in parallel */
      int mid = (beg + end) / 2;
      ivec8 offsets_mid = offsets + t->children[0]->counts;
      __spawn_tasks__;
#if _OPENMP
#pragma omp task shared(bodies, t_bodies)
#endif
      spawn_task2(bodies, t_bodies,
                  spawn moveBodies(bodies, t_bodies, beg, mid, 
                                   t->children[0], offsets, X));
      call_task(spawn moveBodies(bodies, t_bodies, mid, end,
                                 t->children[1], offsets_mid, X));
    
#if _OPENMP
#pragma omp taskwait
#endif
      __sync__;
    }
  }

  XNode * buildNodes(Bodies& bodies, Bodies& t_bodies, int dest,
                     int beg,  int end, 
                     ivec8Tree * t_beg, ivec8Tree * t_end, 
                     vec3 X, int level, int leaf_len) {
    assert(max_ivec8_nodes_to_build(end - beg, leaf_len) <= t_end - t_beg);
    if (beg == end) return NULL;
    if (end - beg <= NCRIT) {
      if (dest)
        for (int i=beg; i<end; i++) t_bodies[i] = bodies[i];
      return makeNode(level, beg, end, X, true);
    }
    XNode * node = makeNode(level, beg, end, X, false);
    ivec8Tree t_root_[1]; 
    ivec8Tree * t_root = countBodies(bodies, beg, end, X, t_root_, t_beg, t_end, leaf_len);
    ivec8 offsets = prefixSum(t_root->counts, beg);
    moveBodies(bodies, t_bodies, beg, end, t_root, offsets, X);
    ivec8Tree * t = t_beg;
    __spawn_tasks__;
    for (int k=0; k<8; k++) {
      int n_nodes = max_ivec8_nodes_to_build(t_root->counts[k], leaf_len);
      assert(t + n_nodes <= t_end);
#if _OPENMP
#pragma omp task shared(bodies, t_bodies)
#endif
      spawn_task2(t_bodies, bodies, {
          vec3 Y = X;
          real_t r = localRadius / (1 << (level + 1));
          for (int d=0; d<3; d++) {
            Y[d] += r * (((k & 1 << d) >> d) * 2 - 1);
          }
          node->CHILD[k] 
            = buildNodes(t_bodies, bodies, 1 - dest,
                         offsets[k], offsets[k] + t_root->counts[k],
                         t, t + n_nodes,
                         Y, level + 1, leaf_len);
        });
      t += n_nodes;
    }
#if _OPENMP
#pragma omp taskwait
#endif
    __sync__;
    for (int k=0; k<8; k++) {
      if (node->CHILD[k]) {
        node->NNODE += node->CHILD[k]->NNODE;
      }
    }
    return node;
  }

  int nodes2cellsRec(XNode * node, int parent, C_iter C, C_iter H, Cells& cells,
                     B_iter B0, int level) {
    assert(C < cells.end());
    C->PARENT = parent;
    C->R      = localRadius / (1 << node->LEVEL);
    C->X      = node->X;
    C->NDLEAF = node->NLEAF;
    C->LEAF   = B0 + node->BODY;
    if (node->NOCHILD) {
      C->CHILD = 0;
      C->NCHILD = 0;
      C->NCLEAF = node->NLEAF;
      assert(C->NCLEAF > 0);
      return node->LEVEL;
    } else {
      C->NCLEAF = 0;
      int nsub=0;
      char child_octants[8];
      for (int octant=0; octant<8; octant++) {
        if (node->CHILD[octant]) {
          child_octants[nsub] = octant;
          ++nsub;
        }
      }
      C_iter Ci = H;
      C->CHILD = Ci - Ci0;
      C->NCHILD = nsub;
      assert(C->NCHILD > 0);
      H += nsub;
      __spawn_tasks__;
      int levels[8];
      for (int k=0; k<nsub; k++) {
        int octant = child_octants[k];
#if _OPENMP
#pragma omp task if(node->NNODE > 1000) shared(levels, cells)
#endif
        spawn_task2_if(node->NNODE > 1000,
                       levels, cells,
                       spawn levels[k] 
                       = nodes2cellsRec(node->CHILD[octant], C - Ci0, Ci, 
                                        H, cells, B0, level + 1));
        Ci += 1;
        H += node->CHILD[octant]->NNODE - 1;
      }
#if _OPENMP
#pragma omp taskwait
#endif
      __sync__;
      int max_level = node->LEVEL;
      for (int k=0; k<nsub; k++) {
        int octant = child_octants[k];
        max_level = std::max(max_level, levels[k]);
        delete node->CHILD[octant];
      }
      return max_level;
    }
  }
  XNode * root_node;

 protected:
  void growTreeRec(Bodies &bodies, Bodies &t_bodies) {
    startTimer("Grow tree");
    int leaf_len = 1000;
    int n_nodes = max_ivec8_nodes_to_build(bodies.size(), leaf_len);
    ivec8Tree * t_beg = new ivec8Tree[n_nodes];
    ivec8Tree * t_end = t_beg + n_nodes;
    root_node = buildNodes(bodies, t_bodies, 0, 0, bodies.size(), 
                           t_beg, t_end, localCenter, 0, leaf_len);
    delete[] t_beg;
    stopTimer("Grow tree",printNow);
  }

  void linkTreeRec(Bodies& bodies, Cells &cells) {
    startTimer("Link tree");
    cells.resize(root_node->NNODE);
    Ci0 = cells.begin();
    MAXLEVEL = nodes2cellsRec(root_node, 0, Ci0, Ci0 + 1, cells, bodies.begin(), 0);
    delete root_node; root_node = NULL;
    stopTimer("Link tree",printNow);
  }

#endif


  void printTreeData(Cells &cells) {
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "Bodies               : " << cells.front().NDLEAF << std::endl;
    std::cout << "Cells                : " << cells.size()         << std::endl;
    std::cout << "Tree depth           : " << MAXLEVEL             << std::endl;
#if COUNT
    std::cout << "P2P calls            : " << NP2P                 << std::endl;
    std::cout << "M2P calls            : " << NM2P                 << std::endl;
    std::cout << "M2L calls            : " << NM2L                 << std::endl;
    std::cout << "traverse calls       : " << NTRAVERSE            << std::endl;
#endif
    std::cout << "-----------------------------------------------" << std::endl;
  }
};

#endif
