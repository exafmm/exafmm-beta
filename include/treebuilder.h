#ifndef treebuilder_h
#define treebuilder_h
#include "evaluator.h"

#if PARALLEL_EVERYTHING

struct BinaryTree {
  ivec8 counts;
  BinaryTree * left;
  BinaryTree * right;
};

struct XNode {
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
  int      MAXLEVEL;

protected:
  real_t localRadius;                                           //!< Radius of local root cell
  vec3 localCenter;                                             //!< Center of local root cell
  vec3 localXmin;                                               //!< Local Xmin for a given rank
  vec3 localXmax;                                               //!< Local Xmax for a given rank

private:
  void init(Node &node) {
    node.NCHILD = 0;
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
    N->NCHILD++;
    N->CHILD[octant] = nodes.size();
    nodes.push_back(child);
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
    if (!nodes[i].NCHILD) {
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
    nodes.reserve(NLEAF);
    Node node;
    init(node);
    node.LEVEL = 0;
    node.X = localCenter;
    nodes.push_back(node);
    MAXLEVEL = 0;
    for (L_iter L=leafs.begin(); L!=leafs.end(); L++) {
      int n = 0;
      while (nodes[n].NCHILD) {
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
    cells.resize(nodes.size());
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
  ivec8 prefixSum(ivec8 a, int begin) {
    ivec8 s;
    int p = begin;
    for (int i=0; i<8; i++) {
      s[i] = p;
      p += a[i];
    }
    return s;
  }

  XNode * makeNode(int level, int begin, int end, vec3 X, bool nochild) {
    XNode * n = new XNode();
    n->LEVEL = level; 
    n->BODY = begin; 
    n->NLEAF = end - begin; 
    n->X = X;
    n->NNODE = 1;
    if (nochild) {
      for (int k=0; k<8; k++) n->CHILD[k] = NULL;
    }
    return n;
  }

  int max_ivec8_nodes_to_count(int n) {
    if (n <= NSPAWN) return 1;
    else return 4 * ((n - 1) / NSPAWN) - 1;
  }

  int max_ivec8_nodes_to_build(int n) {
    return (4 * n) / NSPAWN;
  }

  BinaryTree * countBodies(Bodies& bodies, int begin, int end, vec3 X, 
                          BinaryTree * t_root, BinaryTree * t_begin, BinaryTree * t_end) {
    assert(max_ivec8_nodes_to_count(end - begin) <= t_end - t_begin + 1);
    if (end - begin <= NSPAWN) {
      for (int k=0; k<8; k++) t_root->counts[k] = 0;
      t_root->left = t_root->right = NULL;
      for (int i=begin; i<end; i++) {
        vec3 x = bodies[i].X;
        int oct = (x[0] > X[0]) + ((x[1] > X[1]) << 1) + ((x[2] > X[2]) << 2);
        t_root->counts[oct]++;
      } 
    } else {
      int mid = (begin + end) / 2;
      int n0 = max_ivec8_nodes_to_count(mid - begin);
      int n1 = max_ivec8_nodes_to_count(end - mid);
      BinaryTree * t_mid = t_begin + n0;
      assert(t_end - t_begin >= n0 + n1);
      __init_tasks__;
      spawn_task1(bodies, {
                  t_root->left
                  = countBodies(bodies, begin, mid, X, t_begin, t_begin + 1, t_begin + n0);
      });
      t_root->right
        = countBodies(bodies, mid, end, X, t_mid, t_mid + 1, t_mid + n1);
      __sync_tasks__;
      t_root->counts = t_root->left->counts + t_root->right->counts;
    }
    return t_root;
  }

  void moveBodies(Bodies& bodies, Bodies& buffer, int begin, int end, 
                  BinaryTree * t, ivec8 offsets, vec3 X) {
    if (t->left == NULL) {
      for (int i=begin; i<end; i++) {
        vec3 x = bodies[i].X;
        int oct = (x[0] > X[0]) + ((x[1] > X[1]) << 1) + ((x[2] > X[2]) << 2);
        buffer[offsets[oct]] = bodies[i];
        offsets[oct]++;
      }
    } else {
      int mid = (begin + end) / 2;
      ivec8 offsets_mid = offsets + t->left->counts;
      __init_tasks__;
      spawn_task2(bodies, buffer, {
                  moveBodies(bodies, buffer, begin, mid, 
                                   t->left, offsets, X);
      });
      moveBodies(bodies, buffer, mid, end, t->right, offsets_mid, X);
      __sync_tasks__;
    }
  }

  XNode * buildNodes(Bodies& bodies, Bodies& buffer, int dest,
                     int begin,  int end, 
                     BinaryTree * t_begin, BinaryTree * t_end, 
                     vec3 X, int level) {
    assert(max_ivec8_nodes_to_build(end - begin) <= t_end - t_begin);
    if (begin == end) return NULL;
    if (end - begin <= NCRIT) {
      if (dest)
        for (int i=begin; i<end; i++) buffer[i] = bodies[i];
      return makeNode(level, begin, end, X, true);
    }
    XNode * node = makeNode(level, begin, end, X, false);
    BinaryTree t_root_[1]; 
    BinaryTree * t_root = countBodies(bodies, begin, end, X, t_root_, t_begin, t_end);
    ivec8 offsets = prefixSum(t_root->counts, begin);
    moveBodies(bodies, buffer, begin, end, t_root, offsets, X);
    BinaryTree * t = t_begin;
    __init_tasks__;
    for (int k=0; k<8; k++) {
      int n_nodes = max_ivec8_nodes_to_build(t_root->counts[k]);
      assert(t + n_nodes <= t_end);
      spawn_task2(buffer, bodies, {
          vec3 Y = X;
          real_t r = localRadius / (1 << (level + 1));
          for (int d=0; d<3; d++) {
            Y[d] += r * (((k & 1 << d) >> d) * 2 - 1);
          }
          node->CHILD[k] 
            = buildNodes(buffer, bodies, 1 - dest,
                         offsets[k], offsets[k] + t_root->counts[k],
                         t, t + n_nodes,
                         Y, level + 1);
        });
      t += n_nodes;
    }
    __sync_tasks__;
    for (int k=0; k<8; k++) {
      if (node->CHILD[k]) {
        node->NNODE += node->CHILD[k]->NNODE;
      }
    }
    return node;
  }

  int nodes2cellsRec(XNode * node, int parent, C_iter C, C_iter H,
                     B_iter B0, int level) {
    C->PARENT = parent;
    C->R      = localRadius / (1 << node->LEVEL);
    C->X      = node->X;
    C->NDLEAF = node->NLEAF;
    C->LEAF   = B0 + node->BODY;
    if (node->NNODE == 1) {
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
      __init_tasks__;
      int levels[8];
      for (int k=0; k<nsub; k++) {
        int octant = child_octants[k];
        spawn_task1_if(node->NNODE > 1000,
                       levels,
                       levels[k] 
                       = nodes2cellsRec(node->CHILD[octant], C - Ci0, Ci, 
                                        H, B0, level + 1));
        Ci += 1;
        H += node->CHILD[octant]->NNODE - 1;
      }
      __sync_tasks__;
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
  void growTreeRec(Bodies &bodies) {
    Bodies buffer = bodies;
    startTimer("Grow tree");
    int n_nodes = max_ivec8_nodes_to_build(bodies.size());
    BinaryTree * t_begin = new BinaryTree[n_nodes];
    BinaryTree * t_end = t_begin + n_nodes;
    root_node = buildNodes(bodies, buffer, 0, 0, bodies.size(), 
                           t_begin, t_end, localCenter, 0);
    delete[] t_begin;
    stopTimer("Grow tree",printNow);
  }

  void linkTreeRec(Bodies& bodies, Cells &cells) {
    startTimer("Link tree");
    cells.resize(root_node->NNODE);
    Ci0 = cells.begin();
    MAXLEVEL = nodes2cellsRec(root_node, 0, Ci0, Ci0 + 1, bodies.begin(), 0);
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
    std::cout << "M2L calls            : " << NM2L                 << std::endl;
#endif
    std::cout << "-----------------------------------------------" << std::endl;
  }
};

#endif
