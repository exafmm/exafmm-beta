#ifndef treebuilder_h
#define treebuilder_h
#include "evaluator.h"

class TreeBuilder : public Evaluator {
private:
  struct BinaryTreeNode {
    ivec8            NLEAF;                                     //!< Number of descendant leafs
    BinaryTreeNode * LEFT;                                      //!< Pointer to left child
    BinaryTreeNode * RIGHT;                                     //!< Pointer to right child
  };

  struct OctreeNode {
    int          LEAF;                                          //!< Index offset for first leaf in node
    int          NLEAF;                                         //!< Number of descendant leafs
    int          NNODE;                                         //!< Number of descendant nodes
    OctreeNode * CHILD[8];                                      //!< Pointer to child node
    vec3         X;                                             //!< Coordinate at center
  };

  int          MAXLEVEL;                                        //!< Maximum level of tree
  OctreeNode * N0;                                              //!< Octree root node

protected:
  real_t localRadius;                                           //!< Radius of local root cell
  vec3   localCenter;                                           //!< Center of local root cell
  vec3   localXmin;                                             //!< Local Xmin for a given rank
  vec3   localXmax;                                             //!< Local Xmax for a given rank

private:
  inline int getNumBinNode(int n) const {
    if (n <= NSPAWN) return 1;
    else return 4 * ((n - 1) / NSPAWN) - 1;
  }

  inline int getMaxBinNode(int n) const {
    return (4 * n) / NSPAWN;
  }

  ivec8 prefixSum(ivec8 a, int begin) {
    ivec8 s;
    int p = begin;
    for (int i=0; i<8; i++) {
      s[i] = p;
      p += a[i];
    }
    return s;
  }

  OctreeNode * makeNode(int begin, int end, vec3 X, bool nochild) {
    OctreeNode * octNode = new OctreeNode();
    octNode->LEAF = begin; 
    octNode->NLEAF = end - begin; 
    octNode->NNODE = 1;
    octNode->X = X;
    if (nochild) {
      for (int i=0; i<8; i++) octNode->CHILD[i] = NULL;
    }
    return octNode;
  }

  BinaryTreeNode * countBodies(Bodies& bodies, int begin, int end, vec3 X, 
                          BinaryTreeNode * t_root, BinaryTreeNode * t_begin, BinaryTreeNode * t_end) {
    assert(getNumBinNode(end - begin) <= t_end - t_begin + 1);
    if (end - begin <= NSPAWN) {
      for (int k=0; k<8; k++) t_root->NLEAF[k] = 0;
      t_root->LEFT = t_root->RIGHT = NULL;
      for (int i=begin; i<end; i++) {
        vec3 x = bodies[i].X;
        int oct = (x[0] > X[0]) + ((x[1] > X[1]) << 1) + ((x[2] > X[2]) << 2);
        t_root->NLEAF[oct]++;
      } 
    } else {
      int mid = (begin + end) / 2;
      int numLeftNode = getNumBinNode(mid - begin);
      int numRightNode = getNumBinNode(end - mid);
      BinaryTreeNode * t_mid = t_begin + numLeftNode;
      assert(t_end - t_begin >= numLeftNode + numRightNode);
      __init_tasks__;
      spawn_task1(bodies, {
                  t_root->LEFT
                  = countBodies(bodies, begin, mid, X, t_begin, t_begin + 1, t_begin + numLeftNode);
      });
      t_root->RIGHT
        = countBodies(bodies, mid, end, X, t_mid, t_mid + 1, t_mid + numRightNode);
      __sync_tasks__;
      t_root->NLEAF = t_root->LEFT->NLEAF + t_root->RIGHT->NLEAF;
    }
    return t_root;
  }

  void moveBodies(Bodies& bodies, Bodies& buffer, int begin, int end, 
                  BinaryTreeNode * t, ivec8 offsets, vec3 X) {
    if (t->LEFT == NULL) {
      for (int i=begin; i<end; i++) {
        vec3 x = bodies[i].X;
        int oct = (x[0] > X[0]) + ((x[1] > X[1]) << 1) + ((x[2] > X[2]) << 2);
        buffer[offsets[oct]] = bodies[i];
        offsets[oct]++;
      }
    } else {
      int mid = (begin + end) / 2;
      ivec8 offsets_mid = offsets + t->LEFT->NLEAF;
      __init_tasks__;
      spawn_task2(bodies, buffer, {
                  moveBodies(bodies, buffer, begin, mid, 
                                   t->LEFT, offsets, X);
      });
      moveBodies(bodies, buffer, mid, end, t->RIGHT, offsets_mid, X);
      __sync_tasks__;
    }
  }

  OctreeNode * buildNodes(Bodies& bodies, Bodies& buffer, int dest,
                     int begin,  int end, 
                     BinaryTreeNode * t_begin, BinaryTreeNode * t_end, 
                     vec3 X, int level) {
    assert(getMaxBinNode(end - begin) <= t_end - t_begin);
    if (begin == end) return NULL;
    if (end - begin <= NCRIT) {
      if (dest)
        for (int i=begin; i<end; i++) buffer[i] = bodies[i];
      return makeNode(begin, end, X, true);
    }
    OctreeNode * octNode = makeNode(begin, end, X, false);
    BinaryTreeNode t_root_[1]; 
    BinaryTreeNode * t_root = countBodies(bodies, begin, end, X, t_root_, t_begin, t_end);
    ivec8 offsets = prefixSum(t_root->NLEAF, begin);
    moveBodies(bodies, buffer, begin, end, t_root, offsets, X);
    BinaryTreeNode * t = t_begin;
    __init_tasks__;
    for (int k=0; k<8; k++) {
      int n_nodes = getMaxBinNode(t_root->NLEAF[k]);
      assert(t + n_nodes <= t_end);
      spawn_task2(buffer, bodies, {
          vec3 Y = X;
          real_t r = localRadius / (1 << (level + 1));
          for (int d=0; d<3; d++) {
            Y[d] += r * (((k & 1 << d) >> d) * 2 - 1);
          }
          octNode->CHILD[k] 
            = buildNodes(buffer, bodies, 1 - dest,
                         offsets[k], offsets[k] + t_root->NLEAF[k],
                         t, t + n_nodes,
                         Y, level + 1);
        });
      t += n_nodes;
    }
    __sync_tasks__;
    for (int k=0; k<8; k++) {
      if (octNode->CHILD[k]) {
        octNode->NNODE += octNode->CHILD[k]->NNODE;
      }
    }
    return octNode;
  }

  void nodes2cells(OctreeNode * octNode, int parent, C_iter C, C_iter H,
                  B_iter B0, int level) {
    C->PARENT = parent;
    C->R      = localRadius / (1 << level);
    C->X      = octNode->X;
    C->NDLEAF = octNode->NLEAF;
    C->LEAF   = B0 + octNode->LEAF;
    if (octNode->NNODE == 1) {
      C->CHILD = 0;
      C->NCHILD = 0;
      C->NCLEAF = octNode->NLEAF;
      assert(C->NCLEAF > 0);
      MAXLEVEL = std::max(MAXLEVEL,level);
    } else {
      C->NCLEAF = 0;
      int nsub=0;
      char child_octants[8];
      for (int octant=0; octant<8; octant++) {
        if (octNode->CHILD[octant]) {
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
      for (int k=0; k<nsub; k++) {
        int octant = child_octants[k];
        spawn_task0_if(octNode->NNODE > 1000,
                       nodes2cells(octNode->CHILD[octant], C - Ci0, Ci, 
                                   H, B0, level + 1));
        Ci += 1;
        H += octNode->CHILD[octant]->NNODE - 1;
      }
      __sync_tasks__;
      for (int k=0; k<nsub; k++) {
        int octant = child_octants[k];
        delete octNode->CHILD[octant];
      }
      if(nsub) MAXLEVEL = std::max(MAXLEVEL, level+1);
    }
  }

protected:
  void growTree(Bodies &bodies) {
    Bodies buffer = bodies;
    startTimer("Grow tree");
    int n_nodes = getMaxBinNode(bodies.size());
    BinaryTreeNode * t_begin = new BinaryTreeNode[n_nodes];
    BinaryTreeNode * t_end = t_begin + n_nodes;
    N0 = buildNodes(bodies, buffer, 0, 0, bodies.size(), 
                    t_begin, t_end, localCenter, 0);
    delete[] t_begin;
    stopTimer("Grow tree",printNow);
  }

  void linkTree(Bodies& bodies, Cells &cells) {
    startTimer("Link tree");
    cells.resize(N0->NNODE);
    Ci0 = cells.begin();
    nodes2cells(N0, 0, Ci0, Ci0 + 1, bodies.begin(), 0);
    delete N0; N0 = NULL;
    stopTimer("Link tree",printNow);
  }

  void printTreeData(Cells &cells) {
    std::cout << "----------------------------------" << std::endl;
    std::cout << "Bodies               : " << cells.front().NDLEAF << std::endl;
    std::cout << "Cells                : " << cells.size() << std::endl;
    std::cout << "Tree depth           : " << MAXLEVEL << std::endl;
#if COUNT
    std::cout << "P2P calls            : " << NP2P << std::endl;
    std::cout << "M2L calls            : " << NM2L << std::endl;
#endif
    std::cout << "----------------------------------" << std::endl;
  }
public:
  TreeBuilder() : MAXLEVEL(0) {}

};

#endif
