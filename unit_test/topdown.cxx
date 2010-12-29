#include "body.h"

struct node {
  int NLEAF;
  int NCHILD;
  int LEAF[NCRIT];
  node *CHILD[8];
  vect CENTER;

  void init(vect x) {
    NLEAF = 0;
    NCHILD = 0;
    for( int i=0; i!=NCRIT; ++i )
      LEAF[i] = 0;
    for( int i=0; i!=8; ++i )
      CHILD[i] = 0;
    CENTER = x;
  }

  void add_leaf(int const i) {
    LEAF[NLEAF] = i;
    NLEAF++;
  }

  int find_octant(vect const pos) {
    int I(pos[0] > CENTER[0] ? 1 : 0);
    if   (pos[1] > CENTER[1]) I |= 2;
    if   (pos[2] > CENTER[2]) I |= 4;
    return I;
  }

  void add_child(int const i, real r, node *&N) {
    vect new_center(CENTER);                                    // Initialize new center position
    if(i&1) new_center[0] += r;  else  new_center[0] -= r;      // Shift new center in x direction
    if(i&2) new_center[1] += r;  else  new_center[1] -= r;      // Shift new center in y direction
    if(i&4) new_center[2] += r;  else  new_center[2] -= r;      // Shift new center in z direction
    CHILD[i] = ++N;                                             // Increment node pointer and assign to child
    CHILD[i]->init(new_center);                                 // Initialize child node
    NCHILD |= (1 << i);                                         // Flip bit of octant
  }

  void split_node(real r, bodies *B, node *&N) {
    int octant;
    for( int i=0; i!=NCRIT; ++i ) {                             // Loop over all leafs in parent node
      octant = find_octant(B->pos(LEAF[i]));                    //  Find the octant where the body belongs
      if( !(NCHILD & (1 << octant)) )                           //  If child doesn't exist in this octant
        add_child(octant,r,N);                                  //   Add new child to list
      CHILD[octant]->add_leaf(LEAF[i]);                         //  Add leaf to child
      if( CHILD[octant]->NLEAF >= NCRIT ) {                     //  If there are still too many leafs
        r /= 2;                                                 //   New radius
        CHILD[octant]->split_node(r,B,N);                       //   Split the node into smaller ones
      }                                                         //  End if statement for splitting
    }                                                           // End loop over leafs
  }
};

void traverse(node *N, int &nodes, int &leafs) {
  if( N->NLEAF >= NCRIT ) {
    for( int i=0; i!=8; ++i )
      if( N->NCHILD & (1 << i) )
        traverse(N->CHILD[i],nodes,leafs);
  } else {
    for( int i=0; i!=N->NLEAF; ++i ) {
//      std::cout << nodes << " " << leafs << " " << N->LEAF[i] << std::endl;
      leafs++;
    }
    nodes++;
  }
}

int main(int argc, const char* argv[])
{
  double tic,toc;
  int Nbody=10000000;
  if(argc>1) Nbody = atoi(argv[1]);
  bodies B(Nbody,Nbody);

  tic = get_time();
  for( B=B.begin(); B!=B.end(); ++B ) {                         // Loop over all bodies
    for( int d=0; d!=3; ++d )                                   //  Loop over each dimension
      B.pos()[d] = rand()/(1.+RAND_MAX)*2-1;                    //   Initialize positions
    real r = sqrt(B.pos()[0]*B.pos()[0]+B.pos()[1]*B.pos()[1]+B.pos()[2]*B.pos()[2]);
    for( int d=0; d!=3; ++d )                                   //  Loop over each dimension
      B.pos()[d] /= r*1.1;                                      //   Normalize positions
    B.scal() = 1./B.size();                                     //  Initialize source value
  }                                                             // End loop over all bodies
  toc = get_time();
  std::cout << "Initialize    : " << toc-tic << std::endl;

  tic = get_time();
  real r0(0);                                                   // Root radius
  vect xmin,xmax,x0(0);                                         // Min,Max,Center of domain
  B.begin();                                                    // Reset bodies counter
  xmin = xmax = B.pos();                                        // Initialize xmin,xmax
  for( B=B.begin(); B!=B.end(); ++B ) {                         // Loop over all bodies
    for( int d=0; d!=3; ++d ) {                                 //  Loop over each dimension
      if     (B.pos()[d] < xmin[d]) xmin[d] = B.pos()[d];       //   Determine xmin
      else if(B.pos()[d] < xmax[d]) xmax[d] = B.pos()[d];       //   Determine xmax
    }                                                           //  End loop over each dimension
    x0 += B.pos();                                              //  Sum positions
  }                                                             // End loop over all bodies
  x0 /= B.size();                                               // Calculate average position
  for( int d=0; d!=3; ++d ) {                                   // Loop over each dimension
    x0[d] = int(x0[d]+0.5);                                     //  Shift center to nearest integer
    r0 = std::max(xmax[d]-x0[d],r0);                            //  Calculate max distance from center
    r0 = std::max(x0[d]-xmin[d],r0);                            //  Calculate max distance from center
  }                                                             // End loop over each dimension
  r0 = pow(2.0,int(1.0+log(r0)/M_LN2));                         // Add some leeway to root radius
  toc = get_time();
  std::cout << "Set domain    : " << toc-tic << std::endl;

  tic = get_time();
  int level;                                                    // Level of tree, root = 0
  int octant;                                                   // In which octant is the body located?
  real r;                                                       // Radius of node
  node *N0,*NN,*N;                                              // Declare pointers to nodes
  N0 = new node [Nbody];                                        // Allocate all nodes
  N0->init(x0);                                                 // Initialize root node
  NN = N0;                                                      // Keep copy for node counter
  for( B=B.begin(); B!=B.end(); ++B ) {                         // Loop over all bodies
    level = 0;                                                  //  Always start from root level
    for( N=N0; ; ) {                                            //  Loop over nodes
      if( N->NLEAF >= NCRIT ) {                                 //   If the node is not at the bottom of tree
        level++;                                                //    Increment the level of tree
        N->NLEAF++;                                             //    Increment the cumulative leaf counter
        octant = N->find_octant(B.pos());                       //    Find the octant where the body belongs
        if( !(N->NCHILD & (1 << octant)) ) {                    //    If child doesn't exist in this octant
          r = r0 / (1 << level);                                //     New radius
          N->add_child(octant,r,NN);                            //     Add new child to list
        }                                                       //    Endif for child existence
        N = N->CHILD[octant];                                   //     Update node pointer to child
      } else {                                                  //    If node is at the bottom of tree
        N->add_leaf(B);                                         //    Add body to node as leaf
        if( N->NLEAF >= NCRIT ) {                               //    If there are too many leafs
          level++;                                              //     Increment the level of tree
          r = r0 / (1 << level);                                //     New radius
          N->split_node(r,&B,NN);                               //     Split the node into smaller ones
        }                                                       //    Endif for splitting
        break;                                                  //    Exit the loop for nodes
      }                                                         //   Endif for child nodes
    }                                                           //  End loop over nodes
  }                                                             // End loop over all bodies
  toc = get_time();
  std::cout << "Construct tree: " << toc-tic << std::endl;

  int nodes(0),leafs(0);
  traverse(N0,nodes,leafs);
  std::cout << nodes << " "<< leafs << std::endl;
  delete[] N0;                                                  // Free all nodes
}
