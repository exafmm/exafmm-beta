#include "octree.h"

int main() {
  for( int it=0; it<25; it++ ) {
    uint numBodies = uint(pow(10,(it+24)/8.0));
    octree *tree = new octree(numBodies);
    printf("N     : %d\n",numBodies);
    for( uint i=0; i<numBodies; i++ ) {
      tree->bodyPos[i].w  = 1. / numBodies;
      tree->bodyPos[i].x  = drand48();
      tree->bodyPos[i].y  = drand48();
      tree->bodyPos[i].z  = drand48();
    }
    tree->bodyPos.h2d();
    tree->iterate(); 
    double tic = tree->get_time();
    tree->direct();
    double toc = tree->get_time();
    tree->bodyAcc.d2h();
    tree->bodyAcc2.d2h();
    float diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;
    for( uint i=0; i<numBodies/100; i++ ) {
      float4 fapprox = tree->bodyAcc[i];
      float4 fdirect = tree->bodyAcc2[i];
      diff1 += (fapprox.w - fdirect.w) * (fapprox.w - fdirect.w);
      diff2 += (fapprox.x - fdirect.x) * (fapprox.x - fdirect.x);
      diff2 += (fapprox.y - fdirect.y) * (fapprox.y - fdirect.y);
      diff2 += (fapprox.z - fdirect.z) * (fapprox.z - fdirect.z);
      norm1 += fdirect.w * fdirect.w;
      norm2 += fdirect.x * fdirect.x;
      norm2 += fdirect.y * fdirect.y;
      norm2 += fdirect.z * fdirect.z;
    }
    printf("Direct: %lf\n",toc-tic);
    printf("P Err : %f\n",sqrtf(diff1/norm1));
    printf("F Err : %f\n",sqrtf(diff2/norm2));
    delete tree;
  }
  return 0;
}
