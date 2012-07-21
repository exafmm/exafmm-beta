#ifndef serialfmm_h
#define serialfmm_h
#include "treebuilder.h"

class SerialFMM : public TreeBuilder {
public:
  void setBounds(Bodies &bodies) {
    startTimer("Set bounds");
    if( IMAGES != 0 ) {
      X0 = 0;
      R0 = M_PI;
    } else {
      vect xmin, xmax;
      X0 = 0;
      xmax = xmin = bodies.begin()->X;
      for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
        for( int d=0; d!=3; ++d ) {
          if     (B->X[d] < xmin[d]) xmin[d] = B->X[d];
          else if(B->X[d] > xmax[d]) xmax[d] = B->X[d];
        }
        X0 += B->X;
      }
      X0 /= bodies.size();
      for( int d=0; d!=3; ++d ) {
        R0 = std::max(xmax[d] - X0[d], R0);
        R0 = std::max(X0[d] - xmin[d], R0);
      }
      R0 *= 1.000001;
    }
    stopTimer("Set bounds",printNow);
  }

  void buildTree(Bodies &bodies, Cells &cells) {
    setLeafs(bodies);
    growTree();
    linkTree(bodies,cells);
  }

  void upwardPass(Cells &cells) {
    startTimer("Upward pass");
    setRootCell(cells);
    for( C_iter C=cells.end()-1; C!=cells.begin()-1; --C ) {
      real Rmax = 0;
      setCenter(C);
      C->M = 0;
      C->L = 0;
      P2M(C,Rmax);
      M2M(C,Rmax);
    }
#if Cartesian
    for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {
      for( int i=1; i<MTERM; ++i ) C->M[i] /= C->M[0];
    }
#endif
    setRcrit(cells);
    stopTimer("Upward pass",printNow);
  }

  void evaluate(Cells &cells) {
    setRootCell(cells);
    CellQueue cellQueue;
    pushCell(Ci0,cellQueue);
    Xperiodic = 0;
    startTimer("Traverse");
    traverse(cellQueue);
    stopTimer("Traverse",printNow);
  }

  void evaluate(Cells &icells, Cells &jcells) {
    setRootCell(icells,jcells);
    Pair pair(Ci0,Cj0);
    PairQueue pairQueue;
    startTimer("Traverse");
    if( IMAGES == 0 ) {
      Xperiodic = 0;
      pairQueue.push_back(pair);
      traverse(pairQueue);
    } else {
      for( int ix=-1; ix<=1; ++ix ) {                           //  Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //   Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz ) {                       //    Loop over z periodic direction
            Xperiodic[0] = ix * 2 * R0;                         //     Coordinate offset for x periodic direction
            Xperiodic[1] = iy * 2 * R0;                         //     Coordinate offset for y periodic direction
            Xperiodic[2] = iz * 2 * R0;                         //     Coordinate offset for z periodic direction
            pairQueue.push_back(pair);                          //     Push pair to queue
            traverse(pairQueue);                                //     Traverse a pair of trees
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      traversePeriodic();                                       //  Traverse tree for periodic images
    }
    stopTimer("Traverse",printNow);
  }

  void downwardPass(Cells &cells) {
    startTimer("Downward pass");
    C_iter C0 = cells.begin();
    L2P(C0);
    for( C_iter C=C0+1; C!=cells.end(); ++C ) {
      L2L(C);
      L2P(C);
    }
    stopTimer("Downward pass",printNow);
    if(printNow) printTreeData(cells);
  }

  void direct(Bodies &ibodies, Bodies &jbodies) {
    Cells cells;
    cells.resize(2);
    C_iter Ci = cells.begin(), Cj = cells.begin()+1;
    Ci->LEAF = ibodies.begin();
    Ci->NDLEAF = ibodies.size();
    Cj->LEAF = jbodies.begin();
    Cj->NDLEAF = jbodies.size();
    int prange = 0;                                             // Range of periodic images
    for( int i=0; i<IMAGES; i++ ) {                             // Loop over periodic image sublevels
      prange += int(pow(3,i));                                  //  Accumulate range of periodic images
    }                                                           // End loop over perioidc image sublevels
    for( int ix=-prange; ix<=prange; ++ix ) {                   // Loop over x periodic direction
      for( int iy=-prange; iy<=prange; ++iy ) {                 //  Loop over y periodic direction
        for( int iz=-prange; iz<=prange; ++iz ) {               //   Loop over z periodic direction
          Xperiodic[0] = ix * 2 * R0;                           //    Coordinate offset for x periodic direction
          Xperiodic[1] = iy * 2 * R0;                           //    Coordinate offset for y periodic direction
          Xperiodic[2] = iz * 2 * R0;                           //    Coordinate offset for z periodic direction
          P2P(Ci,Cj,false);                                     //    Evaluate P2P kernel
        }                                                       //   End loop over z periodic direction
      }                                                         //  End loop over y periodic direction
    }                                                           // End loop over x periodic direction
  }

  void normalize(Bodies &bodies) {
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {
      B->TRG /= B->SRC;
    }
  }
};

#endif
