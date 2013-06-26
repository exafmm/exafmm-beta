/*
Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/
#include "evaluator.h"

int main() {
  const int numBodies = 100;                                    // Number of bodies
  IMAGES = 0;                                                   // Level of periodic image tree (0 for non-periodic)
  THETA = 1 / sqrtf(4);                                         // Multipole acceptance criteria
  Bodies ibodies(numBodies);                                    // Define vector of target bodies
  Bodies ibodies2(numBodies);                                   // Define another vector of target bodies
  Bodies jbodies(numBodies);                                    // Define vector of source bodies
  Cells  icells;                                                // Define vector of target cells
  Cells  jcells;                                                // Define vector of source cells
  Evaluator<Laplace> FMM;                                       // Instantiate Evaluator class
  FMM.initialize();                                             // Initialize FMM
  FMM.preCalculation();                                         // Kernel pre-processing

  for( int it=0; it!=5; ++it ) {                                // Loop over kernel iterations
    real_t dist = (1 << it) / 2;                                //  Distance between source and target
    for( B_iter B=ibodies.begin(); B!=ibodies.end(); ++B ) {    //  Loop over target bodies
      for( int d=0; d!=3; ++d ) {                               //   Loop over dimensions
        B->X[d] = -drand48() - dist;                            //    Initialize positions
      }                                                         //   End loop over dimensions
    }                                                           //  End loop over target bodies
    for( B_iter B=jbodies.begin(); B!=jbodies.end(); ++B ) {    //  Loop over source bodies
      for( int d=0; d!=3; ++d ) {                               //   Loop over dimensions
        B->X[d] = drand48();                                    //    Initialize positions
      }                                                         //   End loop over dimensions
    }                                                           //  End loop over source bodies
    FMM.initSource(jbodies);                                    //  Initialize source values
    bool IeqJ = false;                                          //  Target == Source ?
    FMM.initTarget(ibodies,IeqJ);                               //  Initialize target values

    Cell cell;                                                  //  Define cell
    cell.NDLEAF   = numBodies;                                  //  Number of leafs
    cell.LEAF     = jbodies.begin();                            //  Iterator of first leaf
    cell.X        = 0.5;                                        //  Position
    cell.M        = 0;                                          //  Initialize multipole coefficients
    cell.ICELL    = 8;                                          //  Cell index
    cell.NCHILD   = 0;                                          //  Number of child cells
    cell.PARENT   = 1;                                          //  Iterator offset of parent cell
    jcells.push_back(cell);                                     //  Push cell into source cell vector
    FMM.evalP2M(jcells);                                        //  Evaluate P2M kernel
    cell.X        = 1;                                          //  Position
    cell.M        = 0;                                          //  Initialize multipole coefficients
    cell.ICELL    = 0;                                          //  Cell index
    cell.NCHILD   = 1;                                          //  Number of child cells
    cell.CHILD    = 0;                                          //  Iterator offset of child cell
    jcells.push_back(cell);                                     //  Push cell into source cell vector
    FMM.evalM2M(jcells,jcells);                                 //  Evaluate M2M kernel
    jcells.erase(jcells.begin());                               //  Erase first source cell
    cell.X        = -1 - dist;                                  //  Position
    cell.M        = 1;                                          //  Initialize multipole coefficients
    cell.L        = 0;                                          //  Initialize local coefficients
    cell.ICELL    = 0;                                          //  Cell index
    icells.push_back(cell);                                     //  Push cell into target cell vector
    FMM.addM2L(jcells.begin());                                 //  Add source cell to M2L list
    FMM.evalM2L(icells);                                        //  Evaluate M2L kernel
    cell.NDLEAF    = numBodies;                                 //  Number of leafs
    cell.LEAF     = ibodies.begin();                            //  Iterator of first leaf
    cell.X        = -0.5 - dist;                                //  Position
    cell.L        = 0;                                          //  Initialize local coefficients
    cell.ICELL    = 1;                                          //  Cell index
    cell.NCHILD   = 0;                                          //  Number of child cells
    cell.PARENT   = 1;                                          //  Iterator offset of parent cell
    icells.insert(icells.begin(),cell);                         //  Insert cell to begining of target cell vector
    FMM.evalL2L(icells);                                        //  Evaluate L2L kernel
    icells.pop_back();                                          //  Pop back target cell vector
    FMM.evalL2P(icells);                                        //  Evaluate L2P kernel

    ibodies2 = ibodies;                                         //  Copy target bodies
    FMM.initTarget(ibodies2,IeqJ);                              //  Reinitialize target values
    FMM.evalP2P(ibodies2,jbodies);                              //  Evaluate P2P kernel

    real_t diff1 = 0, norm1 = 0, diff2 = 0, norm2 = 0;          //  Initialize accumulators
    FMM.evalError(ibodies,ibodies2,diff1,norm1,diff2,norm2);    //  Evaluate error
    std::cout << "Distance      : " << dist << std::endl;       //  Print distance between target and source
    FMM.printError(diff1,norm1,diff2,norm2);                    //  Print the L2 norm error

    FMM.initTarget(ibodies);                                    //  Reinitialize target values
    FMM.addM2P(jcells.begin());                                 //  Add source cell to M2P list
    FMM.evalM2P(icells);                                        //  Evaluate M2P kernel
    icells.clear();                                             //  Clear target cell vector
    jcells.clear();                                             //  Clear source cell vector
    diff1 = norm1 = diff2 = norm2 = 0;                          //  Reinitialize accumulators
    FMM.evalError(ibodies,ibodies2,diff1,norm1,diff2,norm2);    //  Evaluate error
    FMM.printError(diff1,norm1,diff2,norm2);                    //  Print the L2 norm error
  }                                                             // End loop over kernel iterations
  FMM.postCalculation();                                        // Kernel post-processing
  FMM.finalize();                                               // Finalize FMM
}
