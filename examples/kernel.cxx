#include "kernel.h"
#if VTK
#include "vtk.h"
#endif

int main() {
  Bodies bodies(3);
  Cells cells(4);
  Kernel kernel;

  B_iter Bj = bodies.begin();
  Bj->X = 2;
  Bj->SRC = 1;
  C_iter Cj = cells.begin();
  Cj->X = 1;
  Cj->BODY = Bj;
  Cj->NCBODY = 1;
  Cj->M = 0;
  kernel.P2M(Cj);

#if 0
  C_iter CJ = cells.begin()+1;
  CJ->CHILD = Cj-cells.begin();
  CJ->NCHILD = 1;
  CJ->X = 0;
  CJ->M = 0;
  kernel.M2M(CJ,cells.begin());

  C_iter CI = cells.begin()+2;
  CI->X = 8;
  CI->M = 1;
  CI->L = 0;
  kernel.M2L(CI,CJ,false);

  C_iter Ci = cells.begin()+3;
  Ci->X = 7;
  Ci->PARENT = 2;
  Ci->M = 1;
  Ci->L = 0;
  kernel.L2L(Ci,cells.begin());
#else
  C_iter Ci = cells.begin()+3;
  Ci->X = 7;
  Ci->M = 1;
  Ci->L = 0;
  kernel.M2L(Ci,Cj,false);
#endif

  B_iter Bi = bodies.begin()+1;
  Bi->X = 6;
  Bi->SRC = 1;
  Bi->TRG = 0;
  Ci->BODY = Bi;
  Ci->NCBODY = 1;
  kernel.L2P(Ci);

  B_iter Bi2 = bodies.begin()+2;
  *Bi2 = *Bi;
  Bi2->TRG = 0;
  Cj->NDBODY = 1;
  Ci->NDBODY = 1;
  Ci->BODY = Bi2;
  kernel.P2P(Ci,Cj,false);

  for (int d=0; d<4; d++) {
    std::cout << Bi->TRG[d] << " " << Bi2->TRG[d] << std::endl;
  }

#if VTK
  vtk2DPlot vtk;
  const int numRows = 69;
  real_t * x = new real_t [numRows];
  real_t * y = new real_t [numRows];
  real_t * z = new real_t [numRows];
  real_t inc = 7.5 / (numRows-1);
  for (int i=0; i<numRows; i++) {
    x[i] = i * inc;
    y[i] = cos(i * inc);
    z[i] = sin(i * inc);
  }
  vtk.setName("x");
  vtk.setName("cos");
  vtk.setName("sin");

  vtk.setNumRows(numRows);
  vtk.setData(0,numRows,x);
  vtk.setData(0,numRows,y);
  vtk.setData(0,numRows,z);
  vtk.plot();

  delete[] x;
  delete[] y;
  delete[] z;
#endif
  return 0;
}
