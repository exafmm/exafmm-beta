#include <fstream>
#include "kernel.h"
#include <vector>
#if VTK
#include "vtk.h"
#endif
#if COMkernel
#error Turn off COMkernel for this test
#endif

int main() {
  Bodies bodies(3);
  Cells cells(4);
  Kernel kernel;
  const real_t THETA = 0.5;
  const real_t R = 2 / THETA;

  B_iter Bj = bodies.begin();
  Bj->X = 2;
  Bj->SRC = 1;
  C_iter Cj = cells.begin();
  Cj->X = 1;
  Cj->BODY = Bj;
  Cj->NCBODY = 1;
  Cj->M = 0;
  kernel.P2M(Cj);

#if 1
  C_iter CJ = cells.begin()+1;
  CJ->CHILD = Cj-cells.begin();
  CJ->NCHILD = 1;
  CJ->X = 0;
  CJ->M = 0;
  kernel.M2M(CJ,cells.begin());

  C_iter CI = cells.begin()+2;
  CI->X = R + 4;
  CI->M = 1;
  CI->L = 0;
  kernel.M2L(CI,CJ,false);

  C_iter Ci = cells.begin()+3;
  Ci->X = R + 3;
  Ci->PARENT = 2;
  Ci->M = 1;
  Ci->L = 0;
  kernel.L2L(Ci,cells.begin());
#else
  C_iter Ci = cells.begin()+3;
  Ci->X = R + 3;
  Ci->M = 1;
  Ci->L = 0;
  kernel.M2L(Ci,Cj,false);
#endif

  B_iter Bi = bodies.begin()+1;
  Bi->X = R + 2;
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

  std::fstream file;
  file.open("kernel.dat", std::ios::out | std::ios::app);
  double err = std::abs((Bi->TRG[0] - Bi2->TRG[0])/Bi2->TRG[0]);
  file << P << " " << err << std::endl;
  file.close();

#if VTK
  file.open("kernel.dat", std::ios::in);
  std::string line;
  std::vector<int> order;
  std::vector<double> error;
  std::vector<double> bound;
  while (std::getline(file,line)) {
    int p;
    std::stringstream stream(line);
    stream >> p >>  err;
    order.push_back(p);
    error.push_back(log10(err));
    bound.push_back(log10(std::pow(THETA,p)));
  }
  file.close();

  if (order.size() > 5) {
    vtk2DPlot vtk;
    vtk.setName("order");
    vtk.setName("log10(error)");
    vtk.setName("log10(bound)");
    vtk.setNumRows(order.size());
    vtk.setData(0,order.size(),&order[0]);
    vtk.setData(0,error.size(),&error[0]);
    vtk.setData(0,bound.size(),&bound[0]);
    vtk.plot();
  }
#endif
  return 0;
}
