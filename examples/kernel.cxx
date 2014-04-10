#include <fstream>
#include "kernel.h"
#include <vector>
#if VTK
#include "vtk.h"
#endif

int main() {
  Bodies bodies(16), bodies2(16), jbodies(16);
  Cells cells(4);
  vec3 Xperiodic = 0;
  const real_t theta = 0.5;
  const real_t R = 2 / theta;

  for (B_iter B=jbodies.begin(); B!=jbodies.end(); B++) {
    B->X[0] = 2 * drand48();
    B->X[1] = 2 * drand48();
    B->X[2] = 2 * drand48();
    B->SRC  = drand48();
  }
  C_iter Cj = cells.begin();
  Cj->X = 1;
  Cj->BODY = jbodies.begin();
  Cj->NBODY = jbodies.size();
  Cj->M = 0;
  kernel::P2M(Cj);

#if 1
  C_iter CJ = cells.begin()+1;
  CJ->ICHILD = Cj-cells.begin();
  CJ->NCHILD = 1;
  CJ->X = 0;
  CJ->M = 0;
  kernel::M2M(CJ, cells.begin());

  C_iter CI = cells.begin()+2;
  CI->X = R + 4;
  CI->M = 1;
  CI->L = 0;
#if MASS
  for (int i=1; i<NTERM; i++) CJ->M[i] /= CJ->M[0];
#endif
  kernel::M2L(CI, CJ, Xperiodic, false);

  C_iter Ci = cells.begin()+3;
  Ci->X = R + 3;
  Ci->IPARENT = 2;
  Ci->M = 1;
  Ci->L = 0;
  kernel::L2L(Ci, cells.begin());
#else
  C_iter Ci = cells.begin()+3;
  Ci->X = R + 3;
  Ci->M = 1;
  Ci->L = 0;
#if MASS
  for (int i=1; i<NTERM; i++) Cj->M[i] /= Cj->M[0];
#endif
  kernel::M2L(Ci, Cj, Xperiodic, false);
#endif

  for (B_iter B=bodies.begin(); B!=bodies.end(); B++) {
    B->X[0] = R + 2 + 2 * drand48();
    B->X[1] = R + 2 + 2 * drand48();
    B->X[2] = R + 2 + 2 * drand48();
    B->SRC  = drand48();
    B->TRG  = 0;
  }
  Ci->BODY = bodies.begin();
  Ci->NBODY = bodies.size();
  kernel::L2P(Ci);

  for (B_iter B=bodies2.begin(); B!=bodies2.end(); B++) {
    *B = bodies[B-bodies2.begin()];
    B->TRG = 0;
  }
  Cj->NBODY = jbodies.size();
  Ci->NBODY = bodies2.size();
  Ci->BODY = bodies2.begin();
  kernel::P2P(Ci, Cj, Xperiodic, false);
  for (B_iter B=bodies2.begin(); B!=bodies2.end(); B++) {
    B->TRG /= B->SRC;
  }

  std::fstream file;
  file.open("kernel.dat", std::ios::out | std::ios::app);
  double diff = 0, norm = 0;
  for (B_iter B=bodies.begin(),B2=bodies2.begin(); B!=bodies.end(); B++,B2++) {
    diff += (B->TRG[0] - B2->TRG[0]) * (B->TRG[0] - B2->TRG[0]);
    norm += B2->TRG[0] * B2->TRG[0];
  }
  double err = std::sqrt(diff/norm);
  std::cout << P << " " << err << std::endl;
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
    double e;
    std::stringstream stream(line);
    stream >> p >> e;
    order.push_back(p);
    error.push_back(log10(e));
    bound.push_back(log10(std::pow(theta,p)));
  }
  file.close();

  if (order.size() > 12) {
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
