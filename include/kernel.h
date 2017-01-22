#if EXAFMM_LAPLACE
#include "laplace.h"
namespace exafmm {
  typedef LaplaceKernel<Pmax> Kernel;
}
#elif EXAFMM_HELMHOLTZ
#include "helmholtz.h"
namespace exafmm {
  typedef HelmholtzKernel<Pmax> Kernel;
}
#elif EXAFMM_BIOTSAVART
#include "biot_savart.h"
namespace exafmm {
  typedef BiotSavartKernel<Pmax> Kernel;
}
#endif
