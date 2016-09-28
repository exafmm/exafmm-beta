#if EXAFMM_LAPLACE
#if EXAFMM_CARTESIAN
#include "LaplaceCartesianCPU.h"
typedef exafmm::LaplaceCartesianCPU kernel;
#elif EXAFMM_SPHERICAL
#include "LaplaceSphericalCPU.h"
typedef exafmm::LaplaceSphericalCPU kernel;
#endif
#elif EXAFMM_HELMHOLTZ
#include "HelmholtzSphericalCPU.h"
typedef exafmm::HelmholtzSphericalCPU kernel;
#elif EXAFMM_BIOTSAVART
#include "BiotSavartSphericalCPU.h"
typedef exafmm::BiotSavartSphericalCPU kernel;
#else
#error Please specify the type of kernel (ex: -DEXAFMM_LAPLACE)
#endif
