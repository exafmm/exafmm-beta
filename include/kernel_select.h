#if EXAFMM_CARTESIAN
#include "LaplaceCartesianCPU.h"
#elif EXAFMM_SPHERICAL
#include "LaplaceSphericalCPU.h"
#endif
#include "HelmholtzSphericalCPU.h"
#include "BiotSavartSphericalCPU.h"
#if EXAFMM_LAPLACE
#if EXAFMM_CARTESIAN
typedef exafmm::LaplaceCartesianCPU kernel;
#elif EXAFMM_SPHERICAL
typedef exafmm::LaplaceSphericalCPU kernel;
#endif
#elif EXAFMM_HELMHOLTZ
typedef exafmm::HelmholtzSphericalCPU kernel;
#elif EXAFMM_BIOTSAVART
typedef exafmm::BiotSavartSphericalCPU kernel;
#else
#error Please specify the type of kernel (ex: -DEXAFMM_LAPLACE)
#endif
