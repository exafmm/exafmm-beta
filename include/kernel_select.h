#if EXAFMM_LAPLACE
#if EXAFMM_CARTESIAN
#include "LaplaceCartesianCPU.h"
typedef exafmm::LaplaceCartesianCPU kernel;
#elif EXAFMM_SPHERICAL
#include "LaplaceSphericalCPU.h"
typedef exafmm::LaplaceSphericalCPU kernel;
#endif
#elif EXAFMM_HELMHOLTZ
#if DRR
#include "HelmholtzDrrSphericalCPU.h"
typedef exafmm::HelmholtzDrrSphericalCPU kernel;
#else
#include "HelmholtzSphericalCPU.h"
typedef exafmm::HelmholtzSphericalCPU kernel;
#endif
#elif EXAFMM_BIOTSAVART
#include "BiotSavartSphericalCPU.h"
typedef exafmm::BiotSavartSphericalCPU kernel;
#else
#error No valid kernel flag specified (ex: -DEXAFMM_BIOTSAVART)
#endif
