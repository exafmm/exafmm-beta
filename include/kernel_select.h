#if EXAFMM_CARTESIAN
#include "laplace_cartesian_cpu.h"
#elif EXAFMM_SPHERICAL
#include "laplace_spherical_cpu.h"
#endif
#include "helmholtz_spherical_cpu.h"
#include "biot_savart_spherical_cpu.h"
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
