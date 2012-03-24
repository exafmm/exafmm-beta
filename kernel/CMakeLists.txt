
IF(USE_SPHERICAL)
    set(cpu_kernels CPUSphericalStokes.cxx CPUSphericalLaplace.cxx)
ELSE()
    set(cpu_kernels CPUCartesianLaplace.cxx)
ENDIF()

set(gpu_kernels GPUSphericalStokes.cu GPUSphericalLaplace.cu GPUVanDerWaals.cu)

IF(USE_CUDA)
    cuda_add_library(Kernels STATIC ${gpu_kernels} ${cpu_kernels} CPUP2P.cxx)
    TARGET_LINK_LIBRARIES(Kernels ${CUDA_LIBRARIES})
ELSE()
    ADD_LIBRARY(Kernels ${cpu_kernels} CPUP2P.cxx CPUVanDerWaals.cxx CPUEwaldLaplace.cxx)
ENDIF()
