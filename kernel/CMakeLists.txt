
IF(USE_SPHERICAL)
    set(cpu_kernels CPUP2P.cxx CPUSphericalLaplace.cxx CPUSphericalStokes.cxx)
ELSE()
    set(cpu_kernels CPUP2P.cxx CPUCartesianLaplace.cxx)
ENDIF()

ADD_LIBRARY(cpuKernel ${cpu_kernels} CPUVanDerWaals.cxx CPUEwaldLaplace.cxx)

IF(USE_GPU)
    SET(CUDA_TARGETS_SOURCES GPUSphericalLaplace.cu GPUSphericalStokes.cu GPUVanDerWaals.cu CPUP2P.cxx)
    cuda_add_library(gpuKernel STATIC ${CUDA_TARGETS_SOURCES})
    TARGET_LINK_LIBRARIES(gpuKernel ${CUDA_LIBRARIES})
ENDIF(USE_GPU)
