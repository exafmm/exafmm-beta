include ../Makefile.am.include

AM_CXXFLAGS  += -DDIRECT # Compare with direct summation
#AM_CXXFLAGS += -DWRITE_TIME # Write timings to file

bin_PROGRAMS = laplace_cartesian laplace_spherical tree
laplace_cartesian_SOURCES = fmm.cxx ../kernels/LaplaceCartesianCPU.cxx ../kernels/LaplaceP2PCPU.cxx
laplace_cartesian_CPPFLAGS = $(AM_CPPFLAGS) -DLaplace -DCartesian -DMASS -DEXPANSION=4
laplace_spherical_SOURCES = fmm.cxx ../kernels/LaplaceSphericalCPU.cxx ../kernels/LaplaceP2PCPU.cxx
laplace_spherical_CPPFLAGS = $(AM_CPPFLAGS) -DLaplace -DSpherical -DEXPANSION=10
tree_SOURCES = tree.cxx
tree_CPPFLAGS = $(laplace_cartesian_CPPFLAGS)

bin_PROGRAMS += laplace_cartesian_kernel laplace_spherical_kernel
laplace_cartesian_kernel_SOURCES = kernel.cxx ../kernels/LaplaceCartesianCPU.cxx ../kernels/LaplaceP2PCPU.cxx
laplace_cartesian_kernel_CPPFLAGS = $(laplace_cartesian_CPPFLAGS)
laplace_spherical_kernel_SOURCES = kernel.cxx ../kernels/LaplaceSphericalCPU.cxx ../kernels/LaplaceP2PCPU.cxx
laplace_spherical_kernel_CPPFLAGS = $(laplace_spherical_CPPFLAGS)

if MPI
bin_PROGRAMS += laplace_cartesian_mpi laplace_spherical_mpi ewald_mpi
laplace_cartesian_mpi_SOURCES = fmm_mpi.cxx ../kernels/LaplaceCartesianCPU.cxx ../kernels/LaplaceP2PCPU.cxx
laplace_cartesian_mpi_CPPFLAGS = $(laplace_cartesian_CPPFLAGS)
laplace_spherical_mpi_SOURCES = fmm_mpi.cxx ../kernels/LaplaceSphericalCPU.cxx ../kernels/LaplaceP2PCPU.cxx
laplace_spherical_mpi_CPPFLAGS = $(laplace_spherical_CPPFLAGS)
ewald_mpi_SOURCES = ewald.cxx ../kernels/LaplaceSphericalCPU.cxx ../kernels/LaplaceP2PCPU.cxx
ewald_mpi_CPPFLAGS = $(laplace_spherical_CPPFLAGS)
endif
