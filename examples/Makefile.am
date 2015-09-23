include ../Makefile.am.include

AM_CXXFLAGS  += -DEXPANSION=4 # Specify expansion order
AM_CXXFLAGS  += -DMASS # Use all positive sources
#AM_CXXFLAGS += -DIneJ # Different target and source points
AM_CXXFLAGS  += -DDIRECT # Compare with direct summation
#AM_CXXFLAGS += -DWRITE_TIME # Write timings to file
AM_CXXFLAGS  += -DFP64 # Use double precision

bin_PROGRAMS = laplace_cartesian laplace_spherical
laplace_cartesian_SOURCES = serial.cxx ../kernels/LaplaceCartesianCPU.cxx ../kernels/LaplaceP2PCPU.cxx
laplace_cartesian_CPPFLAGS = -DLaplace -DCartesian
laplace_spherical_SOURCES = serial.cxx ../kernels/LaplaceSphericalCPU.cxx ../kernels/LaplaceP2PCPU.cxx
laplace_spherical_CPPFLAGS = -DLaplace -DSpherical

bin_PROGRAMS += laplace_cartesian_mpi laplace_spherical_mpi
laplace_cartesian_mpi_SOURCES = parallel.cxx ../kernels/LaplaceCartesianCPU.cxx ../kernels/LaplaceP2PCPU.cxx
laplace_cartesian_mpi_CPPFLAGS = -DLaplace -DCartesian
laplace_spherical_mpi_SOURCES = parallel.cxx ../kernels/LaplaceSphericalCPU.cxx ../kernels/LaplaceP2PCPU.cxx
laplace_spherical_mpi_CPPFLAGS = -DLaplace -DSpherical


