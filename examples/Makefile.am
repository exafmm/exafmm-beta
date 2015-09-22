include ../Makefile.am.include

AM_CXXFLAGS  += -DEXPANSION=4 # Specify expansion order
AM_CXXFLAGS  += -DMASS # Use all positive sources
#AM_CXXFLAGS += -DIneJ # Different target and source points
AM_CXXFLAGS  += -DDIRECT # Compare with direct summation
#AM_CXXFLAGS += -DWRITE_TIME # Write timings to file
AM_CXXFLAGS  += -DFP64 # Use double precision

bin_PROGRAMS = serial parallel
serial_SOURCES = serial.cxx ../kernels/LaplaceCartesianCPU.cxx ../kernels/LaplaceP2PCPU.cxx
parallel_SOURCES = parallel.cxx ../kernels/LaplaceCartesianCPU.cxx ../kernels/LaplaceP2PCPU.cxx
