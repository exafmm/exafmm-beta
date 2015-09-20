include ../Makefile.am.include

LFLAGS  += -DEXPANSION=4 # Specify expansion order
LFLAGS  += -DMASS # Use all positive sources
#LFLAGS += -DIneJ # Different target and source points
LFLAGS  += -DDIRECT # Compare with direct summation
#LFLAGS += -DWRITE_TIME # Write timings to file
LFLAGS  += -DFP64 # Use double precision

bin_PROGRAMS = serial parallel
serial_SOURCES = serial.cxx ../kernels/LaplaceCartesianCPU.cxx ../kernels/LaplaceP2PCPU.cxx
parallel_SOURCES = parallel.cxx ../kernels/LaplaceCartesianCPU.cxx ../kernels/LaplaceP2PCPU.cxx
