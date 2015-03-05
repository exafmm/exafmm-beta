AM_CXXFLAGS = -I../include -msse4a -ffast-math -funroll-loops -fopenmp -DCartesian -DEXPANSION=4 -DUSE_SIMD -DTBB -DMASS -DDIRECT
bin_PROGRAMS = serial
serial_SOURCES = serial.cxx ../kernels/LaplaceCartesianCPU.cxx ../kernels/LaplaceP2PCPU.cxx
serial_LDADD = -ltbb
