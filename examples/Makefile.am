AM_CXXFLAGS = -I../include -msse4a -ffast-math -funroll-loops -fopenmp -DCartesian -DEXPANSION=4 -DUSE_SIMD -DMASS -DDIRECT
# -DTBB
bin_PROGRAMS = serial parallel
serial_SOURCES = serial.cxx ../kernels/LaplaceCartesianCPU.cxx ../kernels/LaplaceP2PCPU.cxx
# serial_LDADD = -ltbb
parallel_SOURCES = parallel.cxx ../kernels/LaplaceCartesianCPU.cxx ../kernels/LaplaceP2PCPU.cxx
# parallel_LDADD = -ltbb
