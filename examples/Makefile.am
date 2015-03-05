AM_CPPFLAGS = -I../include -DCartesian -DEXPANSION=4 -DUSE_SIMD -DTBB -DMASS -DDIRECT
bin_PROGRAMS = serial
LDADD = -ltbb
serial_SOURCES = serial.cxx ../kernels/LaplaceCartesianCPU.cxx ../kernels/LaplaceP2PCPU.cxx
