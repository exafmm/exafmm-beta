include ../Makefile.am.include

lib_LIBRARIES = libpetiga.a
libpetiga_a_SOURCES = petiga.cxx ../kernels/LaplaceCartesianCPU.cxx ../kernels/LaplaceP2PCPU.cxx
libpetiga_a_CPPFLAGS = $(AM_CPPFLAGS) -DLaplace -DCartesian -DEXPANSION=4
