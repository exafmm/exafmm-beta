include ../Makefile.am.include

lib_LIBRARIES = libpetiga.a libgromacs.a libcharmm.a libcharmm2.a
libpetiga_a_SOURCES = petiga.cxx ../kernels/LaplaceCartesianCPU.cxx ../kernels/LaplaceP2PCPU.cxx
libpetiga_a_CPPFLAGS = $(AM_CPPFLAGS) -DLaplace -DCartesian -DEXPANSION=4
libgromacs_a_SOURCES = gromacs.cxx ../kernels/LaplaceSphericalCPU.cxx ../kernels/LaplaceP2PCPU.cxx
libgromacs_a_CPPFLAGS = $(AM_CPPFLAGS) -DLaplace -DSpherical -DEXPANSION=9
libcharmm_a_SOURCES = charmm.cxx ../kernels/LaplaceSphericalCPU.cxx ../kernels/LaplaceP2PCPU.cxx
libcharmm_a_CPPFLAGS = $(AM_CPPFLAGS) -DLaplace -DSpherical -DEXPANSION=9
libcharmm2_a_SOURCES = charmm2.cxx ../kernels/LaplaceSphericalCPU.cxx ../kernels/LaplaceP2PCPU.cxx
libcharmm2_a_CPPFLAGS = $(AM_CPPFLAGS) -DLaplace -DSpherical -DEXPANSION=9

bin_PROGRAMS = petiga gromacs charmm charmm2
petiga_SOURCES = test_petiga.cxx
petiga_CPPFLAGS = $(AM_CPPFLAGS)
petiga_LDADD = libpetiga.a
gromacs_SOURCES = test_gromacs.cxx
gromacs_CPPFLAGS = $(AM_CPPFLAGS)
gromacs_LDADD = libgromacs.a
charmm_SOURCES = test_charmm.f90
charmm_CPPFLAGS = $(AM_CPPFLAGS)
charmm_LDADD = libcharmm.a
charmm2_SOURCES = test_charmm.f90
charmm2_CPPFLAGS = $(AM_CPPFLAGS)
charmm2_LDADD = libcharmm2.a


