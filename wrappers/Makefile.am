include ../Makefile.am.include

if EXAFMM_HAVE_MPI
lib_LIBRARIES = libpetiga.a libgromacs.a libcharmm.a libcharmm2.a
libpetiga_a_SOURCES = petiga.cxx ../kernels/LaplaceCartesianCPU.cxx ../kernels/LaplaceP2PCPU.cxx
libpetiga_a_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE -DEXAFMM_CARTESIAN -DEXAFMM_EXPANSION=4
libgromacs_a_SOURCES = gromacs.cxx ../kernels/LaplaceSphericalCPU.cxx ../kernels/LaplaceP2PCPU.cxx
libgromacs_a_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE -DEXAFMM_SPHERICAL -DEXAFMM_EXPANSION=9 -DEXAFMM_NO_P2P
libcharmm_a_SOURCES = charmm.cxx ../kernels/LaplaceSphericalCPU.cxx ../kernels/LaplaceP2PCPU.cxx
libcharmm_a_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE -DEXAFMM_SPHERICAL -DEXAFMM_EXPANSION=9
libcharmm2_a_SOURCES = charmm2.cxx ../kernels/LaplaceSphericalCPU.cxx ../kernels/LaplaceP2PCPU.cxx
libcharmm2_a_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE -DEXAFMM_SPHERICAL -DEXAFMM_EXPANSION=9

bin_PROGRAMS = petiga gromacs charmm charmm2
petiga_SOURCES = test_petiga.cxx
petiga_CPPFLAGS = $(AM_CPPFLAGS)
petiga_LDADD = libpetiga.a
gromacs_SOURCES = test_gromacs.cxx
gromacs_CPPFLAGS = $(AM_CPPFLAGS)
gromacs_LDADD = libgromacs.a
charmm_SOURCES = test_charmm.f90
charmm_CPPFLAGS = $(AM_CPPFLAGS)
charmm_LDADD = libcharmm.a -lstdc++ -lgomp
charmm2_SOURCES = test_charmm.f90
charmm2_CPPFLAGS = $(AM_CPPFLAGS)
charmm2_LDADD = libcharmm2.a -lstdc++ -lgomp

run_charmm: charmm
	./$< water_32.fmm water.out 0
run_charmm2: charmm2
	./$< water_32.fmm water.out 0
run_gromacs: gromacs
	$(MPIRUN) -n 2 ./$<
run_petiga: petiga
	$(MPIRUN) -n 2 ./$<
endif
