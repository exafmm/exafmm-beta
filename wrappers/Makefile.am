include ../Makefile.am.include

if EXAFMM_HAVE_MPI
lib_LIBRARIES = libpetiga_laplace.a libpetiga_helmholtz.a libgromacs.a libcharmm.a libcharmm2.a
libpetiga_laplace_a_SOURCES = petiga.cxx ../kernels/LaplaceCartesianCPU.cxx ../kernels/LaplaceP2PCPU.cxx
libpetiga_laplace_a_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE -DEXAFMM_CARTESIAN -DEXAFMM_EXPANSION=4
libpetiga_helmholtz_a_SOURCES = petiga.cxx ../kernels/HelmholtzSphericalCPU.cxx ../kernels/HelmholtzP2PCPU.cxx
libpetiga_helmholtz_a_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_HELMHOLTZ -DEXAFMM_SPHERICAL -DEXAFMM_EXPANSION=10
libgromacs_a_SOURCES = gromacs.cxx ../kernels/LaplaceSphericalCPU.cxx ../kernels/LaplaceP2PCPU.cxx
libgromacs_a_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE -DEXAFMM_SPHERICAL -DEXAFMM_EXPANSION=9 -DEXAFMM_NO_P2P -fPIC
libcharmm_a_SOURCES = charmm.cxx ../kernels/LaplaceSphericalCPU.cxx ../kernels/LaplaceP2PCPU.cxx
libcharmm_a_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE -DEXAFMM_SPHERICAL -DEXAFMM_EXPANSION=9
libcharmm2_a_SOURCES = charmm2.cxx ../kernels/LaplaceSphericalCPU.cxx ../kernels/LaplaceP2PCPU.cxx
libcharmm2_a_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE -DEXAFMM_SPHERICAL -DEXAFMM_EXPANSION=9

bin_PROGRAMS = petiga_laplace_mpi petiga_helmholtz_mpi gromacs_mpi charmm_mpi charmm2_mpi
petiga_laplace_mpi_SOURCES = test_petiga.cxx
petiga_laplace_mpi_CPPFLAGS = $(libpetiga_laplace_a_CPPFLAGS)
petiga_laplace_mpi_LDADD = libpetiga_laplace.a
petiga_helmholtz_mpi_SOURCES = test_petiga.cxx
petiga_helmholtz_mpi_CPPFLAGS = $(libpetiga_helmholtz_a_CPPFLAGS)
petiga_helmholtz_mpi_LDADD = libpetiga_helmholtz.a
gromacs_mpi_SOURCES = test_gromacs.cxx
gromacs_mpi_CPPFLAGS = $(AM_CPPFLAGS)
gromacs_mpi_LDADD = libgromacs.a
charmm_mpi_SOURCES = test_charmm.f90
charmm_mpi_CPPFLAGS = $(AM_CPPFLAGS)
if EXAFMM_HAVE_FX
nodist_EXTRA_charmm_mpi_SOURCES = test_gromacs.cxx
charmm_mpi_LDADD = libcharmm.a --linkfortran
else
if EXAFMM_HAVE_SX
charmm_mpi_LDADD = libcharmm.a -lC++
else
charmm_mpi_LDADD = libcharmm.a -lstdc++ -lgomp
endif
endif
charmm2_mpi_SOURCES = test_charmm.f90
charmm2_mpi_CPPFLAGS = $(AM_CPPFLAGS)
if EXAFMM_HAVE_FX
nodist_EXTRA_charmm2_mpi_SOURCES = test_gromacs.cxx
charmm2_mpi_LDADD = libcharmm2.a --linkfortran
else
if EXAFMM_HAVE_SX
charmm2_mpi_LDADD = libcharmm2.a -lC++
else
charmm2_mpi_LDADD = libcharmm2.a -lstdc++ -lgomp
endif
endif

run_charmm_mpi: charmm_mpi
	./$< water_32.fmm water.out 0
run_charmm2_mpi: charmm2_mpi
	./$< water_32.fmm water.out 0
run_gromacs_mpi: gromacs_mpi
	$(MPIRUN) -n 2 ./$<
run_petiga_laplace_mpi: petiga_laplace_mpi
	$(MPIRUN) -n 2 ./$<
run_petiga_helmholtz_mpi: petiga_helmholtz_mpi
	$(MPIRUN) -n 1 ./$<
endif
