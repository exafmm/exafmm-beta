include ../Makefile.am.include

if EXAFMM_HAVE_MPI
lib_LIBRARIES = libpetiga_exafmm_laplace.a libpetiga_exafmm_helmholtz.a libgromacs_exafmm.a libcharmm_exafmm.a libcharmm2_exafmm.a libmd_distributed.a
libpetiga_exafmm_laplace_a_SOURCES = petiga.cxx ../kernels/LaplaceCartesianCPU.cxx ../kernels/LaplaceP2PCPU.cxx
libpetiga_exafmm_laplace_a_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE -DEXAFMM_CARTESIAN -DEXAFMM_EXPANSION=4
libpetiga_exafmm_helmholtz_a_SOURCES = petiga.cxx ../kernels/HelmholtzSphericalCPU.cxx ../kernels/HelmholtzP2PCPU.cxx
libpetiga_exafmm_helmholtz_a_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_HELMHOLTZ -DEXAFMM_SPHERICAL -DEXAFMM_EXPANSION=10
libgromacs_exafmm_a_SOURCES = gromacs.cxx ../kernels/LaplaceSphericalCPU.cxx ../kernels/LaplaceP2PCPU.cxx
libgromacs_exafmm_a_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE -DEXAFMM_SPHERICAL -DEXAFMM_EXPANSION=9 -DEXAFMM_NO_P2P
libcharmm_exafmm_a_SOURCES = charmm.cxx ../kernels/LaplaceSphericalCPU.cxx ../kernels/LaplaceP2PCPU.cxx
libcharmm_exafmm_a_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE -DEXAFMM_SPHERICAL -DEXAFMM_EXPANSION=9
libcharmm2_exafmm_a_SOURCES = charmm2.cxx ../kernels/LaplaceSphericalCPU.cxx ../kernels/LaplaceP2PCPU.cxx
libcharmm2_exafmm_a_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE -DEXAFMM_SPHERICAL -DEXAFMM_EXPANSION=9
libmd_distributed_a_SOURCES = md_distributed.cxx ../kernels/LaplaceSphericalCPU.cxx ../kernels/LaplaceP2PCPU.cxx
libmd_distributed_a_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE -DEXAFMM_SPHERICAL -DEXAFMM_EXPANSION=9

bin_PROGRAMS = petiga_laplace_mpi petiga_helmholtz_mpi gromacs_mpi charmm_mpi charmm2_mpi md_distributed_mpi
petiga_laplace_mpi_SOURCES = test_petiga.cxx
petiga_laplace_mpi_CPPFLAGS = $(libpetiga_laplace_a_CPPFLAGS)
petiga_laplace_mpi_LDADD = libpetiga_exafmm_laplace.a
petiga_helmholtz_mpi_SOURCES = test_petiga.cxx
petiga_helmholtz_mpi_CPPFLAGS = $(libpetiga_helmholtz_a_CPPFLAGS)
petiga_helmholtz_mpi_LDADD = libpetiga_exafmm_helmholtz.a
gromacs_mpi_SOURCES = test_gromacs.cxx
gromacs_mpi_CPPFLAGS = $(AM_CPPFLAGS)
gromacs_mpi_LDADD = libgromacs_exafmm.a
charmm_mpi_SOURCES = test_charmm.f90
charmm_mpi_CPPFLAGS = $(AM_CPPFLAGS)
if EXAFMM_HAVE_FX
nodist_EXTRA_charmm_mpi_SOURCES = test_gromacs.cxx
charmm_mpi_LINK = $(CXXLINK) --linkfortran
charmm_mpi_LDADD = libcharmm_exafmm.a
else
charmm_mpi_LINK = $(FCLINK) $(AM_FCFLAGS)
if EXAFMM_HAVE_SX
charmm_mpi_LDADD = libcharmm_exafmm.a -lC++
else
charmm_mpi_LDADD = libcharmm_exafmm.a -lstdc++ -lgomp
endif
endif
charmm2_mpi_SOURCES = test_charmm.f90
charmm2_mpi_CPPFLAGS = $(AM_CPPFLAGS)
if EXAFMM_HAVE_FX
nodist_EXTRA_charmm2_mpi_SOURCES = test_gromacs.cxx
charmm2_mpi_LINK = $(CXXLINK) --linkfortran
charmm2_mpi_LDADD = libcharmm2_exafmm.a
else
charmm2_mpi_LINK = $(FCLINK) $(AM_FCFLAGS)
if EXAFMM_HAVE_SX
charmm2_mpi_LDADD = libcharmm2_exafmm.a -lC++
else
charmm2_mpi_LDADD = libcharmm2_exafmm.a -lstdc++ -lgomp
endif
endif
md_distributed_mpi_SOURCES = test_md_distributed.cxx
md_distributed_mpi_CPPFLAGS = $(AM_CPPFLAGS)
md_distributed_mpi_LDADD = libmd_distributed.a

run_charmm_mpi: charmm_mpi
	$(MPIRUN) -n 2 ./$< water_32.fmm water.out 0
run_charmm2_mpi: charmm2_mpi
	$(MPIRUN) -n 2 ./$< water_32.fmm water.out 0
run_gromacs_mpi: gromacs_mpi
	$(MPIRUN) -n 2 ./$<
run_petiga_laplace_mpi: petiga_laplace_mpi
	$(MPIRUN) -n 2 ./$<
run_petiga_helmholtz_mpi: petiga_helmholtz_mpi
	$(MPIRUN) -n 1 ./$<
run_md_distributed_mpi: md_distributed_mpi
	$(MPIRUN) -n 2 ./$<
endif
