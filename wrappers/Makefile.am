include ../Makefile.am.include

if EXAFMM_HAVE_MPI
lib_LIBRARIES = libpetiga_exafmm_laplace.a libpetiga_exafmm_helmholtz.a libgromacs_exafmm.a libcharmm_exafmm.a libmd_distributed.a
libpetiga_exafmm_laplace_a_SOURCES = petiga.cxx
libpetiga_exafmm_laplace_a_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE
libpetiga_exafmm_helmholtz_a_SOURCES = petiga.cxx
libpetiga_exafmm_helmholtz_a_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_HELMHOLTZ
libgromacs_exafmm_a_SOURCES = gromacs.cxx
libgromacs_exafmm_a_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE
libcharmm_exafmm_a_SOURCES = charmm.cxx
libcharmm_exafmm_a_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE
libmd_distributed_a_SOURCES = md_distributed.cxx
libmd_distributed_a_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE

bin_PROGRAMS = petiga_laplace_mpi petiga_helmholtz_mpi gromacs_mpi charmm_mpi md_distributed_mpi
petiga_laplace_mpi_SOURCES = test_petiga.cxx
petiga_laplace_mpi_CPPFLAGS = $(libpetiga_laplace_a_CPPFLAGS)
petiga_laplace_mpi_LDADD = libpetiga_exafmm_laplace.a
petiga_helmholtz_mpi_SOURCES = test_petiga.cxx
petiga_helmholtz_mpi_CPPFLAGS = $(libpetiga_helmholtz_a_CPPFLAGS)
petiga_helmholtz_mpi_LDADD = libpetiga_exafmm_helmholtz.a
gromacs_mpi_SOURCES = test_gromacs.cxx
gromacs_mpi_CPPFLAGS = $(libgromacs_exafmm_a_CPPFLAGS)
gromacs_mpi_LDADD = libgromacs_exafmm.a
charmm_mpi_SOURCES = test_charmm.f90
charmm_mpi_CPPFLAGS = $(libcharmm_exafmm_a_CPPFLAGS)
if EXAFMM_HAVE_FX
nodist_EXTRA_charmm_mpi_SOURCES = test_gromacs.cxx
charmm_mpi_LINK = $(CXXLINK) --linkfortran
charmm_mpi_LDADD = libcharmm_exafmm.a
else
charmm_mpi_LINK = $(FCLINK) $(AM_FCFLAGS)
charmm_mpi_LDADD = libcharmm_exafmm.a -lstdc++ -lgomp
endif
md_distributed_mpi_SOURCES = test_md_distributed.cxx
md_distributed_mpi_CPPFLAGS = $(libmd_distributed_a_CPPFLAGS)
md_distributed_mpi_LDADD = libmd_distributed.a

run_charmm_mpi: charmm_mpi
	$(MPIRUN) -n 2 ./$< ./ water_32.fmm water.out 0
run_gromacs_mpi: gromacs_mpi
	$(MPIRUN) -n 2 ./$<
run_petiga_laplace_mpi: petiga_laplace_mpi
	$(MPIRUN) -n 2 ./$<
run_petiga_helmholtz_mpi: petiga_helmholtz_mpi
	$(MPIRUN) -n 1 ./$<
run_md_distributed_mpi: md_distributed_mpi
	$(MPIRUN) -n 2 ./$<
endif
