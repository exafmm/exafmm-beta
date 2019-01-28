include ../Makefile.am.include

if EXAFMM_HAVE_MPI
if EXAFMM_WITH_STRUMPACK
bin_PROGRAMS = laplace helmholtz
laplace_SOURCES = laplace.cxx
laplace_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE -I ~/strumpack/src
laplace_LDADD = -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
helmholtz_SOURCES = helmholtz.cxx
helmholtz_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_HELMHOLTZ -I ~/strumpack/src
helmholtz_LDADD = -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5

run_laplace: laplace
	$(MPIRUN) -n 4 ./$< -DgG -n 1024 -P 4
run_helmholtz: helmholtz
	$(MPIRUN) -n 4 ./$< -DgG -n 256 -P 10
endif
endif
