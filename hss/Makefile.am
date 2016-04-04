include ../Makefile.am.include

if EXAFMM_HAVE_MPI
if EXAFMM_WITH_STRUMPACK
bin_PROGRAMS = laplace helmholtz
laplace_SOURCES = laplace.cxx
laplace_CPPFLAGS = $(AM_CPPFLAGS) -I ~/strumpack/src -DEXAFMM_LAPLACE -DEXAFMM_CARTESIAN -DEXAFMM_EXPANSION=4
laplace_LDADD = -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5
helmholtz_SOURCES = helmholtz.cxx
helmholtz_CPPFLAGS = $(AM_CPPFLAGS) -I ~/strumpack/src -DEXAFMM_HELMHOLTZ -DEXAFMM_SPHERICAL -DEXAFMM_EXPANSION=10
helmholtz_LDADD = -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5

run_laplace: laplace
	$(MPIRUN) -n 4 ./$< -DgGmox -n 4096
run_helmholtz: helmholtz
	$(MPIRUN) -n 4 ./$< -DgG -n 4096
endif
endif
