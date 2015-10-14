include ../Makefile.am.include

if MPI
bin_PROGRAMS = fmm fmm_mpi fmm_ijhpca
fmm_SOURCES = fmm.cxx
fmm_CPPFLAGS = $(AM_CPPFLAGS) -DLaplace -DCartesian -DEXPANSION=6 -DSerial
fmm_mpi_SOURCES = fmm.cxx
fmm_mpi_CPPFLAGS = $(AM_CPPFLAGS) -DLaplace -DCartesian -DEXPANSION=6
fmm_ijhpca_SOURCES = fmm.cxx
fmm_ijhpca_CPPFLAGS = $(AM_CPPFLAGS) -DLaplace -DCartesian -DEXPANSION=6 -DIJHPCA

run_fmm:
	./fmm -n 1000 -c 32 -i 3
run_fmm_mpi:
	mpirun -np 8 ./fmm_mpi -n 1000 -c 32 -i 3
run_fmm_ijhpca:
	mpirun -np 8 ./fmm_ijhpca -n 62500 -c 32 -i 1
endif
