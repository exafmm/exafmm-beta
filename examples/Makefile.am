include ../Makefile.am.include

bin_PROGRAMS = fmm tree vec
fmm_SOURCES = fmm.cxx
fmm_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_PMAX=10
tree_SOURCES = tree.cxx
tree_CPPFLAGS = $(fmm_CPPFLAGS)
vec_SOURCES = vec.cxx
vec_CPPFLAGS = $(fmm_CPPFLAGS)

run_laplace: fmm
	./$< -aDgmovx -n 100000 -r 10 -e laplace -P 10
run_helmholtz: fmm
	./$< -aDgv -n 10000 -r 10 -e helmholtz -P 14
run_biotsavart: fmm
	./$< -aDgv -n 10000 -r 10 -e biotsavart -P 10
run_tree: tree
	./$< -v -n 10000000 -r 10
run_vec: vec
	./$<

bin_PROGRAMS += kernel
kernel_SOURCES = kernel.cxx
kernel_CPPFLAGS = $(fmm_CPPFLAGS)

run_laplace_kernel: kernel
	./$< -v -e laplace -P 10
run_helmholtz_kernel: kernel
	./$< -v -e helmholtz -P 14
run_biotsavart_kernel: kernel
	./$< -v -e biotsavart -P 10

if EXAFMM_HAVE_MPI
bin_PROGRAMS += fmm_mpi ewald_mpi
fmm_mpi_SOURCES = fmm_mpi.cxx
fmm_mpi_CPPFLAGS = $(fmm_CPPFLAGS)
ewald_mpi_SOURCES = ewald.cxx
ewald_mpi_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_PMAX=10

run_laplace_mpi: fmm_mpi
	$(MPIRUN) -n 2 ./$< -aDgmovx -n 100000 -r 10 -e laplace -P 10
run_helmholtz_mpi: fmm_mpi
	$(MPIRUN) -n 2 ./$< -aDgv -n 10000 -r 10 -e helmholtz -P 14
run_biotsavart_mpi: fmm_mpi
	$(MPIRUN) -n 2 ./$< -aDgv -n 10000 -r 10 -e biotsavart -P 10
run_ewald_mpi: ewald_mpi
	$(MPIRUN) -n 2 ./$< -aDgmovx -r 10 -P 10
endif
