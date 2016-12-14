include ../Makefile.am.include

bin_PROGRAMS = fmm tree
fmm_SOURCES = fmm.cxx
if EXAFMM_HAVE_FX
fmm_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_PMAX=7
else
bin_PROGRAMS += vec
fmm_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_PMAX=10
endif
tree_SOURCES = tree.cxx
tree_CPPFLAGS = $(fmm_CPPFLAGS)
vec_SOURCES = vec.cxx
vec_CPPFLAGS = $(fmm_CPPFLAGS)

run_laplace_cartesian: fmm
	./$< -DgmMovx -r 10 -e laplace -b cartesian -P 10
run_laplace_spherical: fmm
	./$< -Dgmovx -r 10 -e laplace -P 10
run_helmholtz_spherical: fmm
	./$< -Dgv -n 40000 -r 10 -e helmholtz -P 14
run_biotsavart_spherical: fmm
	./$< -Dgv -n 40000 -r 10 -e biotsavart -P 10
run_tree: tree
	./$< -v -n 10000000 -r 10
run_vec: vec
	./$<

bin_PROGRAMS += kernel
kernel_SOURCES = kernel.cxx
kernel_CPPFLAGS = $(fmm_CPPFLAGS)

run_laplace_cartesian_kernel: kernel
	./$< -v -e laplace -b cartesian -P 10
run_laplace_spherical_kernel: kernel
	./$< -v -e laplace -P 10
run_helmholtz_spherical_kernel: kernel
	./$< -v -e helmholtz -P 14
run_biotsavart_spherical_kernel: kernel
	./$< -v -e biotsavart -P 10

if EXAFMM_HAVE_MPI
bin_PROGRAMS += fmm_mpi ewald_mpi
fmm_mpi_SOURCES = fmm_mpi.cxx
fmm_mpi_CPPFLAGS = $(fmm_CPPFLAGS)
ewald_mpi_SOURCES = ewald.cxx
ewald_mpi_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_PMAX=10

run_laplace_cartesian_mpi: fmm_mpi
	$(MPIRUN) -n 2 ./$< -DgmMovx -r 10 -e laplace -b cartesian -P 10
run_laplace_spherical_mpi: fmm_mpi
	$(MPIRUN) -n 2 ./$< -Dgmovx -r 10 -e laplace -P 10
run_helmholtz_spherical_mpi: fmm_mpi
	$(MPIRUN) -n 2 ./$< -Dgv -n 100000 -r 10 -e helmholtz -P 14
run_biotsavart_spherical_mpi: fmm_mpi
	$(MPIRUN) -n 2 ./$< -Dgv -n 100000 -r 10 -e biotsavart -P 10
run_ewald_mpi: ewald_mpi
	$(MPIRUN) -n 2 ./$< -Dgmovx -r 10 -P 10
endif
