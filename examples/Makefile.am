include ../Makefile.am.include

bin_PROGRAMS = laplace_cartesian laplace_spherical helmholtz_spherical tree hss
laplace_cartesian_SOURCES = fmm.cxx ../kernels/LaplaceCartesianCPU.cxx ../kernels/LaplaceP2PCPU.cxx
laplace_cartesian_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE -DEXAFMM_CARTESIAN -DEXAFMM_MASS -DEXAFMM_EXPANSION=4
laplace_spherical_SOURCES = fmm.cxx ../kernels/LaplaceSphericalCPU.cxx ../kernels/LaplaceP2PCPU.cxx
laplace_spherical_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_LAPLACE -DEXAFMM_SPHERICAL -DEXAFMM_EXPANSION=10
helmholtz_spherical_SOURCES = fmm.cxx ../kernels/HelmholtzSphericalCPU.cxx ../kernels/HelmholtzP2PCPU.cxx
helmholtz_spherical_CPPFLAGS = $(AM_CPPFLAGS) -DEXAFMM_HELMHOLTZ -DEXAFMM_SPHERICAL -DEXAFMM_EXPANSION=10
tree_SOURCES = tree.cxx
tree_CPPFLAGS = $(laplace_cartesian_CPPFLAGS)
hss_SOURCES = hss.cxx ../kernels/LaplaceCartesianCPU.cxx ../kernels/LaplaceP2PCPU.cxx
hss_CPPFLAGS = $(AM_CPPFLAGS) -I ~/strumpack/src -DEXAFMM_LAPLACE -DEXAFMM_CARTESIAN -DEXAFMM_EXPANSION=4
hss_LDADD = -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5

run_laplace_cartesian: laplace_cartesian
	./$< -Dgmovx
run_laplace_spherical: laplace_spherical
	./$< -Dgmovx
run_helmholtz_spherical: helmholtz_spherical
	./$< -Dgv -n 100000
run_tree: tree
	./$< -v -n 10000000
run_hss: hss
	./$< -DgGmovx -n 4096

bin_PROGRAMS += laplace_cartesian_kernel laplace_spherical_kernel helmholtz_spherical_kernel
laplace_cartesian_kernel_SOURCES = kernel.cxx ../kernels/LaplaceCartesianCPU.cxx ../kernels/LaplaceP2PCPU.cxx
laplace_cartesian_kernel_CPPFLAGS = $(laplace_cartesian_CPPFLAGS)
laplace_spherical_kernel_SOURCES = kernel.cxx ../kernels/LaplaceSphericalCPU.cxx ../kernels/LaplaceP2PCPU.cxx
laplace_spherical_kernel_CPPFLAGS = $(laplace_spherical_CPPFLAGS)
helmholtz_spherical_kernel_SOURCES = kernel.cxx ../kernels/HelmholtzSphericalCPU.cxx ../kernels/HelmholtzP2PCPU.cxx
helmholtz_spherical_kernel_CPPFLAGS = $(helmholtz_spherical_CPPFLAGS)

run_laplace_cartesian_kernel: laplace_cartesian_kernel
	./$< -v
run_laplace_spherical_kernel: laplace_spherical_kernel
	./$< -v
run_helmholtz_spherical_kernel: helmholtz_spherical_kernel
	./$< -v

if EXAFMM_HAVE_MPI
bin_PROGRAMS += laplace_cartesian_mpi laplace_spherical_mpi helmholtz_spherical_mpi ewald_mpi
laplace_cartesian_mpi_SOURCES = fmm_mpi.cxx ../kernels/LaplaceCartesianCPU.cxx ../kernels/LaplaceP2PCPU.cxx
laplace_cartesian_mpi_CPPFLAGS = $(laplace_cartesian_CPPFLAGS)
laplace_spherical_mpi_SOURCES = fmm_mpi.cxx ../kernels/LaplaceSphericalCPU.cxx ../kernels/LaplaceP2PCPU.cxx
laplace_spherical_mpi_CPPFLAGS = $(laplace_spherical_CPPFLAGS)
helmholtz_spherical_mpi_SOURCES = fmm_mpi.cxx ../kernels/HelmholtzSphericalCPU.cxx ../kernels/HelmholtzP2PCPU.cxx
helmholtz_spherical_mpi_CPPFLAGS = $(helmholtz_spherical_CPPFLAGS)
ewald_mpi_SOURCES = ewald.cxx ../kernels/LaplaceSphericalCPU.cxx ../kernels/LaplaceP2PCPU.cxx
ewald_mpi_CPPFLAGS = $(laplace_spherical_CPPFLAGS)

run_laplace_cartesian_mpi: laplace_cartesian_mpi
	$(MPIRUN) -n 2 ./$< -Dgmovx
run_laplace_spherical_mpi: laplace_spherical_mpi
	$(MPIRUN) -n 2 ./$< -Dgmovx
run_helmholtz_spherical_mpi: helmholtz_spherical_mpi
	$(MPIRUN) -n 2 ./$< -Dgv -n 100000
run_ewald_mpi: ewald_mpi
	$(MPIRUN) -n 2 ./$< -Dgmovx
endif
