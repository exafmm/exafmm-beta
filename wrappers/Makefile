include ../Makefile.include

all:
	make gromacs
	make charmm

gromacs:
	make libgromacs.a
	$(CXX) test_gromacs.cxx -L. -lfmm -ltbb
	mpirun -np 2 ./a.out

charmm:
	make libcharmm.a
	$(FC) test_charmm.f90 -ggdb3 -Wall -fbounds-check -L. -lfmm -ltbb -lmpi_cxx -lstdc++
#	mpirun -np 2 ./a.out ala3.fmm
	mpirun -np 2 ./a.out water.fmm 10
	make cleandat
	mpirun -np 1 ./a.out water.fmm 10

libgromacs.a: gromacs.o ../kernels/Laplace$(BASIS)$(DEVICE).o ../kernels/LaplaceP2P$(DEVICE).o
	ar ruv libgromacs.a $?
	ranlib libgromacs.a
	mv libgromacs.a libfmm.a

libcharmm.a: charmm.o ../kernels/Laplace$(BASIS)$(DEVICE).o ../kernels/LaplaceP2P$(DEVICE).o
	ar ruv libcharmm.a $?
	ranlib libcharmm.a
	mv libcharmm.a libfmm.a
