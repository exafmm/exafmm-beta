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
	$(FC) test_charmm.f90 -L. -lfmm -ltbb -lstdc++
	mpirun -np 2 ./a.out ala3.fmm

libgromacs.a: gromacs.o ../kernels/Laplace$(BASIS)$(DEVICE).o ../kernels/LaplaceP2P$(DEVICE).o
	ar ruv libgromacs.a $?
	ranlib libgromacs.a
	mv libgromacs.a libfmm.a

libcharmm.a: charmm.o ../kernels/Laplace$(BASIS)$(DEVICE).o ../kernels/LaplaceP2P$(DEVICE).o
	ar ruv libcharmm.a $?
	ranlib libcharmm.a
	mv libcharmm.a libfmm.a