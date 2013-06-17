include ../Makefile.include

all:
	make gromacs
	make charmm

gromacs:
	make libgromacs.a
	make libewald.a
	$(CXX) test_gromacs.cxx -L. -lfmm -ltbb -lewald
	mpirun -np 2 ./a.out

charmm:
	make libcharmm.a
	make libewald.a
	$(FC) test_charmm.f90 -L. -lfmm -ltbb -lewald -lstdc++
	mpirun -np 2 ./a.out

libgromacs.a: gromacs.o ../kernels/Laplace$(BASIS)$(DEVICE).o ../kernels/LaplaceP2P$(DEVICE).o
	ar ruv libgromacs.a $?
	ranlib libgromacs.a
	mv libgromacs.a libfmm.a

libcharmm.a: charmm.o ../kernels/Laplace$(BASIS)$(DEVICE).o ../kernels/LaplaceP2P$(DEVICE).o
	ar ruv libcharmm.a $?
	ranlib libcharmm.a
	mv libcharmm.a libfmm.a

libewald.a: ewald.o
	ar ruv libewald.a $?
	ranlib libewald.a