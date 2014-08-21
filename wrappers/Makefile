include ../Makefile.include

#LFLAGS += -DMASS
LFLAGS += -DIneJ

all:
	make charmm
	make gromacs
	make petiga

petiga:
	make libpetiga.a
	$(CXX) test_petiga.cxx -L. -lfmm -ltbb
	mpirun -np 2 ./a.out

gromacs:
	make libgromacs.a
	$(CXX) test_gromacs.cxx -L. -lfmm -ltbb
	mpirun -np 2 ./a.out

charmm:
	make libcharmm.a
	$(FC) test_charmm.f90 -ggdb3 -Wall -fbounds-check -L. -lfmm -ltbb -lstdc++ -lgomp
#	mpirun -np 2 ./a.out ala3.fmm
	make cleandat
	mpirun -np 2 ./a.out water.fmm 10
	make cleandat
	mpirun -np 1 ./a.out water.fmm 10

libpetiga.a: petiga.o ../kernels/Laplace$(BASIS)$(DEVICE).o ../kernels/LaplaceP2P$(DEVICE).o
	ar ruv libpetiga.a $?
	ranlib libpetiga.a
	mv libpetiga.a libfmm.a

libgromacs.a: gromacs.o ../kernels/Laplace$(BASIS)$(DEVICE).o ../kernels/LaplaceP2P$(DEVICE).o
	ar ruv libgromacs.a $?
	ranlib libgromacs.a
	mv libgromacs.a libfmm.a

libcharmm.a: charmm.o ../kernels/Laplace$(BASIS)$(DEVICE).o ../kernels/LaplaceP2P$(DEVICE).o
	ar ruv libcharmm.a $?
	ranlib libcharmm.a
	mv libcharmm.a libfmm.a
