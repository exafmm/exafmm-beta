include ../Makefile.include

#LFLAGS	+= -DMASS
LFLAGS 	+= -DIneJ
LFLAGS 	+= -DNOWRAP
#LFLAGS	+= -DNO_P2P
LFLAGS	+= -DCluster
#LFLAGS	+= -DSerial

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
	$(CXX) test_gromacs.cxx $(LFLAGS) -lfmm
	mpirun -np 1 ./a.out

charmm: # exafmm
	make libcharmm.a
	$(FC) test_charmm.f90 -L. -lfmm -ltbb -lstdc++ -lgomp
	make cleandat
	mpirun -np 1 ./a.out water_32.fmm water.out 1
	make cleandat
	mpirun -np 2 ./a.out water_32.fmm water.out 1

charmm2: # uniform
	make libcharmm2.a
	$(FC) test_charmm.f90 -L. -lfmm -ltbb -lstdc++ -lgomp
	make cleandat
	mpirun -np 2 ./a.out water_32.fmm water.out 1

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

libcharmm2.a: charmm2.o ../kernels/Laplace$(BASIS)$(DEVICE).o ../kernels/LaplaceP2P$(DEVICE).o
	ar ruv libcharmm2.a $?
	ranlib libcharmm2.a
	mv libcharmm2.a libfmm.a
