include ../Makefile.include

#LFLAGS	+= -DMASS
LFLAGS 	+= -DIneJ
LFLAGS 	+= -DNOWRAP
#LFLAGS	+= -DNO_P2P
LFLAGS	+= -DSerial

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
	mpirun -np 1 ./a.out

charmm: # exafmm-dev
	make libcharmm.a
	$(FC) test_charmm.f90 -L. -lfmm -ltbb -lstdc++ -lgomp
	mpirun -np 1 ./a.out ala3.fmm
	make cleandat
	mpirun -np 1 ./a.out water.fmm water.out 10
	make cleandat
	mpirun -np 2 ./a.out water.fmm water.out 10

charmm2: # uniform
	make libcharmm2.a
	$(FC) test_charmm.f90 -L. -lfmm -ltbb -lstdc++ -lgomp
	make cleandat
	OMP_NUM_THREADS=1 mpirun -np 1 ./a.out water.fmm water.out 1

charmm2spike: # uniform
	rm -f spike* uniform.data
	cp water.fmm spike.0
	touch uniform.data
	make libcharmm2.a
	$(FC) test_charmm.f90 -L. -lfmm -ltbb -lstdc++ -lgomp -DSPIKE
	it=0; while [[ $$it -lt 100000 ]]; do \
	((it1 = it + 1)); \
	make cleandat; \
	OMP_NUM_THREADS=1 mpirun -np 1 ./a.out spike.$$it spike.$$it1 1; \
        ((it = it + 1)); done

charmm3: # uniform-serial (Sakura interface)
	make libcharmm3.a
	$(FC) test_charmm.f90 -L. -lfmm -ltbb -lstdc++ -lgomp
	make cleandat
	mpirun -np 1 ./a.out water.fmm water.out 10

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

libcharmm3.a: charmm3.o
	ar ruv libcharmm3.a $?
	ranlib libcharmm3.a
	mv libcharmm3.a libfmm.a
