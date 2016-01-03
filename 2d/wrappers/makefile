include ../Makefile.include

# Test for Octave wrapper
octave:
	mkoctfile -v --mex -I../include -DFP64 -DEXPANSION=20 exafmm2d.cpp

# Test for Matlab wrapper
matlab:
	mex -I../include -DFP64 -DEXPANSION=20 exafmm2d.cpp

petiga:
	make libpetiga.a
	$(CXX) test_petiga.cxx -L. -lfmm -ltbb
	mpirun -np 8 ./a.out

nek5000:
	make libnek5000.a
	$(FC) test_nek5000.f90 -L. -lfmm -ltbb -lmpi_cxx -lstdc++ -lgomp
	mpirun -np 8 ./a.out

libpetiga.a: petiga.o
	ar ruv libpetiga.a $?
	ranlib libpetiga.a
	mv libpetiga.a libfmm.a

libnek5000.a: nek5000.o
	ar ruv libnek5000.a $?
	ranlib libnek5000.a
	mv libnek5000.a libfmm.a
