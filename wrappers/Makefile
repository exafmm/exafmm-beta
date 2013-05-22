include ../Makefile.include
FC = mpif90 -ggdb3 -O3 -msse4a -ffast-math -funroll-loops

laplace:
	make liblaplace.a
	$(CXX) test_laplace.cxx -L. -lfmm -ltbb
	mpirun -np 8 ./a.out
	$(FC) test_laplace.f90 -L. -lfmm -ltbb -lstdc++ -lmpi_cxx
	mpirun -np 8 ./a.out

liblaplace.a: laplace.o ../kernels/Laplace$(BASIS)$(DEVICE).o ../kernels/LaplaceP2P$(DEVICE).o
	ar ruv liblaplace.a $?
	ranlib liblaplace.a
	mv liblaplace.a libfmm.a

liblaplace.so: laplace.o
	$(CXX) $? -shared -o liblaplace.so $(LFLAGS)
	mv liblaplace.so libfmm.so
