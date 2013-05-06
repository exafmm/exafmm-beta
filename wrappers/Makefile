include ../Makefile.include

laplace: test_laplace.o ../kernels/Laplace$(BASIS)$(DEVICE).o ../kernels/LaplaceP2P$(DEVICE).o
	make liblaplace.a
	$(CXX) $? -L. -llaplace $(LFLAGS)
	mpirun -np 4 ./a.out

liblaplace.a: laplace.o
	ar ruv liblaplace.a $?
	ranlib liblaplace.a

liblaplace.so: laplace.o
	$(CXX) $? -shared -o liblaplace.so $(LFLAGS)

matern: test_matern.o ../kernels/Matern$(BASIS)$(DEVICE).o ../kernels/MaternP2P$(DEVICE).o
	make libmatern.a
	$(CXX) $? -L. -lmatern $(LFLAGS)
	mpirun -np 4 ./a.out

libmatern.a: matern.o
	ar ruv libmatern.a $?
	ranlib libmatern.a

libmatern.so: matern.o
	$(CXX) $? -shared -o libmatern.so $(LFLAGS)
