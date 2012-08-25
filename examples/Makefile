include ../Makefile.include

#LFLAGS	+= -DMANY

serial	: serial.cxx $(KERNELS)
	$(CXX) $? $(LFLAGS)
	./a.out

parallel: parallel.cxx $(KERNELS)
	$(CXX) $? $(LFLAGS)
	mpirun -np 32 ./a.out

wrapper: wrapper.cxx $(KERNELS)
	make -C ../wrappers libcoulomb.a
	$(CXX) $? $(LFLAGS) -lcoulomb
	mpirun -np 4 ./a.out
