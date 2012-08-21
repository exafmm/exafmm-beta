include ../Makefile.include

#LFLAGS	+= -DMANY

serial	: serial.cxx $(KERNELS)
	$(CXX) $? $(LFLAGS)
	./a.out

serial_argv : serial_argv.cxx $(KERNELS)
	$(CXX) $? $(LFLAGS) $(VFLAGS)
	./a.out

parallel: parallel.cxx $(KERNELS)
	$(CXX) $? $(LFLAGS) $(VFLAGS)
	mpirun -np 32 ./a.out

wrapper: wrapper.cxx $(KERNELS)
	make -C ../wrapper libcoulomb.a
	$(CXX) $? $(LFLAGS) -lcoulomb
	mpirun -np 4 ./a.out
