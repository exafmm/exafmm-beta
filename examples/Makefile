include ../Makefile.include

#LFLAGS	+= -DMANY

serial	: serial.cxx $(KERNELS)
	$(CXX) $? $(LFLAGS) -DEVAL_ERROR_PARTIAL_ACCUMULATE
	mpirun -np 1 ./a.out --numBodies 10000
#	numactl --interleave=all ./a.out

parallel: parallel.cxx $(KERNELS)
	$(CXX) $? $(LFLAGS)
	mpirun -np 8 ./a.out

wrapper: wrapper.cxx $(KERNELS)
	make -C ../wrappers libcoulomb.a
	$(CXX) $? -lcoulomb $(LFLAGS)
	mpirun -np 4 ./a.out
