include ../Makefile.include

#LFLAGS	+= -DMANY

serial	: serial.cxx $(KERNELS)
	$(CXX) $? $(LFLAGS) -DEVAL_ERROR_PARTIAL_ACCUMULATE
#	./a.out --numBodies 20
	./a.out --numBodies 100000 --theta 0.5
#	numactl --interleave=all ./a.out

parallel: parallel.cxx $(KERNELS)
	$(CXX) $? $(LFLAGS)
	mpirun -np 8 ./a.out

wrapper: wrapper.cxx $(KERNELS)
	make -C ../wrappers libcoulomb.a
	$(CXX) $? -lcoulomb $(LFLAGS)
	mpirun -np 4 ./a.out
