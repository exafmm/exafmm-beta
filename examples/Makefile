include ../Makefile.include

#LFLAGS	+= -DMANY

kernel	: kernel.cxx $(KERNELS)
	$(CXX) $? $(LFLAGS)
	./a.out

serial	: serial.cxx $(KERNELS)
	$(CXX) $? $(LFLAGS)
	./a.out
#	./a.out --numBodies 100 --mutual 0
#	numactl --interleave=all ./a.out

parallel: parallel.cxx $(KERNELS)
	$(CXX) $? $(LFLAGS)
	mpirun -np 8 ./a.out

wrapper: wrapper.cxx $(KERNELS)
	make -C ../wrappers libfmm.a
	$(CXX) $? -L../wrappers -lfmm $(LFLAGS)
	mpirun -np 4 ./a.out
