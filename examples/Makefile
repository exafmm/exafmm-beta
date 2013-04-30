include ../Makefile.include

#LFLAGS	+= -DMANY

serial: serial.cxx $(OBJECTS)
	$(CXX) $? $(LFLAGS)
#	./a.out
	./a.out --numBodies 100 --mutual 0
#	numactl --interleave=all ./a.out

parallel: parallel.cxx $(OBJECTS)
	$(CXX) $? $(LFLAGS)
	mpirun -np 8 ./a.out

wrapper: wrapper.cxx $(OBJECTS)
	make -C ../wrappers libfmm.a
	$(CXX) $? -L../wrappers -lfmm $(LFLAGS)
	mpirun -np 4 ./a.out

kernel: kernel.cxx $(OBJECTS)
	$(CXX) $? $(LFLAGS)
	./a.out

convergence: kernel.cxx $(SOURCES)
	rm -f kernel.dat
	for P in 3 4 5 6 7 8; do echo P = $$P compiling && \
	$(CXX) $? $(LFLAGS) -DEXPANSION=$$P && \
	echo P = $$P executing && ./a.out; done