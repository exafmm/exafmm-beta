include ../Makefile.include

#LFLAGS	+= -DMANY

serial: serial.o $(OBJECTS)
	$(CXX) $? $(LFLAGS)
#	./a.out
	./a.out --theta 0.5 --numBodies 100
#	numactl --interleave=all ./a.out

parallel: parallel.o $(OBJECTS)
	$(CXX) $? $(LFLAGS)
	mpirun -np 8 ./a.out

wrapper: wrapper.o $(OBJECTS)
	make -C ../wrappers libfmm.a
	$(CXX) $? -L../wrappers -lfmm $(LFLAGS)
	mpirun -np 4 ./a.out

kernel: kernel.o $(OBJECTS)
	$(CXX) $? $(LFLAGS)
	./a.out

convergence: kernel.cxx $(SOURCES)
	rm -f kernel.dat
	for P in 3 4 5 6 7 8; do echo P = $$P compiling && \
	$(CXX) $? $(LFLAGS) -DEXPANSION=$$P && \
	echo P = $$P executing && ./a.out; done