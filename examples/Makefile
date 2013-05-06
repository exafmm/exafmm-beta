include ../Makefile.include

#LFLAGS	+= -DMANY

serial: serial.o $(OBJECTS)
	$(CXX) $? $(LFLAGS)
	./a.out --numBodies 10000
#	numactl --interleave=all ./a.out

parallel: parallel.o $(OBJECTS)
	$(CXX) $? $(LFLAGS)
	mpirun -np 8 ./a.out

kernel: kernel.o $(OBJECTS)
	$(CXX) $? $(LFLAGS)
	./a.out

convergence: kernel.cxx $(SOURCES)
	rm -f kernel.dat
	for P in 3 4 5 6 7 8; do echo P = $$P compiling && \
	$(CXX) $? $(LFLAGS) -DEXPANSION=$$P && \
	echo P = $$P executing && ./a.out; done