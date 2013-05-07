include ../Makefile.include

#LFLAGS	+= -DIneJ # Different target and source points

serial: serial.o $(OBJECTS)
	$(CXX) $? $(LFLAGS)
	./a.out

parallel: parallel.o $(OBJECTS)
	$(CXX) $? $(LFLAGS)
	mpirun -np 8 ./a.out

scalability: serial.o $(OBJECTS)
	$(CXX) $? $(LFLAGS)
	for N in 10000 20000 50000 100000 200000 500000 1000000; do \
        ./a.out --numBodies $$N --verbose 0; echo; \
	done

kernel: kernel.o $(OBJECTS)
	$(CXX) $? $(LFLAGS)
	./a.out

convergence: kernel.cxx $(SOURCES)
	rm -f kernel.dat
	for P in 3 4 5 6 7 8; do echo P = $$P compiling && \
	$(CXX) $? $(LFLAGS) -DEXPANSION=$$P && \
	echo P = $$P executing && ./a.out; done