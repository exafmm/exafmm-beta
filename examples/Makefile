include ../Makefile.include

#LFLAGS	+= -DMASS # Use all positive sources
#LFLAGS	+= -DIneJ # Different target and source points

# Test for serial FMM
serial: serial.o $(OBJECTS)
	$(CXX) $? $(LFLAGS)
	./a.out

# Test for parallel FMM
parallel: parallel.o $(OBJECTS)
	$(CXX) $? $(LFLAGS)
	mpirun -np 8 ./a.out

# Checking O(N) complexity
complexity: serial.o $(OBJECTS)
	$(CXX) $? $(LFLAGS)
	for N in 10000 20000 50000 100000 200000 500000 1000000; do \
        ./a.out --numBodies $$N --verbose 0; echo; \
	done

# Test for kernels only
kernel: kernel.o $(OBJECTS)
	$(CXX) $? $(LFLAGS)
	./a.out

# Checking O(theta^p) convergence (Turn on VTK to see a graph at the end)
convergence: kernel.cxx $(SOURCES)
	rm -f kernel.dat
	for P in 3 4 5 6 7 8 9 10; do echo P = $$P compiling && \
	$(CXX) $? $(LFLAGS) -DEXPANSION=$$P && \
	echo P = $$P executing && ./a.out; done

# Ewald vs. periodic FMM
ewald: ewald.o $(OBJECTS)
	$(CXX) $? $(LFLAGS)
	mpirun -np 1 ./a.out --numBodies 1000 --theta .35 --ncrit 32 --images 3
