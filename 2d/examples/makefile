include ../Makefile.include

#LFLAGS	+= -DMASS # Use all positive sources
LFLAGS	+= -DIneJ # Different target and source points

# Test for serial FMM
serial: serial.o
	$(CXX) $? $(LFLAGS)
	./a.out

# Test for parallel FMM
parallel: parallel.o
	$(CXX) $? $(LFLAGS)
	mpirun -np 8 ./a.out --numBodies 10000000

# Test for hss FMM
hss: hss.o
	$(CXX) $? $(LFLAGS)
	mpirun -n 4 ./a.out --numBodies 4096

# Checking O(N) complexity
complexity: serial.o
	$(CXX) $? $(LFLAGS)
	for N in 10000 20000 50000 100000 200000 500000 1000000; do \
        ./a.out --numBodies $$N --verbose 0; echo; \
	done

# Test for kernels only
kernel: kernel.o
	$(CXX) $? $(LFLAGS)
	./a.out

# Checking O(theta^p) convergence
convergence: kernel.cxx
	rm -f kernel.dat
	for P in 3 4 5 6 7 8 9 10; do echo P = $$P compiling && \
	$(CXX) $? $(LFLAGS) -DEXPANSION=$$P && \
	echo P = $$P executing && ./a.out; done
