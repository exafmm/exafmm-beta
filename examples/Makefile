include ../Makefile.include

LFLAGS	+= -DMASS # Use all positive sources
#LFLAGS	+= -DIneJ # Different target and source points
LFLAGS	+= -DDIRECT # Compare with direct summation
#LFLAGS	+= -DWRITE_TIME # Write timings to file

# Test for tree construction
tree: tree.o $(OBJECTS)
	$(CXX) $? $(LFLAGS)
	rm -f time.dat
	for N in 10000 20000 50000 100000 200000 500000 1000000 2000000 5000000 10000000 20000000 50000000 100000000; do \
	  for T in 1 2 4 8 16 32 64; do \
	    echo N = $$N && ./a.out -n $$N -T $$T -c 32 -r 10 -d cube; \
	  done; \
	done
	mv time.dat stelletto-18-10-exafmm-cube.txt
	for N in 10000 20000 50000 100000 200000 500000 1000000 2000000 5000000 10000000 20000000 50000000 100000000; do \
	  for T in 1 2 4 8 16 32 64; do \
	    echo N = $$N && ./a.out -n $$N -T $$T -c 32 -r 10 -d octant; \
	  done; \
	done
	mv time.dat stelletto-18-10-exafmm-octant.txt
	for N in 10000 20000 50000 100000 200000 500000 1000000 2000000 5000000 10000000 20000000 50000000 100000000; do \
	  for T in 1 2 4 8 16 32 64; do \
	    echo N = $$N && ./a.out -n $$N -T $$T -c 32 -r 10 -d plummer; \
	  done; \
	done
	mv time.dat stelletto-18-10-exafmm-plummer.txt

# Test for serial FMM
serial: serial.o $(OBJECTS)
	$(CXX) $? $(LFLAGS)
	./a.out

Nserial:
	@for N in 1000 1333 1778 2371 3162 4216 5623 7498 10000 13335 17782 23713 31622 42169 56234 74989 100000 \
	133352 177827 237137 316227 421696 562341 749894 1000000 1333521 1778279 2371373 3162277 4216965 5623413 7498942 10000000; do \
	./a.$(P) -n $$N -v 0 -T $(T) -t $(t) -c $(c); echo; \
	done

# Test for parallel FMM
parallel: parallel.o $(OBJECTS)
	$(CXX) $? $(LFLAGS)
	mpirun -np 2 ./a.out
#	mpirun -np 9 ./a.out --ncrit 256 --distribution plummer

# Checking O(N) complexity
complexity: serial.o $(OBJECTS)
	$(CXX) $? $(LFLAGS)
	for N in 10000 20000 50000 100000 200000 500000 1000000; do \
        ./a.out --numBodies $$N; echo; \
	done

# Test for kernels only
kernel: kernel.o $(OBJECTS)
	$(CXX) $? $(LFLAGS)
	./a.out

# Checking O(theta^p) convergence
convergence: kernel.cxx $(SOURCES)
	rm -f kernel.dat
	for P in 3 4 5 6 7 8 9 10 11 12 13 14 15; do echo P = $$P compiling && \
	$(CXX) $? $(LFLAGS) -DEXPANSION=$$P && \
	echo P = $$P executing && ./a.out; done

# Ewald vs. periodic FMM
ewald: ewald.o $(OBJECTS)
	$(CXX) $? $(LFLAGS)
	mpirun -np 2 ./a.out --numBodies 1000 --ncrit 32 --images 3
