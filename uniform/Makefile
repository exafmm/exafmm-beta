include ../Makefile.include

serial	: main.cxx $(OBJECTS)
	$(CXX) $? $(LFLAGS) -DSerial
	OMP_NUM_THREADS=1 ./a.out -n 1000 -c 32 -i 3

parallel: main.cxx $(OBJECTS)
	$(CXX) $? $(LFLAGS)
	mpirun -np 8 ./a.out -n 1000 -c 32 -i 3
#	$(CXX) $? $(LFLAGS) -DIJHPCA
#	mpirun -np 8 ./a.out -n 62500 -c 32 -i 1
