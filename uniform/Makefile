include ../Makefile.include

serial	: main.cxx $(OBJECTS)
	$(CXX) $? $(LFLAGS) -DSerial
	./a.out -n 10000 -c 32 -i 3

parallel: main.cxx $(OBJECTS)
	$(CXX) $? $(LFLAGS)
	mpirun -np 8 ./a.out -n 10000 -c 32 -i 3
