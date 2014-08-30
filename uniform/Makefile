include ../Makefile.include

serial	: main.cxx $(OBJECTS)
	$(CXX) $? $(LFLAGS) -DSerial
	./a.out -n 100000 -c 100 -i 0

parallel: main.cxx $(OBJECTS)
	$(CXX) $? $(LFLAGS)
	mpirun -np 8 ./a.out -n 100000 -c 100 -i 0
