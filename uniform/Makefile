include ../Makefile.include

serial	: main.cxx
	$(CXX) $? $(LFLAGS) -DSerial
	./a.out

parallel: main.cxx
	$(CXX) $? $(LFLAGS)
	mpirun -np 2 ./a.out
